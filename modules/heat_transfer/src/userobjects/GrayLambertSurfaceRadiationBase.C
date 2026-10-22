//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GrayLambertSurfaceRadiationBase.h"
#include "ViewFactorMatrix.h"
#include "SideSetPatch.h"
#include "MathUtils.h"
#include "Function.h"
#include "MooseMesh.h"
#include "libmesh/quadrature.h"
#include "libmesh/petsc_solver_exception.h"
#include "libmesh/utility.h"

#include <petscksp.h>
#include <cmath>

InputParameters
GrayLambertSurfaceRadiationBase::validParams()
{
  InputParameters params = SideUserObject::validParams();
  params.addParam<Real>(
      "stefan_boltzmann_constant",
      5.670367e-8,
      "The Stefan-Boltzmann constant. Default value is in units of [W / m^2 K^4].");
  params.addRequiredCoupledVar("temperature", "The coupled temperature variable.");
  params.addRequiredParam<std::vector<FunctionName>>("emissivity",
                                                     "Emissivities for each boundary.");
  params.addParam<std::vector<BoundaryName>>(
      "fixed_temperature_boundary",
      {},
      "The list of boundary IDs from the mesh with fixed temperatures.");
  params.addParam<std::vector<FunctionName>>(
      "fixed_boundary_temperatures", {}, "The temperatures of the fixed boundary.");
  params.addParam<std::vector<BoundaryName>>(
      "adiabatic_boundary", {}, "The list of boundary IDs from the mesh that are adiabatic.");

  params.addParam<std::vector<unsigned int>>(
      "n_patches",
      {},
      "Number of patches per boundary. If empty, every boundary listed in 'boundary' is a single "
      "surface. If zero, one patch is created for every face of the boundary. The patches then "
      "replace the patches materialized as sidesets and are indexed internally.");

  MultiMooseEnum partitioning(
      "default=-3 metis=-2 parmetis=-1 linear=0 centroid hilbert_sfc morton_sfc");
  partitioning.addValidName("grid");
  params.addParam<MultiMooseEnum>(
      "partitioners",
      partitioning,
      "The partitioner to use for splitting each boundary into patches. Empty (the default) "
      "uses the default partitioner. Only used when 'n_patches' is provided.");

  MultiMooseEnum direction("x y z radial");
  params.addParam<MultiMooseEnum>("centroid_partitioner_directions",
                                  direction,
                                  "The sort direction if using the centroid partitioner.");

  MooseEnum environment("vacuum black_body", "vacuum");
  params.addParam<MooseEnum>(
      "environment",
      environment,
      "The environment surrounding the surfaces. 'vacuum' assumes a black body environment at "
      "zero temperature, so radiation that escapes the participating surfaces is lost. "
      "'black_body' "
      "assumes a black body environment with the temperature given by 'environment_temperature'.");
  params.addParam<FunctionName>(
      "environment_temperature",
      "The temperature of the black body environment. Only used if 'environment' is "
      "'black_body'.");

  params.addClassDescription(
      "This object implements the exchange of heat by radiation between sidesets.");
  return params;
}

GrayLambertSurfaceRadiationBase::GrayLambertSurfaceRadiationBase(const InputParameters & parameters)
  : SideUserObject(parameters),
    _sigma_stefan_boltzmann(getParam<Real>("stefan_boltzmann_constant")),
    _n_sides(0),
    _temperature(coupledValue("temperature")),
    _internal_patches(false),
    _environment(getParam<MooseEnum>("environment")),
    _environment_temperature(isParamSetByUser("environment_temperature")
                                 ? &getFunctionByName(
                                       getParam<FunctionName>("environment_temperature"))
                                 : nullptr),
    _view_factors(nullptr),
    _mat(nullptr),
    _ksp(nullptr),
    _pc(nullptr),
    _rhs(nullptr),
    _radiosity_vec(nullptr)
{
  const std::vector<BoundaryName> boundary_names = getParam<std::vector<BoundaryName>>("boundary");
  for (const auto & name : boundary_names)
    if (name == "ANY_BOUNDARY_ID")
      paramError("boundary", "boundary must be explicitly provided.");

  // the patches are either sidesets of the mesh or internal; when they are internal the number of
  // patches per boundary is given by the action
  std::vector<unsigned int> n_patches = getParam<std::vector<unsigned int>>("n_patches");
  _internal_patches = !n_patches.empty();
  if (_internal_patches && n_patches.size() != boundary_names.size())
    paramError("n_patches", "The number of entries must match the number of boundary entries.");
  if (!_internal_patches)
    n_patches.assign(boundary_names.size(), 1);

  const MultiMooseEnum & partitioners = getParam<MultiMooseEnum>("partitioners");
  const MultiMooseEnum & directions = getParam<MultiMooseEnum>("centroid_partitioner_directions");
  if (_internal_patches)
  {
    if (partitioners.size() != 0 && partitioners.size() != boundary_names.size())
      paramError("partitioners",
                 "The number of entries must match the number of boundary entries.");
    if (directions.size() != 0 && directions.size() != boundary_names.size())
      paramError("centroid_partitioner_directions",
                 "The number of entries must match the number of boundary entries.");
  }

  // split the boundaries into patches and build the boundary/side to surface maps
  unsigned int offset = 0;
  for (const auto j : index_range(boundary_names))
  {
    const BoundaryID id = _mesh.getBoundaryID(boundary_names[j]);
    _boundary_ids.push_back(id);
    _boundary_index[id] = j;
    _side_id_index[id] = offset;
    _surface_offsets.push_back(offset);

    if (_internal_patches)
    {
      const SideSetPatch patches(_mesh.getMesh(),
                                 id,
                                 n_patches[j],
                                 partitioners.size() ? partitioners[j] : "metis",
                                 directions.size() ? directions[j] : "");
      _n_patches_per_boundary.push_back(patches.nPatches());
      for (const auto & [side, patch] : patches.patchIndices())
        _side_surface[side] = offset + patch;
    }
    else
      _n_patches_per_boundary.push_back(1);

    _surface_boundary.insert(_surface_boundary.end(), _n_patches_per_boundary.back(), j);
    offset += _n_patches_per_boundary.back();
  }
  _n_sides = offset;

  // size the per-surface data
  _fixed_side_temperature.assign(_n_sides, nullptr);
  _radiosity.resize(_n_sides);
  _heat_flux_density.resize(_n_sides);
  _side_temperature.resize(_n_sides);
  _side_type.assign(_n_sides, VARIABLE_TEMPERATURE);
  _areas.resize(_n_sides);
  _beta.resize(_n_sides);
  _surface_irradiation.resize(_n_sides);
  _escape_fractions.assign(_n_sides, 0);

  // get emissivity functions, one per surface
  const auto & eps_names = getParam<std::vector<FunctionName>>("emissivity");
  if (eps_names.size() != boundary_names.size())
    paramError("emissivity", "The number of entries must match the number of boundary entries.");

  if (_internal_patches)
    for (const auto j : index_range(eps_names))
    {
      const auto & function = getFunctionByName(eps_names[j]);
      _emissivity.insert(_emissivity.end(), _n_patches_per_boundary[j], &function);
    }
  else
    for (const auto & name : eps_names)
      _emissivity.push_back(&getFunctionByName(name));

  // get the fixed boundaries of the system if any are provided
  if (isParamValid("fixed_temperature_boundary"))
  {
    if (!isParamValid("fixed_boundary_temperatures"))
      paramError("fixed_boundary_temperatures",
                 "fixed_temperature_boundary is provided, but fixed_boundary_temperatures is not.");

    const auto & fixed_names = getParam<std::vector<BoundaryName>>("fixed_temperature_boundary");
    const auto & fixed_functions =
        getParam<std::vector<FunctionName>>("fixed_boundary_temperatures");

    if (fixed_names.size() != fixed_functions.size())
      paramError(
          "fixed_boundary_temperatures",
          "fixed_boundary_temperatures and fixed_temperature_boundary must have the same length.");

    for (const auto k : index_range(fixed_names))
    {
      const auto it = _boundary_index.find(_mesh.getBoundaryID(fixed_names[k]));
      if (it == _boundary_index.end())
        paramError("fixed_temperature_boundary",
                   "fixed_temperature_boundary must be a subset of boundary.");

      const auto [begin, end] = surfaceRange(it->second);
      const Function * const function = &getFunctionByName(fixed_functions[k]);
      for (const auto i : make_range(begin, end))
      {
        _side_type[i] = FIXED_TEMPERATURE;
        _fixed_side_temperature[i] = function;
      }
    }
  }

  // get the adiabatic boundaries of the system if any are provided
  if (isParamValid("adiabatic_boundary"))
  {
    for (const auto & name : getParam<std::vector<BoundaryName>>("adiabatic_boundary"))
    {
      const BoundaryID id = _mesh.getBoundaryID(name);
      const auto it = _boundary_index.find(id);
      if (it == _boundary_index.end())
        paramError("adiabatic_boundary", "adiabatic_boundary must be a subset of boundary.");

      _adiabatic_side_ids.insert(id);

      const auto [begin, end] = surfaceRange(it->second);
      for (const auto i : make_range(begin, end))
      {
        if (_side_type[i] == FIXED_TEMPERATURE)
          paramError("adiabatic_boundary",
                     "Isothermal boundary cannot also be adiabatic boundary.");
        _side_type[i] = ADIABATIC;
      }
    }
  }

  // check that the environment parameters are consistent
  if (_environment == "black_body" && !_environment_temperature)
    paramError("environment_temperature",
               "The 'environment_temperature' function must be provided if 'environment' is "
               "'black_body'.");
}

GrayLambertSurfaceRadiationBase::~GrayLambertSurfaceRadiationBase()
{
  if (_ksp)
    LibmeshPetscCall(KSPDestroy(&_ksp));
  if (_mat)
    LibmeshPetscCall(MatDestroy(&_mat));
  if (_rhs)
    LibmeshPetscCall(VecDestroy(&_rhs));
  if (_radiosity_vec)
    LibmeshPetscCall(VecDestroy(&_radiosity_vec));
}

unsigned int
GrayLambertSurfaceRadiationBase::surfaceIndex(const libMesh::Elem & elem,
                                              const unsigned int side) const
{
  if (_internal_patches)
    return libmesh_map_find(_side_surface, std::make_pair(elem.id(), side));

  // the surfaces are sidesets, so the surface is the boundary the side belongs to
  std::vector<BoundaryID> ids;
  _mesh.getMesh().get_boundary_info().boundary_ids(&elem, side, ids);
  for (const auto id : ids)
  {
    const auto it = _side_id_index.find(id);
    if (it != _side_id_index.end())
      return it->second;
  }

  mooseError("Side ",
             side,
             " of element ",
             elem.id(),
             " does not belong to any of the participating boundaries.");
}

void
GrayLambertSurfaceRadiationBase::execute()
{
  const unsigned int index =
      _internal_patches ? libmesh_map_find(_side_surface,
                                           std::make_pair(_current_elem->id(), _current_side))
                        : getSideIDIndex(_current_boundary_id);

  for (unsigned int qp = 0; qp < _qrule->n_points(); qp++)
  {
    _areas[index] += _JxW[qp] * _coord[qp];

    Real temp = 0;
    if (_side_type[index] == ADIABATIC)
      continue;
    else if (_side_type[index] == VARIABLE_TEMPERATURE)
      temp = _temperature[qp];
    else if (_side_type[index] == FIXED_TEMPERATURE)
    {
      const Function * const function = _fixed_side_temperature[index];
      mooseAssert(function, "Fixed temperature surface without a temperature function.");
      temp = function->value(_t, _q_point[qp]);
    }

    _beta[index] += _JxW[qp] * _coord[qp] * _sigma_stefan_boltzmann *
                    _emissivity[index]->value(temp, Point()) * MathUtils::pow(temp, 4);
    _side_temperature[index] += _JxW[qp] * _coord[qp] * temp;
  }
}

void
GrayLambertSurfaceRadiationBase::initialize()
{
  // view factors are obtained here to make sure that another object had
  // time to compute them on exec initial
  _view_factors = setViewFactors();
  _escape_fractions = setEscapeFractions();

  // create the sparse matrix and the GMRES/ILU(0) solver once, the first time we
  // see the view factors. The pattern is the nonzero structure of the view factor
  // matrix plus the diagonal; the entries are restamped in finalize() at every
  // nonlinear iteration.
  if (!_mat)
  {
    LibmeshPetscCall(MatCreate(PETSC_COMM_SELF, &_mat));
    LibmeshPetscCall(MatSetSizes(_mat, _n_sides, _n_sides, _n_sides, _n_sides));
    LibmeshPetscCall(MatSetUp(_mat));

    for (const auto i : make_range(_n_sides))
    {
      // a zero value still registers the (i, j) location in the pattern
      std::vector<PetscInt> cols;
      std::vector<PetscReal> vals;
      cols.push_back(i);
      vals.push_back(0);
      // the diagonal is registered above, and the pattern must not list a column twice
      for (const auto & [j, value] : _view_factors->row(i))
        if (j != i)
        {
          cols.push_back(j);
          vals.push_back(0);
        }
      const PetscInt row = i;
      LibmeshPetscCall(MatSetValues(
          _mat, 1, &row, (PetscInt)cols.size(), cols.data(), vals.data(), INSERT_VALUES));
    }
    LibmeshPetscCall(MatAssemblyBegin(_mat, MAT_FLUSH_ASSEMBLY));
    LibmeshPetscCall(MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY));

    LibmeshPetscCall(KSPCreate(PETSC_COMM_SELF, &_ksp));
    LibmeshPetscCall(KSPSetType(_ksp, KSPGMRES));
    LibmeshPetscCall(KSPSetOperators(_ksp, _mat, _mat));
    // the PC is owned by the KSP, so it must not be destroyed separately
    LibmeshPetscCall(KSPGetPC(_ksp, &_pc));
    LibmeshPetscCall(PCSetType(_pc, PCILU));
    // 1e-13 is tight enough that the iterative solution matches the previous direct
    // solve within the last digit of the regression test outputs; the rtol test is a
    // relative residual test by default
    LibmeshPetscCall(KSPSetTolerances(_ksp, 1e-13, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT));

    LibmeshPetscCall(VecCreateSeq(PETSC_COMM_SELF, _n_sides, &_rhs));
    LibmeshPetscCall(VecCreateSeq(PETSC_COMM_SELF, _n_sides, &_radiosity_vec));
  }

  // initialize areas, beta, side temps
  for (unsigned int j = 0; j < _n_sides; ++j)
  {
    _areas[j] = 0;
    _beta[j] = 0;
    _side_temperature[j] = 0;
  }
}

void
GrayLambertSurfaceRadiationBase::finalize()
{
  // need to do some parallel communiction here
  gatherSum(_areas);
  gatherSum(_beta);
  gatherSum(_side_temperature);

  // first compute averages from the totals
  for (unsigned int j = 0; j < _n_sides; ++j)
  {
    _beta[j] /= _areas[j];
    _side_temperature[j] /= _areas[j];
  }

  // black body radiation of the environment surrounding the surfaces
  const Real environment_black_body =
      _environment == "black_body"
          ? _sigma_stefan_boltzmann *
                MathUtils::pow(_environment_temperature->value(_t, Point()), 4)
          : 0;

  // stamp the entries of A = I - F o diag(c) at the current emissivities; c is the
  // reflectivity 1 - emissivity, or 1 for adiabatic surfaces. The pattern was
  // created in initialize().
  LibmeshPetscCall(MatZeroEntries(_mat));
  for (const auto i : make_range(_n_sides))
  {
    const Real c =
        _side_type[i] == ADIABATIC ? 1 : 1 - _emissivity[i]->value(_side_temperature[i], Point());
    std::vector<PetscInt> cols;
    std::vector<PetscReal> vals;
    // the diagonal is 1 even if the view factor matrix has no diagonal entry
    cols.push_back(i);
    vals.push_back(1 - c * _view_factors->get(i, i));
    for (const auto & [j, value] : _view_factors->row(i))
      if (j != i)
      {
        cols.push_back(j);
        vals.push_back(-c * value);
      }
    const PetscInt row = i;
    LibmeshPetscCall(MatSetValues(
        _mat, 1, &row, (PetscInt)cols.size(), cols.data(), vals.data(), INSERT_VALUES));
  }
  LibmeshPetscCall(MatAssemblyBegin(_mat, MAT_FLUSH_ASSEMBLY));
  LibmeshPetscCall(MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY));

  // right-hand side of the radiosity system
  LibmeshPetscCall(VecZeroEntries(_rhs));
  for (const auto i : make_range(_n_sides))
  {
    Real rhs_i = _beta[i];
    // radiation that escapes to the environment is re-emitted according to its
    // reflectivity (or fully absorbed and balanced for adiabatic surfaces)
    if (_side_type[i] == ADIABATIC)
      rhs_i += _escape_fractions[i] * environment_black_body;
    else
      rhs_i += (1 - _emissivity[i]->value(_side_temperature[i], Point())) * _escape_fractions[i] *
               environment_black_body;
    LibmeshPetscCall(VecSetValue(_rhs, i, rhs_i, ADD_VALUES));
  }
  LibmeshPetscCall(VecAssemblyBegin(_rhs));
  LibmeshPetscCall(VecAssemblyEnd(_rhs));

  // the ILU(0) factors are stale because the entries changed; PETSc does not
  // detect this, so re-setup the preconditioner explicitly
  LibmeshPetscCall(PCSetUp(_pc));

  // compute the radiosityes
  LibmeshPetscCall(KSPSolve(_ksp, _rhs, _radiosity_vec));
  KSPConvergedReason converged_reason = KSP_CONVERGED_ITERATING;
  LibmeshPetscCallQ(KSPGetConvergedReason(_ksp, &converged_reason));
  if (converged_reason < 0)
    mooseError("The radiosity solve did not converge (PETSc KSP reason: ",
               static_cast<int>(converged_reason),
               ").");

  const PetscScalar * radiosity = nullptr;
  LibmeshPetscCallQ(VecGetArrayRead(_radiosity_vec, &radiosity));

  // store the radiosity, temperatures and heat flux density for each surface
  for (const auto i : make_range(_n_sides))
  {
    _radiosity[i] = radiosity[i];

    // _heat_flux_density is obtained from a somewhat cumbersome relation
    // but it has the advantage that we do not divide by 1 - emissivity
    // which blows up for black bodies
    _heat_flux_density[i] = radiosity[i];
    for (const auto & [j, value] : _view_factors->row(i))
      _heat_flux_density[i] -= value * radiosity[j];
    _heat_flux_density[i] -= _escape_fractions[i] * environment_black_body;

    if (_side_type[i] == ADIABATIC)
    {
      Real e = _emissivity[i]->value(_side_temperature[i], Point());
      _side_temperature[i] = std::pow(
          (radiosity[i] + (1 - e) / e * _heat_flux_density[i]) / _sigma_stefan_boltzmann, 0.25);
    }

    // compute the surface irradiation into i from the radiosities
    _surface_irradiation[i] = 0;
    for (const auto & [j, value] : _view_factors->row(i))
      _surface_irradiation[i] += value * radiosity[j];
    _surface_irradiation[i] += _escape_fractions[i] * environment_black_body;
  }

  LibmeshPetscCallQ(VecRestoreArrayRead(_radiosity_vec, &radiosity));
}

std::vector<Real>
GrayLambertSurfaceRadiationBase::setEscapeFractions()
{
  // by default, no radiation escapes to the environment (closed geometries)
  return std::vector<Real>(_n_sides, 0);
}

void
GrayLambertSurfaceRadiationBase::threadJoin(const UserObject & y)
{
  const auto & pps = cast_ref<const GrayLambertSurfaceRadiationBase &>(y);

  for (unsigned int j = 0; j < _n_sides; ++j)
  {
    _areas[j] += pps._areas[j];
    _side_temperature[j] += pps._side_temperature[j];
    _beta[j] += pps._beta[j];
  }
}

std::set<BoundaryID>
GrayLambertSurfaceRadiationBase::getSurfaceIDs() const
{
  std::set<BoundaryID> surface_ids;
  for (auto & p : _side_id_index)
    surface_ids.insert(p.first);
  return surface_ids;
}

unsigned int
GrayLambertSurfaceRadiationBase::getSideIDIndex(BoundaryID id) const
{
  const auto it = _side_id_index.find(id);
  if (it == _side_id_index.end())
    mooseError("Sideset ID ", id, " (", _mesh.getBoundaryName(id), ") not found in map.");
  return it->second;
}

unsigned int
GrayLambertSurfaceRadiationBase::getBoundaryIndex(const BoundaryID id) const
{
  const auto it = _boundary_index.find(id);
  if (it == _boundary_index.end())
    mooseError("Sideset ID ", id, " (", _mesh.getBoundaryName(id), ") not found in map.");
  return it->second;
}

std::pair<unsigned int, unsigned int>
GrayLambertSurfaceRadiationBase::surfaceRange(const unsigned int j) const
{
  return {_surface_offsets[j], _surface_offsets[j] + _n_patches_per_boundary[j]};
}

Real
GrayLambertSurfaceRadiationBase::getSurfaceIrradiation(const BoundaryID id) const
{
  if (!_internal_patches)
    return _surface_irradiation[getSideIDIndex(id)];

  const auto [begin, end] = surfaceRange(getBoundaryIndex(id));
  Real value = 0;
  Real area = 0;
  for (const auto i : make_range(begin, end))
  {
    value += _surface_irradiation[i] * _areas[i];
    area += _areas[i];
  }
  return value / area;
}

Real
GrayLambertSurfaceRadiationBase::getSurfaceHeatFluxDensity(const BoundaryID id) const
{
  if (!_internal_patches)
    return _heat_flux_density[getSideIDIndex(id)];

  const auto [begin, end] = surfaceRange(getBoundaryIndex(id));
  Real value = 0;
  Real area = 0;
  for (const auto i : make_range(begin, end))
  {
    value += _heat_flux_density[i] * _areas[i];
    area += _areas[i];
  }
  return value / area;
}

Real
GrayLambertSurfaceRadiationBase::getSurfaceTemperature(const BoundaryID id) const
{
  if (!_internal_patches)
    return _side_temperature[getSideIDIndex(id)];

  const auto [begin, end] = surfaceRange(getBoundaryIndex(id));
  Real value = 0;
  Real area = 0;
  for (const auto i : make_range(begin, end))
  {
    value += _side_temperature[i] * _areas[i];
    area += _areas[i];
  }
  return value / area;
}

Real
GrayLambertSurfaceRadiationBase::getSurfaceRadiosity(const BoundaryID id) const
{
  if (!_internal_patches)
    return _radiosity[getSideIDIndex(id)];

  const auto [begin, end] = surfaceRange(getBoundaryIndex(id));
  Real value = 0;
  Real area = 0;
  for (const auto i : make_range(begin, end))
  {
    value += _radiosity[i] * _areas[i];
    area += _areas[i];
  }
  return value / area;
}

Real
GrayLambertSurfaceRadiationBase::getSurfaceEmissivity(const BoundaryID id) const
{
  if (!_internal_patches)
  {
    const auto ind = getSideIDIndex(id);
    return _emissivity[ind]->value(_side_temperature[ind], Point());
  }

  const auto [begin, end] = surfaceRange(getBoundaryIndex(id));
  Real value = 0;
  Real area = 0;
  for (const auto i : make_range(begin, end))
  {
    value += _emissivity[i]->value(_side_temperature[i], Point()) * _areas[i];
    area += _areas[i];
  }
  return value / area;
}

Real
GrayLambertSurfaceRadiationBase::getViewFactor(const BoundaryID from_id,
                                              const BoundaryID to_id) const
{
  if (!_internal_patches)
    return _view_factors->get(getSideIDIndex(from_id), getSideIDIndex(to_id));

  // the view factor between two boundaries is the area-weighted average of the view factors
  // between their patches
  const auto [from_begin, from_end] = surfaceRange(getBoundaryIndex(from_id));
  const auto [to_begin, to_end] = surfaceRange(getBoundaryIndex(to_id));
  Real value = 0;
  Real area = 0;
  for (const auto i : make_range(from_begin, from_end))
  {
    area += _areas[i];
    for (const auto j : make_range(to_begin, to_end))
      value += _view_factors->get(i, j) * _areas[i];
  }
  return value / area;
}

Real
GrayLambertSurfaceRadiationBase::getSurfaceIrradiation(const libMesh::Elem & elem,
                                                       const unsigned int side) const
{
  return _surface_irradiation[surfaceIndex(elem, side)];
}

Real
GrayLambertSurfaceRadiationBase::getSurfaceHeatFluxDensity(const libMesh::Elem & elem,
                                                           const unsigned int side) const
{
  return _heat_flux_density[surfaceIndex(elem, side)];
}

Real
GrayLambertSurfaceRadiationBase::getSurfaceEmissivity(const libMesh::Elem & elem,
                                                      const unsigned int side) const
{
  const auto index = surfaceIndex(elem, side);
  return _emissivity[index]->value(_side_temperature[index], Point());
}
