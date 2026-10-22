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
#include "MathUtils.h"
#include "Function.h"
#include "libmesh/quadrature.h"
#include "libmesh/petsc_solver_exception.h"

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
    _n_sides(boundaryIDs().size()),
    _temperature(coupledValue("temperature")),
    _radiosity(_n_sides),
    _heat_flux_density(_n_sides),
    _side_temperature(_n_sides),
    _side_type(_n_sides),
    _areas(_n_sides),
    _beta(_n_sides),
    _surface_irradiation(_n_sides),
    _escape_fractions(_n_sides, 0),
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
  // get emissivity functions
  auto & eps_names = getParam<std::vector<FunctionName>>("emissivity");
  _emissivity.resize(eps_names.size());
  for (unsigned int j = 0; j < _emissivity.size(); ++j)
    _emissivity[j] = &getFunctionByName(eps_names[j]);

  // set up the map from the side id to the local index & check
  // note that boundaryIDs is not in the right order anymore!
  {
    std::vector<BoundaryName> boundary_names = getParam<std::vector<BoundaryName>>("boundary");

    for (unsigned int j = 0; j < boundary_names.size(); ++j)
    {
      if (boundary_names[j] == "ANY_BOUNDARY_ID")
        paramError("boundary", "boundary must be explicitly provided.");

      _side_id_index[_mesh.getBoundaryID(boundary_names[j])] = j;
      _side_type[j] = VARIABLE_TEMPERATURE;
    }

    // consistency check on emissivity, must be as many entries as boundary
    if (boundary_names.size() != _emissivity.size())
      paramError("emissivity", "The number of entries must match the number of boundary entries.");
  }

  // get the fixed boundaries of the system if any are provided
  if (isParamValid("fixed_temperature_boundary"))
  {
    // if fixed_temperature_boundary is valid, then fixed_side_temperatures must be
    // valid, too
    if (!isParamValid("fixed_boundary_temperatures"))
      paramError("fixed_boundary_temperatures",
                 "fixed_temperature_boundary is provided, but fixed_boundary_temperatures is not.");

    auto fst_fn = getParam<std::vector<FunctionName>>("fixed_boundary_temperatures");
    for (auto & fn : fst_fn)
      _fixed_side_temperature.push_back(&getFunctionByName(fn));

    // get the fixed boundaries and temperatures
    std::vector<BoundaryName> boundary_names =
        getParam<std::vector<BoundaryName>>("fixed_temperature_boundary");

    if (boundary_names.size() != _fixed_side_temperature.size())
      paramError(
          "fixed_boundary_temperatures",
          "fixed_boundary_temperatures and fixed_temperature_boundary must have the same length.");

    unsigned int index = 0;
    for (auto & name : boundary_names)
    {
      _fixed_side_id_index[_mesh.getBoundaryID(name)] = index;
      index++;
    }

    // check that fixed side ids is a subset of boundary ids
    // and update _side_type info
    for (auto & p : _fixed_side_id_index)
    {
      if (_side_id_index.find(p.first) == _side_id_index.end())
        paramError("fixed_temperature_boundary",
                   "fixed_temperature_boundary must be a subset of boundary.");
      _side_type[_side_id_index.find(p.first)->second] = FIXED_TEMPERATURE;
    }
  }

  // get the fixed boundaries of the system if any are provided
  if (isParamValid("adiabatic_boundary"))
  {
    // get the adiabatic boundaries and temperatures
    std::vector<BoundaryName> boundary_names =
        getParam<std::vector<BoundaryName>>("adiabatic_boundary");

    for (auto & name : boundary_names)
      _adiabatic_side_ids.insert(_mesh.getBoundaryID(name));

    // check that adiabatic side ids is a subset of boundary ids
    // and update _side_type info
    for (auto & id : _adiabatic_side_ids)
    {
      if (_side_id_index.find(id) == _side_id_index.end())
        paramError("adiabatic_boundary", "adiabatic_boundary must be a subset of boundary.");
      _side_type[_side_id_index.find(id)->second] = ADIABATIC;
    }

    // make sure that adiabatic boundaries are not already fixed boundaries
    for (auto & id : _adiabatic_side_ids)
      if (_fixed_side_id_index.find(id) != _fixed_side_id_index.end())
        paramError("adiabatic_boundary", "Isothermal boundary cannot also be adiabatic boundary.");
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

void
GrayLambertSurfaceRadiationBase::execute()
{
  mooseAssert(_side_id_index.find(_current_boundary_id) != _side_id_index.end(),
              "Current boundary id not in _side_id_index.");
  unsigned int index = _side_id_index.find(_current_boundary_id)->second;

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
      unsigned int iso_index = _fixed_side_id_index.find(_current_boundary_id)->second;
      temp = _fixed_side_temperature[iso_index]->value(_t, _q_point[qp]);
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

Real
GrayLambertSurfaceRadiationBase::getSurfaceIrradiation(BoundaryID id) const
{
  return _surface_irradiation[getSideIDIndex(id)];
}

Real
GrayLambertSurfaceRadiationBase::getSurfaceHeatFluxDensity(BoundaryID id) const
{
  return _heat_flux_density[getSideIDIndex(id)];
}

Real
GrayLambertSurfaceRadiationBase::getSurfaceTemperature(BoundaryID id) const
{
  return _side_temperature[getSideIDIndex(id)];
}

Real
GrayLambertSurfaceRadiationBase::getSurfaceRadiosity(BoundaryID id) const
{
  return _radiosity[getSideIDIndex(id)];
}

Real
GrayLambertSurfaceRadiationBase::getSurfaceEmissivity(BoundaryID id) const
{
  const auto ind = getSideIDIndex(id);
  return _emissivity[ind]->value(_side_temperature[ind], Point());
}

Real
GrayLambertSurfaceRadiationBase::getViewFactor(BoundaryID from_id, BoundaryID to_id) const
{
  return _view_factors->get(getSideIDIndex(from_id), getSideIDIndex(to_id));
}
