//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ViewFactorBase.h"
#include "libmesh/quadrature.h"

#include <limits>
#include <set>

InputParameters
ViewFactorBase::validParams()
{
  InputParameters params = SideUserObject::validParams();
  params.addParam<Real>("view_factor_tol",
                        std::numeric_limits<Real>::max(),
                        "Tolerance for checking view factors. Default is to allow everything.");
  params.addParam<bool>(
      "print_view_factor_info", false, "Flag to print information about computed view factors.");
  params.addParam<bool>("normalize_view_factor",
                        true,
                        "Determines if view factors are normalized to sum to one (consistent with "
                        "their definition).");

  MooseEnum normalization_method("lagrange_multiplier inverse_row_sum", "lagrange_multiplier");
  normalization_method.addDocumentation("lagrange_multiplier", "Uses Lagrange multiplier method");
  normalization_method.addDocumentation("inverse_row_sum",
                                        "Divides each view factor by its row sum");
  params.addParam<MooseEnum>("normalization_method",
                             normalization_method,
                             "Normalization method to use if 'normalize_view_factor' is 'true'.");

  params.addClassDescription(
      "A base class for automatic computation of view factors between sidesets.");
  return params;
}

ViewFactorBase::ViewFactorBase(const InputParameters & parameters)
  : SideUserObject(parameters),
    _n_sides(boundaryIDs().size()),
    _areas(_n_sides),
    _view_factor_tol(getParam<Real>("view_factor_tol")),
    _normalize_view_factor(getParam<bool>("normalize_view_factor")),
    _normalization_method(
        getParam<MooseEnum>("normalization_method").getEnum<NormalizationMethod>()),
    _print_view_factor_info(getParam<bool>("print_view_factor_info")),
    _view_factors(_n_sides)
{
  // set up the map from the side id to the local index & side name to local index
  std::vector<BoundaryName> boundary_names = getParam<std::vector<BoundaryName>>("boundary");
  for (unsigned int j = 0; j < boundary_names.size(); ++j)
    _side_name_index[boundary_names[j]] = j;
}

unsigned int
ViewFactorBase::getSideNameIndex(std::string name) const
{
  auto it = _side_name_index.find(name);
  if (it == _side_name_index.end())
    mooseError("Boundary ", name, " does not exist.");
  return it->second;
}

Real
ViewFactorBase::getViewFactor(BoundaryID from_id, BoundaryID to_id) const
{
  auto from_name = _mesh.getBoundaryName(from_id);
  auto to_name = _mesh.getBoundaryName(to_id);

  return getViewFactor(from_name, to_name);
}

Real
ViewFactorBase::getViewFactor(BoundaryName from_name, BoundaryName to_name) const
{
  // the view factors of an object with internal patches are indexed by patch, not by boundary;
  // returning the entry of the first patch of each boundary would be misleading
  if (_n_sides != boundaryNames().size())
    mooseError("The view factors of this object are indexed by patch, not by boundary. Use "
               "GrayLambertSurfaceRadiationBase::getViewFactor() to retrieve the area-weighted "
               "view factor between two boundaries.");

  auto from = _side_name_index.find(from_name);
  auto to = _side_name_index.find(to_name);
  if (from == _side_name_index.end())
    mooseError("Boundary id ",
               _mesh.getBoundaryID(from_name),
               " with name ",
               from_name,
               " not listed in boundary parameter.");

  if (to == _side_name_index.end())
    mooseError("Boundary id ",
               _mesh.getBoundaryID(to_name),
               " with name ",
               to_name,
               " not listed in boundary parameter.");

  return _view_factors.get(from->second, to->second);
}

const ViewFactorMatrix &
ViewFactorBase::viewFactorMatrix() const
{
  return _view_factors;
}

void
ViewFactorBase::initialize()
{
  std::fill(_areas.begin(), _areas.end(), 0);
}

void
ViewFactorBase::execute()
{
  // compute areas; internal patches are looked up by the side, because all patches of a
  // boundary share the same boundary id
  unsigned int index;
  if (usesInternalPatches())
    index = surfaceIndex(*_current_elem, _current_side);
  else
    index = getSideNameIndex(_mesh.getBoundaryName(_current_boundary_id));

  _areas[index] += _current_side_volume;
}

unsigned int
ViewFactorBase::surfaceIndex(const libMesh::Elem & elem, const unsigned int side) const
{
  // the default surface index is the index of the boundary the side belongs to
  std::vector<boundary_id_type> ids;
  _mesh.getMesh().get_boundary_info().boundary_ids(&elem, side, ids);
  for (const auto id : ids)
  {
    const auto it = _side_name_index.find(_mesh.getBoundaryName(id));
    if (it != _side_name_index.end())
      return it->second;
  }

  mooseError("Side ",
             side,
             " of element ",
             elem.id(),
             " does not belong to any of the boundaries of this object.");
}

std::string
ViewFactorBase::surfaceName(const unsigned int surface) const
{
  mooseAssert(surface < _n_sides, "Surface index out of range.");
  return boundaryNames()[surface];
}

BoundaryID
ViewFactorBase::surfaceBoundaryID(const unsigned int surface) const
{
  return _mesh.getBoundaryID(surfaceName(surface));
}

void
ViewFactorBase::finalize()
{
  gatherSum(_areas);

  finalizeViewFactor();
  checkAndNormalizeViewFactor();
}

void
ViewFactorBase::threadJoin(const UserObject & y)
{
  const auto & pps = cast_ref<const ViewFactorBase &>(y);
  for (unsigned int i = 0; i < _n_sides; ++i)
    _areas[i] += pps._areas[i];

  threadJoinViewFactor(y);
}

Real
ViewFactorBase::devReciprocity(unsigned int i, unsigned int j) const
{
  return _view_factors.get(i, j) - _areas[j] / _areas[i] * _view_factors.get(j, i);
}

Real
ViewFactorBase::maxDevReciprocity() const
{
  // only pairs of which at least one entry is stored can have a nonzero
  // deviation from reciprocity
  std::set<std::pair<unsigned int, unsigned int>> pairs;
  for (const auto i : make_range(_n_sides))
    for (const auto & [j, value] : _view_factors.row(i))
      if (j != i)
        pairs.insert({std::min(i, j), std::max(i, j)});

  Real v = 0;
  for (const auto & [i, j] : pairs)
  {
    Real s = std::abs(devReciprocity(i, j));
    if (s > v)
      v = s;
  }
  return v;
}

Real
ViewFactorBase::viewFactorRowSum(unsigned int i) const
{
  Real s = 0;
  for (const auto & [to, value] : _view_factors.row(i))
    s += value;
  return s;
}

Real
ViewFactorBase::maxDevRowSum() const
{
  Real v = 0;
  for (unsigned int i = 0; i < _n_sides; ++i)
  {
    // the sum of the view factors and the escape fraction must equal one
    Real s = std::abs(1 - viewFactorRowSum(i) - getEscapeFraction(i));
    if (s > v)
      v = s;
  }
  return v;
}

void
ViewFactorBase::checkAndNormalizeViewFactor()
{
  // collect statistics
  Real max_sum_deviation = maxDevRowSum();
  Real max_reciprocity_deviation = maxDevReciprocity();
  Real max_correction = std::numeric_limits<Real>::lowest();
  Real min_correction = std::numeric_limits<Real>::max();

  if (_print_view_factor_info)
    _console << "\nSum of all view factors from side i to all other faces before normalization.\n"
             << std::endl;

  // check view factors
  if (_print_view_factor_info)
    for (unsigned int from = 0; from < _n_sides; ++from)
    {
      auto from_name = surfaceName(from);
      _console << "View factors from surface " << from_name << " sum to "
               << viewFactorRowSum(from) + getEscapeFraction(from) << " ("
               << viewFactorRowSum(from) << " of which is a view factor, and "
               << getEscapeFraction(from) << " is the escape fraction)." << std::endl;
    }

  if (max_sum_deviation > _view_factor_tol)
    mooseError("Maximum deviation of view factor row sum is ",
               max_sum_deviation,
               " exceeding the set tolerance of ",
               _view_factor_tol);

  // normalize view factors
  if (_normalize_view_factor)
  {
    if (_print_view_factor_info)
      _console << "\nNormalizing view factors.\n" << std::endl;

    const auto view_factors_original = _view_factors;

    switch (_normalization_method)
    {
      case NormalizationMethod::LAGRANGE_MULTIPLIER:
        normalizeUsingLagrangeMultiplier();
        break;
      case NormalizationMethod::INVERSE_ROW_SUM:
        normalizeUsingInverseRowSum();
        break;
      default:
        mooseError("Invalid normalization method.");
    }

    // compute max and min correction for reporting
    for (const auto i : make_range(_n_sides))
      for (const auto j : make_range(_n_sides))
      {
        const Real correction = _view_factors.get(i, j) - view_factors_original.get(i, j);
        min_correction = std::min(correction, min_correction);
        max_correction = std::max(correction, max_correction);
      }
  }

  if (_print_view_factor_info)
  {
    _console << "\nFinal view factors.\n" << std::endl;
    for (const auto from : make_range(_n_sides))
    {
      const auto from_id = surfaceBoundaryID(from);
      // only the stored (nonzero) entries are printed
      for (const auto & [to, value] : _view_factors.row(from))
      {
        const auto to_id = surfaceBoundaryID(to);
        _console << surfaceName(from) << " (" << from_id << ") -> " << surfaceName(to) << " ("
                 << to_id << ") = " << value << std::endl;
      }
    }
  }

  // check sum of view factors and maximum deviation on reciprocity
  Real max_sum_deviation_after_normalization = maxDevRowSum();
  Real max_reciprocity_deviation_after_normalization = maxDevReciprocity();

  // print summary
  _console << std::endl;
  _console << COLOR_CYAN << "Summary of the view factor computation"
           << "\n"
           << COLOR_DEFAULT << std::endl;
  if (_normalize_view_factor)
    _console << "Normalization is performed." << std::endl;
  else
    _console << "Normalization is skipped." << std::endl;
  _console << std::setw(60) << std::left << "Number of patches: " << _n_sides << std::endl;
  _console << std::setw(60) << std::left
           << "Max deviation sum = 1 before normalization: " << max_sum_deviation << std::endl;
  _console << std::setw(60) << std::left
           << "Max deviation from reciprocity before normalization: " << max_reciprocity_deviation
           << std::endl;
  if (_normalize_view_factor)
  {
    _console << std::setw(60) << std::left << "Maximum correction: " << max_correction << std::endl;
    _console << std::setw(60) << std::left << "Minimum correction: " << min_correction << std::endl;
    _console << std::setw(60) << "Max deviation sum = 1 after normalization: "
             << max_sum_deviation_after_normalization << std::endl;
    _console << std::setw(60) << std::left << "Max deviation from reciprocity after normalization: "
             << max_reciprocity_deviation_after_normalization << std::endl;
    _console << std::endl;
  }
}

void
ViewFactorBase::normalizeUsingLagrangeMultiplier()
{
  // allocate space
  DenseVector<Real> rhs(_n_sides);
  DenseVector<Real> lmults(_n_sides);
  DenseMatrix<Real> matrix(_n_sides, _n_sides);

  // equations for the Lagrange multiplier
  for (unsigned int i = 0; i < _n_sides; ++i)
  {
    // set the right hand side
    rhs(i) = 1 - viewFactorRowSum(i);

    // contribution from the delta_ii element
    matrix(i, i) = -0.5;

    // contributions from the upper diagonal
    for (unsigned int j = i + 1; j < _n_sides; ++j)
    {
      Real ar = _areas[i] / _areas[j];
      Real f = 2 * (1 + ar * ar);
      matrix(i, i) += -1 / f;
      matrix(i, j) += -1 / f * ar;
      rhs(i) += ar * ar / (1 + ar * ar) * devReciprocity(i, j);
    }

    // contributions from the lower diagonal
    for (unsigned int j = 0; j < i; ++j)
    {
      Real ar = _areas[j] / _areas[i];
      Real f = 2 * (1 + ar * ar);
      matrix(i, i) += -1 / f * ar * ar;
      matrix(i, j) += -1 / f * ar;
      rhs(i) -= ar * devReciprocity(j, i) * (1 - ar * ar / (1 + ar * ar));
    }
  }

  // solve the linear system
  matrix.lu_solve(rhs, lmults);

  // perform corrections but store view factors to temporary array
  // because we cannot modify _view_factors as it's used in this calc
  // the temporary is dense because the correction touches every entry
  std::vector<std::vector<Real>> vf_temp(_n_sides, std::vector<Real>(_n_sides, 0));
  for (const auto i : make_range(_n_sides))
    for (const auto & [j, value] : _view_factors.row(i))
      vf_temp[i][j] = value;
  for (unsigned int i = 0; i < _n_sides; ++i)
    for (unsigned int j = 0; j < _n_sides; ++j)
    {
      Real correction;
      if (i == j)
        correction = -0.5 * lmults(i);
      else
      {
        Real ar = _areas[i] / _areas[j];
        Real f = 2 * (1 + ar * ar);
        correction = -(lmults(i) + lmults(j) * ar + 2 * ar * ar * devReciprocity(i, j)) / f;
      }

      // update the temporary view factor
      vf_temp[i][j] += correction;
    }

  _view_factors.clear();
  for (const auto i : make_range(_n_sides))
    for (const auto j : make_range(_n_sides))
      if (vf_temp[i][j] != 0)
        _view_factors.add(i, j, vf_temp[i][j]);
}

void
ViewFactorBase::normalizeUsingInverseRowSum()
{
  for (unsigned int i = 0; i < _n_sides; ++i)
  {
    const Real row_sum = viewFactorRowSum(i);
    for (auto & entry : _view_factors.row(i))
      entry.second /= row_sum;
  }
}

unsigned int
ViewFactorBase::indexHelper(unsigned int i, unsigned int j) const
{
  mooseAssert(i <= j, "indexHelper requires i <= j");
  if (i == j)
    return _n_sides + i;
  unsigned int pos = 2 * _n_sides;
  for (unsigned int l = 0; l < _n_sides; ++l)
    for (unsigned int m = l + 1; m < _n_sides; ++m)
    {
      if (l == i && m == j)
        return pos;
      else
        ++pos;
    }
  mooseError("Should never get here");
  return 0;
}
