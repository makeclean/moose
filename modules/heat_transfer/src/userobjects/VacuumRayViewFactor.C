//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VacuumRayViewFactor.h"

#ifdef MOOSE_XDG_ENABLED
#include "VacuumRayViewFactorStudy.h"

registerMooseObject("HeatTransferApp", VacuumRayViewFactor);

InputParameters
VacuumRayViewFactor::validParams()
{
  InputParameters params = ViewFactorBase::validParams();
  params.addRequiredParam<UserObjectName>("ray_study_name",
                                          "Name of the view factor ray study UO.");
  // the escape fraction accounts for the radiation leaving the open geometry;
  // normalization would redistribute it among the participating surfaces
  params.set<bool>("normalize_view_factor") = false;
  params.addClassDescription(
      "Computes view factors and escape fractions for surfaces radiating into an open "
      "(vacuum) environment using ray tracing.");
  return params;
}

VacuumRayViewFactor::VacuumRayViewFactor(const InputParameters & parameters)
  : ViewFactorBase(parameters), _study(getUserObject<VacuumRayViewFactorStudy>("ray_study_name"))
{
  if (getParam<bool>("normalize_view_factor"))
    paramError("normalize_view_factor",
               "Normalization of view factors is not supported with this object, because the "
               "radiation that escapes to the environment must not be redistributed among the "
               "participating surfaces.");
  if (_mesh.dimension() != 3)
    mooseError("The vacuum view factor calculation only supports three-dimensional meshes.");

  // the surfaces are the patches the study splits the participating boundaries into
  _n_sides = _study.nSurfaces();
  _areas.assign(_n_sides, 0);
  _view_factors = ViewFactorMatrix(_n_sides);
  _escape_fractions.assign(_n_sides, 0);
}

void
VacuumRayViewFactor::threadJoinViewFactor(const UserObject & y)
{
  const auto & vf = cast_ref<const VacuumRayViewFactor &>(y);
  _view_factors.merge(vf._view_factors);
  for (const auto i : make_range(_n_sides))
    _escape_fractions[i] += vf._escape_fractions[i];
}

void
VacuumRayViewFactor::finalizeViewFactor()
{
  // the ray study has already summed the ray weights over the processors and
  // threads; the weights are normalized by the area of the leaving surface times
  // pi, the integral of the angular weight function over the hemisphere
  _view_factors.clear();
  for (const auto i : make_range(_n_sides))
  {
    const Real normalization = 1. / (_areas[i] * libMesh::pi);
    _escape_fractions[i] = _study.escapeWeight(i) * normalization;
    for (const auto & [j, weight] : _study.viewFactorRow(i))
      if (weight != 0.0)
        _view_factors.add(i, j, weight * normalization);
  }
}

Real
VacuumRayViewFactor::getEscapeFraction(BoundaryID id) const
{
  // the escape fractions of an object with internal patches are indexed by patch, not by
  // boundary; returning the entry of the first patch of the boundary would be misleading
  if (_n_sides != boundaryNames().size())
    mooseError("The escape fractions of this object are indexed by patch, not by boundary.");

  return _escape_fractions[getSideNameIndex(_mesh.getBoundaryName(id))];
}

Real
VacuumRayViewFactor::getEscapeFraction(const unsigned int surface) const
{
  return _escape_fractions[surface];
}

unsigned int
VacuumRayViewFactor::surfaceIndex(const libMesh::Elem & elem, const unsigned int side) const
{
  return _study.surfaceIndex(elem, side);
}

std::string
VacuumRayViewFactor::surfaceName(const unsigned int surface) const
{
  return _study.surfaceName(surface);
}

BoundaryID
VacuumRayViewFactor::surfaceBoundaryID(const unsigned int surface) const
{
  return _study.surfaceBoundaryID(surface);
}

#endif // MOOSE_XDG_ENABLED
