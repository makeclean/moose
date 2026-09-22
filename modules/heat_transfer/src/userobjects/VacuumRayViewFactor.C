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
  : ViewFactorBase(parameters),
    _study(getUserObject<VacuumRayViewFactorStudy>("ray_study_name")),
    _escape_fractions(_n_sides, 0)
{
  if (getParam<bool>("normalize_view_factor"))
    paramError("normalize_view_factor",
               "Normalization of view factors is not supported with this object, because the "
               "radiation that escapes to the environment must not be redistributed among the "
               "participating surfaces.");
  if (_mesh.dimension() != 3)
    mooseError("The vacuum view factor calculation only supports three-dimensional meshes.");
}

void
VacuumRayViewFactor::threadJoinViewFactor(const UserObject & y)
{
  const auto & vf = cast_ref<const VacuumRayViewFactor &>(y);
  for (unsigned int i = 0; i < _n_sides; ++i)
  {
    _escape_fractions[i] += vf._escape_fractions[i];
    for (unsigned int j = 0; j < _n_sides; ++j)
      _view_factors[i][j] += vf._view_factors[i][j];
  }
}

void
VacuumRayViewFactor::finalizeViewFactor()
{
  // the ray study has already summed the ray weights over the processors and
  // threads; the weights are normalized by the area of the leaving surface times
  // pi, the integral of the angular weight function over the hemisphere
  for (const auto i : make_range(_n_sides))
  {
    const BoundaryID from_id = boundaryIDs()[i];
    const Real normalization = 1. / (_areas[i] * libMesh::pi);
    _escape_fractions[i] = _study.escapeWeight(from_id) * normalization;
    for (const auto j : make_range(_n_sides))
      _view_factors[i][j] = _study.viewFactorWeight(from_id, boundaryIDs()[j]) * normalization;
  }
}

Real
VacuumRayViewFactor::getEscapeFraction(BoundaryID id) const
{
  return _escape_fractions[getSideNameIndex(_mesh.getBoundaryName(id))];
}

#endif // MOOSE_XDG_ENABLED