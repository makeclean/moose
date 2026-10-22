//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RayTracingViewFactor.h"

#include "ViewFactorRayStudy.h"

registerMooseObject("HeatTransferApp", RayTracingViewFactor);

InputParameters
RayTracingViewFactor::validParams()
{
  InputParameters params = ViewFactorBase::validParams();
  params.addRequiredParam<UserObjectName>("ray_study_name",
                                          "Name of the view factor ray study UO.");
  params.addClassDescription("Computes view factors for arbitrary geometries using raytracing.");
  return params;
}

RayTracingViewFactor::RayTracingViewFactor(const InputParameters & parameters)
  : ViewFactorBase(parameters), _ray_study(getUserObject<ViewFactorRayStudy>("ray_study_name"))
{
  if (_mesh.dimension() == 1)
    mooseError("View factor calculations do not support 1D");
}

void
RayTracingViewFactor::threadJoinViewFactor(const UserObject & y)
{
  const auto & vf = cast_ref<const RayTracingViewFactor &>(y);
  _view_factors.merge(vf._view_factors);
}

void
RayTracingViewFactor::finalizeViewFactor()
{
  // the view factors are replaced by the values of the ray study below, which has
  // already summed its contributions over the processors and threads
  _view_factors.clear();

  // get the _view_factors from ray study
  for (const auto & from_name : boundaryNames())
    for (const auto & to_name : boundaryNames())
    {
      const auto from_index = getSideNameIndex(from_name);
      const auto to_index = getSideNameIndex(to_name);
      const BoundaryID from_id = _mesh.getBoundaryID(from_name);
      const BoundaryID to_id = _mesh.getBoundaryID(to_name);
      const Real value = _ray_study.viewFactorInfo(from_id, to_id);
      if (value != 0.0)
        _view_factors.add(from_index, to_index, value);
    }

  // divide view_factor Fij by Ai and pi
  const Real divisor = _mesh.dimension() == 2 ? 2 : libMesh::pi;
  for (const auto i : make_range(_n_sides))
  {
    const Real factor = 1. / (_areas[i] * divisor);
    for (auto & [to, value] : _view_factors.row(i))
      value *= factor;
  }
}
