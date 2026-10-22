//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ViewFactorObjectSurfaceRadiation.h"
#include "ViewFactorBase.h"
#include "MooseUtils.h"

#include <map>

registerMooseObject("HeatTransferApp", ViewFactorObjectSurfaceRadiation);

InputParameters
ViewFactorObjectSurfaceRadiation::validParams()
{
  InputParameters params = GrayLambertSurfaceRadiationBase::validParams();
  params.addRequiredParam<UserObjectName>("view_factor_object_name",
                                          "Name of the ViewFactor userobjects.");
  params.addClassDescription(
      "ViewFactorObjectSurfaceRadiation computes radiative heat transfer between side sets and the "
      "view factors are computed by a ViewFactor object");
  return params;
}

ViewFactorObjectSurfaceRadiation::ViewFactorObjectSurfaceRadiation(
    const InputParameters & parameters)
  : GrayLambertSurfaceRadiationBase(parameters)
{
}

const ViewFactorMatrix *
ViewFactorObjectSurfaceRadiation::setViewFactors()
{
  const ViewFactorBase & view_factor_uo = getUserObject<ViewFactorBase>("view_factor_object_name");

  // when the patches are internal, the surfaces of the view factor object and of this object are
  // ordered the same way (boundary by boundary, patch by patch), so no permutation is needed
  if (view_factor_uo.usesInternalPatches())
    return &view_factor_uo.viewFactorMatrix();

  // this object and the view factor object may list their boundaries in different orders, so
  // the view factor entries are permuted into this object's ordering
  const auto & vf_names = view_factor_uo.boundaryNames();
  const auto & my_names = boundaryNames();

  std::map<BoundaryName, unsigned int> vf_index;
  for (const auto j : index_range(vf_names))
    vf_index[vf_names[j]] = j;

  std::map<BoundaryName, unsigned int> my_index;
  for (const auto i : index_range(my_names))
    my_index[my_names[i]] = i;

  const auto & vf = view_factor_uo.viewFactorMatrix();
  _vf_matrix = ViewFactorMatrix(_n_sides);
  for (const auto i : index_range(my_names))
    for (const auto & [j, value] : vf.row(libmesh_map_find(vf_index, my_names[i])))
      _vf_matrix.add(i, libmesh_map_find(my_index, vf_names[j]), value);

  return &_vf_matrix;
}

std::vector<Real>
ViewFactorObjectSurfaceRadiation::setEscapeFractions()
{
  const ViewFactorBase & view_factor_uo = getUserObject<ViewFactorBase>("view_factor_object_name");
  std::vector<Real> esc(_n_sides, 0);

  if (view_factor_uo.usesInternalPatches())
    for (const auto i : make_range(_n_sides))
      esc[i] = view_factor_uo.getEscapeFraction(i);
  else
  {
    std::vector<BoundaryName> boundary_names = getParam<std::vector<BoundaryName>>("boundary");
    for (unsigned int i = 0; i < _n_sides; ++i)
      esc[i] = view_factor_uo.getEscapeFraction(_mesh.getBoundaryID(boundary_names[i]));
  }

  return esc;
}
