//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ConstantViewFactorSurfaceRadiation.h"

registerMooseObjectDeprecated("HeatTransferApp",
                              ConstantViewFactorSurfaceRadiation,
                              "08/30/2026 24:00");

InputParameters
ConstantViewFactorSurfaceRadiation::validParams()
{
  InputParameters params = GrayLambertSurfaceRadiationBase::validParams();
  params.addRequiredParam<std::vector<std::vector<Real>>>(
      "view_factors", "The view factors from sideset i to sideset j.");
  params.addClassDescription(
      "ConstantViewFactorSurfaceRadiation computes radiative heat transfer between side sets and "
      "the view factors are provided in the input file");
  return params;
}

ConstantViewFactorSurfaceRadiation::ConstantViewFactorSurfaceRadiation(
    const InputParameters & parameters)
  : GrayLambertSurfaceRadiationBase(parameters), _vf_matrix(_n_sides)
{
  mooseDeprecated("ConstantViewFactorSurfaceRadiation is deprecated. Please use "
                  "ViewFactorObjectSurfaceRadiation in conjunction with SpecifiedViewFactor to "
                  "achieve the same effect as this object.");

  const auto & vf = getParam<std::vector<std::vector<Real>>>("view_factors");

  // check that the input has the right format
  if (vf.size() != _n_sides)
    paramError("view_factors",
               "Leading dimension of view_factors must be equal to number of side sets.");

  for (const auto i : index_range(vf))
    if (vf[i].size() != _n_sides)
      paramError("view_factors",
                 "view_factors must be provided as square array. Row ",
                 i,
                 " has ",
                 vf[i].size(),
                 " entries.");

  for (const auto i : index_range(vf))
    for (const auto j : index_range(vf[i]))
      if (vf[i][j] != 0.0)
        _vf_matrix.add(i, j, vf[i][j]);
}

const ViewFactorMatrix *
ConstantViewFactorSurfaceRadiation::setViewFactors()
{
  return &_vf_matrix;
}

void
ConstantViewFactorSurfaceRadiation::initialize()
{
  GrayLambertSurfaceRadiationBase::initialize();

  // check row-sum and normalize if necessary
  for (const auto i : make_range(_n_sides))
  {
    Real sum = 0;
    for (const auto & [j, v] : _vf_matrix.row(i))
      sum += v;

    // an error of 5% is acceptable, but more indicates an error in the
    // problem setup
    if (std::abs(sum - 1) > 0.05)
      mooseError("view_factors row ", i, " sums to ", sum);

    // correct view factors
    const Real factor = 1. / sum;
    for (auto & [j, v] : _vf_matrix.row(i))
      v *= factor;
  }
}
