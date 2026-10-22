//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SpecifiedViewFactor.h"

registerMooseObject("HeatTransferApp", SpecifiedViewFactor);

InputParameters
SpecifiedViewFactor::validParams()
{
  InputParameters params = ViewFactorBase::validParams();
  params.addRequiredParam<std::vector<std::vector<Real>>>(
      "view_factors", "The view factors from sideset i to sideset j.");
  params.addClassDescription("View factors specified directly in the input file");
  return params;
}

SpecifiedViewFactor::SpecifiedViewFactor(const InputParameters & parameters)
  : ViewFactorBase(parameters)
{
  const auto & vf = getParam<std::vector<std::vector<Real>>>("view_factors");
  checkViewFactors(vf);

  for (const auto i : index_range(vf))
    for (const auto j : index_range(vf[i]))
      if (vf[i][j] != 0.0)
        _view_factors.add(i, j, vf[i][j]);
}

void
SpecifiedViewFactor::checkViewFactors(const std::vector<std::vector<Real>> & view_factors) const
{
  // check that the input has the right format
  if (view_factors.size() != _n_sides)
    paramError("view_factors",
               "Leading dimension of view_factors must be equal to number of side sets.");

  for (const auto i : index_range(view_factors))
    if (view_factors[i].size() != _n_sides)
      paramError("view_factors",
                 "view_factors must be provided as square array. Row ",
                 i,
                 " has ",
                 view_factors[i].size(),
                 " entries.");
}
