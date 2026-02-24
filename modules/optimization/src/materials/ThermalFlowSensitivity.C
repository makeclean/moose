//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ThermalFlowSensitivity.h"

registerMooseObject("OptimizationApp", ThermalFlowSensitivity);

InputParameters
ThermalFlowSensitivity::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Computes thermal compliance sensitivity with respect to design "
                             "density for topology optimization of heat exchangers.");
  params.addRequiredCoupledVar("temperature", "Temperature variable");
  params.addRequiredCoupledVar("design_density", "Design density variable name.");
  params.addRequiredParam<MaterialPropertyName>(
      "thermal_conductivity",
      "DerivativeParsedMaterial for thermal conductivity as function of design density.");
  return params;
}

ThermalFlowSensitivity::ThermalFlowSensitivity(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _grad_temperature(coupledGradient("temperature")),
    _design_density(coupledValue("design_density")),
    _design_density_name(coupledName("design_density", 0)),
    _thermal_conductivity(
        getMaterialProperty<Real>(getParam<MaterialPropertyName>("thermal_conductivity"))),
    _dthermal_conductivity_ddesign_density(
        getMaterialPropertyDerivativeByName<Real>(
            getParam<MaterialPropertyName>("thermal_conductivity"), _design_density_name)),
    _sensitivity(declareProperty<Real>("thermal_flow_sensitivity"))
{
}

void
ThermalFlowSensitivity::computeQpProperties()
{
  Real grad_T_squared = _grad_temperature[_qp] * _grad_temperature[_qp];
  _sensitivity[_qp] = _dthermal_conductivity_ddesign_density[_qp] * grad_T_squared;
}
