//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"
#include "DerivativeMaterialInterface.h"

class ThermalFlowSensitivity : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  ThermalFlowSensitivity(const InputParameters & parameters);

  virtual void computeQpProperties() override;

protected:
  const VariableGradient & _grad_temperature;
  const VariableValue & _design_density;
  const MaterialPropertyName _design_density_name;
  const MaterialProperty<Real> & _thermal_conductivity;
  const MaterialProperty<Real> & _dthermal_conductivity_ddesign_density;
  MaterialProperty<Real> & _sensitivity;
};
