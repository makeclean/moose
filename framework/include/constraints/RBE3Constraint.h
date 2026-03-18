//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "NodalConstraint.h"

class RBE3Constraint : public NodalConstraint
{
public:
  static InputParameters validParams();

  RBE3Constraint(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::ConstraintType type) override;
  virtual Real computeQpJacobian(Moose::ConstraintJacobianType type) override;

  void computeWeights();

  std::vector<dof_id_type> _master_node_ids;
  std::string _master_node_set_id;
  std::string _slave_node_set_id;
  Real _penalty;
  MooseEnum _weighting_type;
  int _weight_power;

  std::vector<Point> _master_positions;
  std::vector<std::vector<Real>> _slave_weights;
};
