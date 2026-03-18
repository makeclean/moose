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

class RBE2Constraint : public NodalConstraint
{
public:
  static InputParameters validParams();

  RBE2Constraint(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::ConstraintType type) override;
  virtual Real computeQpJacobian(Moose::ConstraintJacobianType type) override;

  Point getMasterNodePosition() const;

  std::vector<dof_id_type> _master_node_ids;
  std::vector<dof_id_type> _slave_node_ids;
  std::string _master_node_set_id;
  std::string _slave_node_set_id;
  Real _penalty;

  std::vector<std::vector<Real>> _slave_positions;
  Point _master_position;

  bool _use_automatic_weights;
};
