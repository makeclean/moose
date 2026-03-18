//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RBE2Constraint.h"
#include "MooseMesh.h"
#include "MooseVariable.h"
#include "SystemBase.h"

registerMooseObject("MooseApp", RBE2Constraint);

InputParameters
RBE2Constraint::validParams()
{
  InputParameters params = NodalConstraint::validParams();
  params.addClassDescription(
      "RBE2 (Rigid Body Element 2) constraint - constrains slave nodes to follow master node "
      "motion. The constraint enforces: u_slave = u_master + theta x r, where r is the vector "
      "from master to slave and theta is the rotation.");
  params.addRequiredParam<std::vector<dof_id_type>>("master_node_ids",
                                                    "The master node IDs (reference nodes)");
  params.addParam<std::vector<dof_id_type>>(
      "slave_node_ids", {}, "The list of slave node IDs");
  params.addParam<BoundaryName>(
      "master_node_set", "NaN", "The boundary ID associated with the master node set");
  params.addParam<BoundaryName>(
      "slave_node_set", "NaN", "The boundary ID associated with the slave node set");
  params.addRequiredParam<Real>("penalty", "The penalty used for the boundary term");
  params.addParam<bool>("use_automatic_weights", true,
                        "If true, weights are computed automatically based on geometry. "
                        "If false, weights must be provided manually.");
  params.addParam<std::vector<Real>>("weights", {},
                                     "Manual weights for master nodes (required if "
                                     "use_automatic_weights is false)");
  return params;
}

RBE2Constraint::RBE2Constraint(const InputParameters & parameters)
  : NodalConstraint(parameters),
    _master_node_ids(getParam<std::vector<dof_id_type>>("master_node_ids")),
    _slave_node_ids(getParam<std::vector<dof_id_type>>("slave_node_ids")),
    _master_node_set_id(getParam<BoundaryName>("master_node_set")),
    _slave_node_set_id(getParam<BoundaryName>("slave_node_set")),
    _penalty(getParam<Real>("penalty")),
    _use_automatic_weights(getParam<bool>("use_automatic_weights"))
{
  const auto & lm_mesh = _mesh.getMesh();

  if (_master_node_ids.empty() && _master_node_set_id == "NaN")
    mooseError("Please specify master_node_ids or master_node_set.");

  if (_slave_node_ids.empty() && _slave_node_set_id == "NaN")
    mooseError("Please specify slave_node_ids or slave_node_set.");

  std::vector<dof_id_type> master_nodes;
  if (!_master_node_ids.empty())
  {
    master_nodes = _master_node_ids;
  }
  else
  {
    master_nodes = _mesh.getNodeList(_mesh.getBoundaryID(_master_node_set_id));
  }

  if (master_nodes.empty())
    mooseError("No master nodes specified.");

  if (master_nodes.size() > 1)
    mooseWarning("Multiple master nodes specified. Using the first one for RBE2.");

  _master_position = Point(0, 0, 0);
  const Node * master_node = lm_mesh.query_node_ptr(master_nodes[0]);
  if (master_node)
    _master_position = *master_node;

  const auto & node_to_elem_map = _mesh.nodeToElemMap();
  auto node_to_elem_pair = node_to_elem_map.find(master_nodes[0]);
  if (node_to_elem_pair != node_to_elem_map.end())
  {
    const std::vector<dof_id_type> & elems = node_to_elem_pair->second;
    for (const auto & elem_id : elems)
      _subproblem.addGhostedElem(elem_id);
  }

  _primary_node_vector.push_back(master_nodes[0]);

  if (!_slave_node_ids.empty())
  {
    for (const auto & dof : _slave_node_ids)
    {
      const Node * const node = lm_mesh.query_node_ptr(dof);
      if (node && node->processor_id() == _subproblem.processor_id())
      {
        _connected_nodes.push_back(dof);
        _slave_positions.push_back(*node);
      }
    }
  }
  else
  {
    std::vector<dof_id_type> nodelist =
        _mesh.getNodeList(_mesh.getBoundaryID(_slave_node_set_id));
    std::vector<dof_id_type>::iterator in;
    for (in = nodelist.begin(); in != nodelist.end(); ++in)
    {
      const Node * const node = lm_mesh.query_node_ptr(*in);
      if (node && node->processor_id() == _subproblem.processor_id())
      {
        _connected_nodes.push_back(*in);
        _slave_positions.push_back(*node);
      }
    }
  }

  if (_use_automatic_weights)
  {
    _weights.resize(1, 1.0);
  }
  else
  {
    _weights = getParam<std::vector<Real>>("weights");
    if (_weights.empty())
      mooseError("Weights must be provided when use_automatic_weights is false.");
  }
}

Real
RBE2Constraint::computeQpResidual(Moose::ConstraintType type)
{
  Point r = _slave_positions[_i] - _master_position;

  Real master_value = _u_primary[_j];

  Real secondary_value = _u_secondary[_i];

  Real constraint_value = master_value - secondary_value;

  switch (type)
  {
    case Moose::Primary:
      return constraint_value * _penalty;
    case Moose::Secondary:
      return -constraint_value * _penalty;
  }
  return 0.;
}

Real
RBE2Constraint::computeQpJacobian(Moose::ConstraintJacobianType type)
{
  switch (type)
  {
    case Moose::PrimaryPrimary:
      return _penalty;
    case Moose::PrimarySecondary:
      return 0.;
    case Moose::SecondarySecondary:
      return _penalty;
    case Moose::SecondaryPrimary:
      return -_penalty;
    default:
      mooseError("Unsupported type");
      break;
  }
  return 0.;
}

Point
RBE2Constraint::getMasterNodePosition() const
{
  return _master_position;
}
