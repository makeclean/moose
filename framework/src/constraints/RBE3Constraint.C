//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RBE3Constraint.h"
#include "MooseMesh.h"
#include "MooseVariable.h"
#include "SystemBase.h"

registerMooseObject("MooseApp", RBE3Constraint);

InputParameters
RBE3Constraint::validParams()
{
  InputParameters params = NodalConstraint::validParams();
  params.addClassDescription(
      "RBE3 (Rigid Body Element 3) constraint - constrains slave nodes to move as a weighted "
      "combination of master nodes. Supports constant weighting, inverse distance weighting, "
      "or manual weight specification.");
  params.addRequiredParam<std::vector<dof_id_type>>("master_node_ids",
                                                    "The master node IDs");
  params.addParam<BoundaryName>(
      "master_node_set", "NaN", "The boundary ID associated with the master node set");
  params.addParam<BoundaryName>(
      "slave_node_set", "NaN", "The boundary ID associated with the slave node set");
  params.addRequiredParam<Real>("penalty", "The penalty used for the boundary term");
  MooseEnum weighting_type("constant inverse_distance manual");
  params.addParam<MooseEnum>("weighting_type",
                             weighting_type,
                             "constant",
                             "Type of weighting to use: 'constant' gives equal weights to all "
                             "master nodes, 'inverse_distance' uses inverse distance weighting "
                             "from each slave to each master, 'manual' uses user-provided weights.");
  params.addParam<std::vector<Real>>(
      "weights", {},
      "Manual weights for master nodes (required when weighting_type='manual'). "
      "Should be of the same size as master_node_ids.");
  params.addParam<int>("weight_power", 2,
                       "Power for inverse distance weighting (typically 1-3). "
                       "Higher values give more weight to closer nodes.");
  return params;
}

RBE3Constraint::RBE3Constraint(const InputParameters & parameters)
  : NodalConstraint(parameters),
    _master_node_ids(getParam<std::vector<dof_id_type>>("master_node_ids")),
    _master_node_set_id(getParam<BoundaryName>("master_node_set")),
    _slave_node_set_id(getParam<BoundaryName>("slave_node_set")),
    _penalty(getParam<Real>("penalty")),
    _weighting_type(getParam<MooseEnum>("weighting_type")),
    _weight_power(getParam<int>("weight_power"))
{
  if (_master_node_ids.empty() && _master_node_set_id == "NaN")
    mooseError("Please specify master_node_ids or master_node_set.");

  if (_slave_node_set_id == "NaN")
    mooseError("Please specify slave_node_set.");

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

  const auto & node_to_elem_map = _mesh.nodeToElemMap();

  for (const auto & master_id : master_nodes)
  {
    const Node * master_node = _mesh.getMesh().query_node_ptr(master_id);
    if (master_node)
      _master_positions.push_back(*master_node);

    auto node_to_elem_pair = node_to_elem_map.find(master_id);
    if (node_to_elem_pair != node_to_elem_map.end())
    {
      const std::vector<dof_id_type> & elems = node_to_elem_pair->second;
      for (const auto & elem_id : elems)
        _subproblem.addGhostedElem(elem_id);
    }

    _primary_node_vector.push_back(master_id);
  }

  std::vector<dof_id_type> slave_nodes =
      _mesh.getNodeList(_mesh.getBoundaryID(_slave_node_set_id));
  for (const auto & node_id : slave_nodes)
  {
    const Node * const node = _mesh.getMesh().query_node_ptr(node_id);
    if (node && node->processor_id() == _subproblem.processor_id())
      _connected_nodes.push_back(node_id);
  }

  if (_weighting_type == "manual")
  {
    _weights = getParam<std::vector<Real>>("weights");
    if (_weights.size() != _master_node_ids.size())
      mooseError("Weights size must match master_node_ids size when using manual weights.");
  }
  else
  {
    computeWeights();
  }
}

void
RBE3Constraint::computeWeights()
{
  if (_weighting_type == "constant")
  {
    _weights.resize(_master_node_ids.size(), 1.0 / _master_node_ids.size());
    _slave_weights.clear();
    return;
  }

  if (_master_node_ids.size() == 1)
  {
    _weights.resize(1, 1.0);
    _slave_weights.clear();
    return;
  }

  _weights.resize(_master_node_ids.size(), 0.0);
  _slave_weights.resize(_connected_nodes.size());

  for (unsigned int i = 0; i < _connected_nodes.size(); ++i)
  {
    const Node * slave_node = _mesh.getMesh().query_node_ptr(_connected_nodes[i]);
    if (!slave_node)
      continue;

    Point slave_pos = *slave_node;

    Real sum_weights = 0.0;
    std::vector<Real> raw_weights(_master_node_ids.size());

    for (unsigned int j = 0; j < _master_node_ids.size(); ++j)
    {
      Real dist = (slave_pos - _master_positions[j]).norm();
      if (dist < 1e-10)
      {
        raw_weights[j] = 1e10;
        sum_weights = 1e10;
        break;
      }
      raw_weights[j] = 1.0 / std::pow(dist, _weight_power);
      sum_weights += raw_weights[j];
    }

    if (sum_weights > 0)
    {
      for (unsigned int j = 0; j < _master_node_ids.size(); ++j)
      {
        _slave_weights[i].push_back(raw_weights[j] / sum_weights);
      }
    }
  }

  if (_connected_nodes.empty())
    return;

  for (unsigned int j = 0; j < _master_node_ids.size(); ++j)
  {
    Real avg_weight = 0.0;
    for (unsigned int i = 0; i < _slave_weights.size(); ++i)
    {
      if (j < _slave_weights[i].size())
        avg_weight += _slave_weights[i][j];
    }
    if (!_slave_weights.empty())
      avg_weight /= _slave_weights.size();
    _weights[j] = avg_weight;
  }
}

Real
RBE3Constraint::computeQpResidual(Moose::ConstraintType type)
{
  unsigned int primary_size = _primary_node_vector.size();

  Real slave_value = _u_secondary[_i];

  Real weight;
  if (_slave_weights.empty() || _i >= _slave_weights.size() ||
      _j >= _slave_weights[_i].size())
  {
    weight = _weights[_j];
  }
  else
  {
    weight = _slave_weights[_i][_j];
  }

  Real master_contribution = _u_primary[_j] * weight;

  switch (type)
  {
    case Moose::Primary:
      return (master_contribution - slave_value / primary_size) * _penalty;
    case Moose::Secondary:
      return (slave_value / primary_size - master_contribution) * _penalty;
  }
  return 0.;
}

Real
RBE3Constraint::computeQpJacobian(Moose::ConstraintJacobianType type)
{
  unsigned int primary_size = _primary_node_vector.size();

  Real weight;
  if (_slave_weights.empty() || _i >= _slave_weights.size() ||
      _j >= _slave_weights[_i].size())
  {
    weight = _weights[_j];
  }
  else
  {
    weight = _slave_weights[_i][_j];
  }

  switch (type)
  {
    case Moose::PrimaryPrimary:
      return _penalty * weight;
    case Moose::PrimarySecondary:
      return -_penalty / primary_size;
    case Moose::SecondarySecondary:
      return _penalty / primary_size;
    case Moose::SecondaryPrimary:
      return -_penalty * weight;
    default:
      mooseError("Unsupported type");
      break;
  }
  return 0.;
}
