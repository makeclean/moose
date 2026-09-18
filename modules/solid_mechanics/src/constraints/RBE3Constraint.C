#include "RBE3Constraint.h"
#include "MooseTypes.h"
#include "Assembly.h"
#include "SubProblem.h"
#include "FEProblem.h"
#include "MooseMesh.h"
#include "libmesh/point.h"
#include "libmesh/elem.h"

registerMooseObject("SolidMechanicsApp", RBE3Constraint);

template <>
InputParameters
validParams<RBE3Constraint>()
{
  InputParameters params = validParams<NodalConstraint>();
  params.addClassDescription("RBE3 constraint for rigid body connections between nodes");
  
  // Node specification parameters
  params.addParam<std::vector<dof_id_type>>("primary_nodes", "Primary node IDs");
  params.addParam<BoundaryName>("primary_sideset", "Primary sideset name");
  params.addRequiredParam<std::vector<dof_id_type>>("secondary_nodes", "Secondary node IDs");
  
  // Weight parameters
  MooseEnum weight_methods("explicit equal distance", "equal");
  params.addParam<MooseEnum>("weight_method", weight_methods, "Method for calculating weights");
  params.addParam<std::vector<Real>>("weights", "Explicit constraint weights for primary nodes");
  
  // Other parameters
  params.addParam<unsigned int>("ndof", 3, "Number of degrees of freedom (typically 3)");
  
  // Validate that either primary_nodes or primary_sideset is specified, but not both
  params.addCoupledVar("primary_nodes", "Primary node IDs");
  params.addCoupledVar("primary_sideset", "Primary sideset name");
  
  return params;
}

RBE3Constraint::RBE3Constraint(const InputParameters & parameters)
  : NodalConstraint(parameters),
    _primary_nodes(getParam<std::vector<dof_id_type>>("primary_nodes")),
    _secondary_nodes(getParam<std::vector<dof_id_type>>("secondary_nodes")),
    _primary_sideset(getParam<BoundaryName>("primary_sideset")),
    _weight_method(getParam<MooseEnum>("weight_method")),
    _weights(getParam<std::vector<Real>>("weights")),
    _ndof(getParam<unsigned int>("ndof")),
    _using_sideset(isParamValid("primary_sideset"))
{
  // Validate node specification parameters
  bool has_nodes = isParamValid("primary_nodes");
  bool has_sideset = isParamValid("primary_sideset");
  
  if (has_nodes && has_sideset)
    mooseError("Cannot specify both 'primary_nodes' and 'primary_sideset'. Choose one or the other.");
  
  if (!has_nodes && !has_sideset)
    mooseError("Must specify either 'primary_nodes' or 'primary_sideset'.");
  
  if (_secondary_nodes.empty())
    mooseError("RBE3 constraint must have at least one secondary node");
  
  // Initialize weights if not explicitly provided
  if (_weights.empty() && _weight_method != "explicit")
  {
    _weights.resize(_primary_nodes.size(), 1.0);
  }
}

void
RBE3Constraint::initialize()
{
  // Derive primary nodes from sideset if needed
  if (_using_sideset)
  {
    deriveNodesFromSideset();
  }
  
  // Validate we have primary nodes
  if (_primary_nodes.empty())
    mooseError("No primary nodes found. Check the primary_sideset or primary_nodes parameters.");
  
  // Set up the constraint connectivity for the base class
  _primary_node_vector = _primary_nodes;
  _connected_nodes = _secondary_nodes;
  
  // Calculate weights based on selected method
  calculateWeights();
}

void
RBE3Constraint::execute()
{
  // Execute base class functionality
  NodalConstraint::execute();
}

void
RBE3Constraint::finalize()
{
  // Finalize base class functionality
  NodalConstraint::finalize();
}

void
RBE3Constraint::updateConnectivity()
{
  // Update connectivity for the constraint
  NodalConstraint::updateConnectivity();
}

void
RBE3Constraint::deriveNodesFromSideset()
{
  // Get the mesh and boundary information
  const MooseMesh & mesh = _mesh;
  
  // Clear existing nodes
  _primary_nodes.clear();
  
  // Get the boundary ID
  dof_id_type boundary_id = mesh.getBoundaryID(_primary_sideset);
  
  // Get all nodes on the specified boundary
  std::vector<dof_id_type> sideset_nodes = mesh.getNodeList(boundary_id);
  
  // Only add nodes that belong to this processor
  for (const auto node_id : sideset_nodes)
  {
    if (mesh.nodeRef(node_id).processor_id() == _subproblem.processor_id())
    {
      _primary_nodes.push_back(node_id);
    }
  }
  
  // If we have no nodes, issue an error
  if (_primary_nodes.empty())
  {
    mooseWarning("No nodes found on sideset ", _primary_sideset);
  }
}

void
RBE3Constraint::calculateWeights()
{
  // Calculate weights based on the selected method
  if (_weight_method == "explicit")
  {
    // Use explicit weights provided by user
    if (_weights.size() != _primary_nodes.size())
      mooseError("Number of weights must match number of primary nodes when using 'explicit' method");
  }
  else if (_weight_method == "equal")
  {
    // Assign equal weights to all primary nodes
    _weights.assign(_primary_nodes.size(), 1.0 / _primary_nodes.size());
  }
  else if (_weight_method == "distance")
  {
    // Calculate distance-based weights
    calculateDistanceWeights();
  }
}

void
RBE3Constraint::calculateDistanceWeights()
{
  // For distance-based weighting, we need to calculate distances from secondary node to primary nodes
  if (_secondary_nodes.empty())
    mooseError("Need at least one secondary node to calculate distance-based weights");
  
  // Use the first secondary node as reference point
  const Node & secondary_node = _mesh.nodeRef(_secondary_nodes[0]);
  _reference_point = secondary_node;
  
  // Clear existing weights and prepare for new calculation
  _weights.clear();
  _weights.resize(_primary_nodes.size(), 0.0);
  
  // Calculate distances and weights for each primary node
  std::vector<Real> distances(_primary_nodes.size(), 0.0);
  Real total_distance = 0.0;
  
  // Get coordinates of primary nodes
  _primary_node_coords.clear();
  _primary_node_coords.resize(_primary_nodes.size());
  
  for (unsigned int i = 0; i < _primary_nodes.size(); ++i)
  {
    const Node & primary_node = _mesh.nodeRef(_primary_nodes[i]);
    _primary_node_coords[i] = primary_node;
    
    // Calculate distance from reference point to primary node
    Real distance = (_primary_node_coords[i] - _reference_point).norm();
    distances[i] = distance;
    total_distance += distance;
  }
  
  // Calculate weights (inverse of distance - closer nodes get higher weights)
  // Handle case where distances are zero
  for (unsigned int i = 0; i < _primary_nodes.size(); ++i)
  {
    if (distances[i] > 0.0)
      _weights[i] = 1.0 / distances[i];
    else
      _weights[i] = 1.0;  // If distance is zero, set to 1.0 for normalization
  }
  
  // Normalize weights
  Real weight_sum = 0.0;
  for (const auto weight : _weights)
    weight_sum += weight;
  
  if (weight_sum > 0.0)
  {
    for (auto & weight : _weights)
      weight /= weight_sum;
  }
  else
  {
    // If all distances are zero, assign equal weights
    for (auto & weight : _weights)
      weight = 1.0 / _weights.size();
  }
}

Real
RBE3Constraint::computeQpResidual(Moose::ConstraintType type)
{
  // For RBE3 constraints, we enforce that the secondary node's displacement
  // is a linear combination of the primary nodes' displacements
  // Residual = u_secondary - sum(w_i * u_primary_i)
  
  // The constraint equation is: u_secondary = sum(w_i * u_primary_i)
  // So the residual is: u_secondary - sum(w_i * u_primary_i) = 0
  
  // Get current solution values
  Real secondary_value = 0.0;
  Real primary_sum = 0.0;
  
  if (type == Moose::Secondary)
  {
    // For secondary nodes, compute the residual for the constraint
    secondary_value = _u_secondary[_i];
    for (unsigned int j = 0; j < _primary_nodes.size(); ++j)
    {
      if (j < _weights.size())
        primary_sum += _weights[j] * _u_primary[_j];
    }
    return secondary_value - primary_sum;
  }
  else if (type == Moose::Primary)
  {
    // For primary nodes, we're essentially computing the negative of the constraint
    secondary_value = _u_secondary[_i];
    for (unsigned int j = 0; j < _primary_nodes.size(); ++j)
    {
      if (j < _weights.size())
        primary_sum += _weights[j] * _u_primary[_j];
    }
    return primary_sum - secondary_value;
  }
  
  return 0.0;
}

Real
RBE3Constraint::computeQpJacobian(Moose::ConstraintJacobianType type)
{
  // For Jacobian, we compute derivatives of the constraint equation:
  // J = d/dx(u_secondary - sum(w_i * u_primary_i))
  
  switch (type)
  {
    case Moose::SecondarySecondary:
      return 1.0;
      
    case Moose::PrimaryPrimary:
      // This represents d/dx(sum(w_i * u_primary_i)) = w_i
      if (_j < _weights.size())
        return -_weights[_j];
      else
        return 0.0;
      
    case Moose::SecondaryPrimary:
      return -1.0;
      
    case Moose::PrimarySecondary:
      return 1.0;
      
    default:
      return 0.0;
  }
}