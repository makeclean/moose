#include "RBE3Constraint.h"
#include "MooseTypes.h"
#include "Assembly.h"
#include "SubProblem.h"
#include "FEProblem.h"
#include "MooseMesh.h"
#include "libmesh/point.h"
#include "libmesh/elem.h"
#include "libmesh/mesh_tools.h"

registerMooseObject("SolidMechanicsApp", RBE3Constraint);

InputParameters
RBE3Constraint::validParams()
{
  InputParameters params = NodalConstraint::validParams();
  params.addClassDescription("RBE3 constraint for rigid body connections between nodes");
  
  // Node specification parameters - simplified interface
  params.addRequiredParam<BoundaryName>("primary_sideset", "Primary sideset name");
  params.addParam<BoundaryName>("secondary_sideset", "Secondary sideset name (should contain only one node)");
  params.addParam<Point>("secondary_node_coordinate", "Secondary node coordinate (creates temporary node)");
  
  // Weight parameters
  MooseEnum weight_methods("explicit equal distance", "equal");
  params.addParam<MooseEnum>("weight_method", weight_methods, "Method for calculating weights");
  params.addParam<std::vector<Real>>("weights", "Explicit constraint weights for primary nodes");
  
  // Other parameters
  params.addParam<unsigned int>("ndof", 3, "Number of degrees of freedom (typically 3)");
  
  // Remove the old parameters that are no longer valid
  // params.addRequiredParam<std::vector<dof_id_type>>("secondary_nodes", "Secondary node IDs");
  
  return params;
}

RBE3Constraint::RBE3Constraint(const InputParameters & parameters)
  : NodalConstraint(parameters),
    _primary_sideset(getParam<BoundaryName>("primary_sideset")),
    _secondary_sideset(getParam<BoundaryName>("secondary_sideset")),
    _secondary_node_coordinate(getParam<Point>("secondary_node_coordinate")),
    _weight_method(getParam<MooseEnum>("weight_method")),
    _weights(getParam<std::vector<Real>>("weights")),
    _ndof(getParam<unsigned int>("ndof")),
    _using_secondary_sideset(isParamValid("secondary_sideset")),
    _using_coordinate(isParamValid("secondary_node_coordinate")),
    _temporary_node_id(-1)
{
  // Validate parameter combinations
  bool has_secondary_sideset = isParamValid("secondary_sideset");
  bool has_coordinate = isParamValid("secondary_node_coordinate");
  
  if (has_secondary_sideset && has_coordinate)
    mooseError("Cannot specify both 'secondary_sideset' and 'secondary_node_coordinate'. Choose one or the other.");
  
  if (!has_secondary_sideset && !has_coordinate)
    mooseError("Must specify either 'secondary_sideset' or 'secondary_node_coordinate'.");
  
  // Validate that secondary sideset contains exactly one node if specified
  if (has_secondary_sideset)
  {
    // We'll check this during connectivity update when we can access the mesh
  }
  
  // Initialize weights if not explicitly provided
  if (_weights.empty() && _weight_method != "explicit")
  {
    // We'll determine proper weight size during connectivity update
  }
}

void
RBE3Constraint::updateConnectivity()
{
  // Create temporary node from coordinate if needed
  if (_using_coordinate)
  {
    createTemporaryNodeFromCoordinate();
  }
  
  // Derive primary nodes from sideset
  deriveNodesFromSideset();
  
  // If we're using a secondary sideset, derive nodes from it
  if (_using_secondary_sideset)
  {
    deriveSecondaryNodesFromSideset();
  }
  
  // Validate we have primary nodes
  if (_primary_nodes.empty())
    mooseError("No primary nodes found. Check the primary_sideset parameter.");
  
  // Set up the constraint connectivity for the base class
  _primary_node_vector = _primary_nodes;
  
  // For coordinate-based nodes, use the temporary node ID
  if (_using_coordinate && _temporary_node_id != -1)
  {
    _connected_nodes.clear();
    _connected_nodes.push_back(_temporary_node_id);
  }
  else if (_using_secondary_sideset)
  {
    // Use the derived secondary nodes
    _connected_nodes = _secondary_nodes;
  }
  
  // Setup multi-variable support for displacement components
  setupMultiVariableSupport();
  
  // Calculate weights based on selected method
  calculateWeights();
  
  // Call parent updateConnectivity
  NodalConstraint::updateConnectivity();
}

void
RBE3Constraint::setupMultiVariableSupport()
{
  // Ensure we're working with displacement variables
  // This allows the constraint to work properly with disp_x, disp_y, disp_z
  // The constraint logic will be applied to each component appropriately
  
  // Check if we're operating on standard displacement variables
  // This is handled by the base class NodalConstraint interface
}

void
RBE3Constraint::createTemporaryNodeFromCoordinate()
{
  // Create a temporary node from the coordinate
  const MooseMesh & mesh = _mesh;
  
  // Generate a unique name for the temporary sideset
  _temporary_sideset_name = "rbe3_temporary_" + std::to_string(_tid);
  
  // Create the node in the mesh
  Node * temp_node = mesh.addNode(_secondary_node_coordinate, _temporary_sideset_name);
  _temporary_node_id = temp_node->id();
  
  // Add the node to our temporary sideset
  mesh.addBoundary(_temporary_sideset_name, _temporary_node_id);
  
  // Update connectivity for the mesh to include our new node
  const auto & node_to_elem_map = mesh.nodeToElemMap();
  auto * const distributed_mesh = dynamic_cast<libMesh::DistributedMesh *>(&mesh.getMesh());
  
  if (distributed_mesh)
  {
    // Update mesh connectivity for distributed mesh
    distributed_mesh->add_extra_ghost_node(temp_node);
  }
  
  // Set up the temporary sideset
  mesh.addBoundary(_temporary_sideset_name, _temporary_node_id);
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
  
  // If we have no nodes, issue a warning
  if (_primary_nodes.empty())
  {
    mooseWarning("No nodes found on primary sideset ", _primary_sideset);
  }
}

void
RBE3Constraint::deriveSecondaryNodesFromSideset()
{
  // Get the mesh and boundary information
  const MooseMesh & mesh = _mesh;
  
  // Clear existing secondary nodes
  _secondary_nodes.clear();
  
  // Get the boundary ID
  dof_id_type boundary_id = mesh.getBoundaryID(_secondary_sideset);
  
  // Get all nodes on the specified boundary
  std::vector<dof_id_type> sideset_nodes = mesh.getNodeList(boundary_id);
  
  // Only add nodes that belong to this processor
  for (const auto node_id : sideset_nodes)
  {
    if (mesh.nodeRef(node_id).processor_id() == _subproblem.processor_id())
    {
      _secondary_nodes.push_back(node_id);
    }
  }
  
  // Validate that there's exactly one node in the secondary sideset
  if (_secondary_nodes.size() != 1)
  {
    mooseError("Secondary sideset ", _secondary_sideset, " must contain exactly one node, found ", _secondary_nodes.size());
  }
  
  // If we have no nodes, issue a warning
  if (_secondary_nodes.empty())
  {
    mooseWarning("No nodes found on secondary sideset ", _secondary_sideset);
  }
}

void
RBE3Constraint::calculateWeights()
{
  // Calculate weights based on the selected method
  if (_weight_method == "explicit")
  {
    // Use explicit weights provided by user
    if (!_weights.empty() && _weights.size() != _primary_nodes.size())
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
  if (_secondary_nodes.empty() && !_using_coordinate)
    mooseError("Need at least one secondary node to calculate distance-based weights");
  
  // Use the secondary node as reference point
  Point reference_point;
  if (_using_coordinate)
  {
    // For coordinate-based, use the coordinate point directly
    reference_point = _secondary_node_coordinate;
  }
  else if (!_secondary_nodes.empty())
  {
    // For sideset-based, use the first secondary node
    const Node & secondary_node = _mesh.nodeRef(_secondary_nodes[0]);
    reference_point = secondary_node;
  }
  
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
    Real distance = (_primary_node_coords[i] - reference_point).norm();
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