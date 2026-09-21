#pragma once

// MOOSE includes
#include "NodalConstraint.h"
#include "MooseTypes.h"

class RBE3Constraint : public NodalConstraint
{
public:
  static InputParameters validParams();

  RBE3Constraint(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::ConstraintType type) override;
  virtual Real computeQpJacobian(Moose::ConstraintJacobianType type) override;
  virtual void updateConnectivity() override;

  /// Primary sideset name
  BoundaryName _primary_sideset;
  
  /// Secondary sideset name (should contain only one node)
  BoundaryName _secondary_sideset;
  
  /// Secondary node coordinate (when provided)
  Point _secondary_node_coordinate;
  
  /// Weight calculation method
  MooseEnum _weight_method;
  
  /// Explicit weights (if specified)
  std::vector<Real> _weights;
  
  /// Number of degrees of freedom (typically 3)
  unsigned int _ndof;
  
  /// Whether we're using secondary sideset
  bool _using_secondary_sideset;
  
  /// Whether we're using coordinate-based node creation
  bool _using_coordinate;
  
  /// Temporary node ID created from coordinate
  dof_id_type _temporary_node_id;
  
  /// Temporary sideset name for coordinate-based node
  std::string _temporary_sideset_name;
  
  /// Cached node coordinates for distance-based weighting
  std::vector<Point> _primary_node_coords;
  
  /// Reference point for distance-based weighting (secondary node location)
  Point _reference_point;
  
  /// Calculate weights based on the selected method
  void calculateWeights();
  
  /// Derive primary nodes from sideset if specified
  void deriveNodesFromSideset();
  
  /// Derive secondary nodes from sideset if specified
  void deriveSecondaryNodesFromSideset();
  
  /// Create temporary node from coordinate if needed
  void createTemporaryNodeFromCoordinate();
  
  /// Calculate distance-based weights
  void calculateDistanceWeights();
};