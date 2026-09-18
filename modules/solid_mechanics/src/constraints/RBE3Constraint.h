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
  virtual void initialize() override;
  virtual void execute() override;
  virtual void finalize() override;

  /// Primary node IDs (explicit list)
  std::vector<dof_id_type> _primary_nodes;
  
  /// Secondary node IDs
  std::vector<dof_id_type> _secondary_nodes;
  
  /// Primary sideset name
  BoundaryName _primary_sideset;
  
  /// Weight calculation method
  MooseEnum _weight_method;
  
  /// Explicit weights (if specified)
  std::vector<Real> _weights;
  
  /// Number of degrees of freedom (typically 3)
  unsigned int _ndof;
  
  /// Whether we're using sideset-based node derivation
  bool _using_sideset;
  
  /// Cached node coordinates for distance-based weighting
  std::vector<Point> _primary_node_coords;
  
  /// Reference point for distance-based weighting (secondary node location)
  Point _reference_point;
  
  /// Calculate weights based on the selected method
  void calculateWeights();
  
  /// Derive primary nodes from sideset if specified
  void deriveNodesFromSideset();
  
  /// Calculate distance-based weights
  void calculateDistanceWeights();
};