//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "SideUserObject.h"

#ifdef MOOSE_XDG_ENABLED
#include "RayTracingAngularQuadrature.h"
#include "xdg/xdg.h"

#include "libmesh/fe_base.h"
#include "libmesh/quadrature.h"
#include "libmesh/replicated_mesh.h"

#include <map>
#include <memory>
#include <vector>

/**
 * Computes the view factors and escape fractions between surfaces whose view of
 * one another is unobstructed by an enclosure: the surfaces radiate into an open
 * environment (vacuum) through the XDG ray tracing library.
 *
 * The mesh (which must be three-dimensional and serial/replicated) is passed to
 * XDG as a single-volume model. XDG builds a surface for each sideset and an
 * implicit complement volume for the void around the mesh. For every face and
 * angular quadrature direction of a participating boundary, a ray is fired into
 * the void; the first surface it hits receives the associated weight
 * (JxW * (normal . direction) * angular weight), a miss contributes the weight
 * to the escape fraction of the emitting boundary. The weights are normalized
 * by the boundary area times pi in VacuumRayViewFactor.
 *
 * This object is not a RayTracingStudy and does not use the ray tracing module
 * of MOOSE; ray tracing is performed by XDG/Embree.
 */
class VacuumRayViewFactorStudy : public SideUserObject
{
public:
  static InputParameters validParams();

  VacuumRayViewFactorStudy(const InputParameters & parameters);

  void initialSetup() override;
  void initialize() override;
  void execute() override;
  void finalize() override;

  /// The MPI-summed weight of rays leaving \p from_id that arrive at \p to_id
  Real viewFactorWeight(BoundaryID from_id, BoundaryID to_id) const;

  /// The MPI-summed weight of rays leaving \p from_id that escape to the environment
  Real escapeWeight(BoundaryID from_id) const;

protected:
  void threadJoin(const UserObject & y) override;

  /// The local index of a participating boundary id
  unsigned int boundaryIndex(BoundaryID id) const;

  /// the participating boundary ids, in the order given by the boundary parameter
  std::vector<BoundaryID> _boundary_ids;

  /// map from the participating boundary id to its local index
  std::map<BoundaryID, unsigned int> _boundary_index;

  /// number of participating boundaries
  unsigned int _n_boundaries;

  /// face finite element used to integrate the leaving radiation
  std::unique_ptr<FEBase> _fe_face;

  /// face quadrature rule
  std::unique_ptr<QBase> _q_face;

  /// angular quadrature of the leaving hemisphere, one per threaded copy
  std::unique_ptr<RayTracingAngularQuadrature> _aq;

  /// number of angular directions
  std::size_t _num_dir;

  /// the view factor weights accumulated on this thread: [from][to]
  std::vector<std::vector<Real>> _vf_weights;

  /// the escape weights accumulated on this thread: [from]
  std::vector<Real> _escape_weights;

  /// the XDG model shared by all threads; built in initialSetup on thread 0
  std::shared_ptr<xdg::XDG> _xdg;

  /// the implicit complement volume of the model, i.e. the void the rays travel through
  xdg::MeshID _void_volume;

  /// map from an XDG surface id to the participating boundary id
  std::map<xdg::MeshID, BoundaryID> _surface_to_boundary;

  /// map from the sorted node ids of a face to its XDG face id, used to exclude
  /// the emitting face from the ray trace
  std::map<std::vector<dof_id_type>, xdg::MeshID> _face_id_by_node_ids;

  /// distance by which ray origins are nudged off the emitting surface
  Real _bump;
};

#endif // MOOSE_XDG_ENABLED