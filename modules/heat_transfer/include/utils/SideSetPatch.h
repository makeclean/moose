//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "libmesh/boundary_info.h"
#include "libmesh/mesh_base.h"

#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <utility>

/**
 * Splits the faces of one sideset into patches.
 *
 * The faces are put in a replicated lower dimensional boundary mesh which is partitioned into
 * \c n_patches pieces, exactly like PatchSidesetGenerator does when it materializes the patches as
 * sidesets. This class is shared by PatchSidesetGenerator and by the radiation objects, which
 * index the patches internally instead of adding one sideset per patch to the mesh. The latter
 * matters because libMesh boundary ids are 16 bit by default, so a sideset with many faces (facet
 * mode, \c n_patches == 0) would exhaust the boundary id space.
 *
 * \c n_patches == 0 selects facet mode, i.e. one patch per face.
 */
class SideSetPatch
{
public:
  /// A face of the source mesh: (element id, local side, sideset)
  using Face = std::tuple<libMesh::dof_id_type, unsigned short int, libMesh::boundary_id_type>;

  SideSetPatch(const libMesh::MeshBase & mesh,
               const libMesh::boundary_id_type sideset,
               const unsigned int n_patches,
               const std::string & partitioner,
               const std::string & centroid_partitioner_direction);

  /// Number of patches. This can be smaller than the requested number if partitions are empty.
  unsigned int nPatches() const { return _n_patches; }

  /// Patch index of the given side, or libMesh::invalid_uint if the side is not in the sideset
  unsigned int patchIndex(const libMesh::dof_id_type elem_id, const unsigned int side) const;

  /// The lower dimensional boundary mesh the patches were computed on
  libMesh::MeshBase & boundaryMesh() const { return *_boundary_mesh; }

  /// Map from the id of a boundary mesh element to the face it was built from
  const std::map<libMesh::dof_id_type, Face> & boundaryElemToMeshElem() const
  {
    return _boundary_elem_to_mesh_elem;
  }

  /// Map from (element id, side) to patch index for every face of the sideset
  const std::map<std::pair<libMesh::dof_id_type, unsigned int>, unsigned int> & patchIndices() const
  {
    return _patch_index;
  }

protected:
  /// Build a replicated boundary mesh from the faces of the sideset
  void buildBoundaryMesh(const libMesh::MeshBase & mesh);

  /// Partition the boundary mesh into _n_patches pieces
  void partition();

  /// Compress empty partitions and renumber the remaining ones
  void checkPartitionAndCompress();

  /// Partition the boundary mesh with the grid partitioner
  void partitionWithGrid();

  const libMesh::boundary_id_type _sideset;
  const std::string _partitioner;
  const std::string _centroid_partitioner_direction;

  /// dimensionality of the sideset
  const unsigned int _dim;

  /// number of patches, updated in facet mode and when partitions are empty
  unsigned int _n_patches;

  std::unique_ptr<libMesh::MeshBase> _boundary_mesh;

  std::map<libMesh::dof_id_type, Face> _boundary_elem_to_mesh_elem;

  /// (element id, side) to patch index
  std::map<std::pair<libMesh::dof_id_type, unsigned int>, unsigned int> _patch_index;
};
