//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SideSetPatch.h"
#include "MooseError.h"
#include "MooseTypes.h"

#include "libmesh/boundary_info.h"
#include "libmesh/centroid_partitioner.h"
#include "libmesh/elem.h"
#include "libmesh/hilbert_sfc_partitioner.h"
#include "libmesh/linear_partitioner.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/morton_sfc_partitioner.h"
#include "libmesh/replicated_mesh.h"
#include "libmesh/utility.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <set>

SideSetPatch::SideSetPatch(const libMesh::MeshBase & mesh,
                           const libMesh::boundary_id_type sideset,
                           const unsigned int n_patches,
                           const std::string & partitioner,
                           const std::string & centroid_partitioner_direction)
  : _sideset(sideset),
    _partitioner(partitioner),
    _centroid_partitioner_direction(centroid_partitioner_direction),
    _dim(mesh.mesh_dimension() - 1),
    _n_patches(n_patches)
{
  buildBoundaryMesh(mesh);
  partition();

  // record the patch of every face; the processor id of a boundary element is its patch index
  for (const auto & elem : _boundary_mesh->active_element_ptr_range())
  {
    const auto it = _boundary_elem_to_mesh_elem.find(elem->id());
    mooseAssert(it != _boundary_elem_to_mesh_elem.end(),
                "Boundary mesh element not found in boundary_elem_to_mesh_elem.");
    const auto & [elem_id, side, _] = it->second;
    mooseAssert(elem->processor_id() < _n_patches, "Processor id larger than number of patches.");
    _patch_index.insert({{elem_id, side}, elem->processor_id()});
  }
}

void
SideSetPatch::buildBoundaryMesh(const libMesh::MeshBase & mesh)
{
  // create a dim - 1 dimensional mesh
  _boundary_mesh = std::make_unique<libMesh::ReplicatedMesh>(mesh.comm(), _dim);
  _boundary_mesh->set_mesh_dimension(_dim);
  _boundary_mesh->set_spatial_dimension(mesh.mesh_dimension());

  // nodes in the new mesh by boundary_node_id (index)
  std::vector<libMesh::Node *> boundary_nodes;
  // a map from the node numbering on the volumetric mesh to the numbering on the boundary mesh
  std::map<libMesh::dof_id_type, libMesh::dof_id_type> mesh_node_id_to_boundary_node_id;

  libMesh::dof_id_type boundary_node_id = 0;
  for (const auto & side : mesh.get_boundary_info().build_active_side_list())
    if (std::get<2>(side) == _sideset)
    {
      // the original volumetric mesh element
      const libMesh::Elem * elem = mesh.elem_ptr(std::get<0>(side));

      // the boundary element
      std::unique_ptr<const libMesh::Elem> boundary_elem = elem->side_ptr(std::get<1>(side));

      // the boundary node ids of this element in the right order
      std::vector<libMesh::dof_id_type> bnd_elem_node_ids(boundary_elem->n_nodes());

      for (const auto j : make_range(boundary_elem->n_nodes()))
      {
        const libMesh::Node * node = boundary_elem->node_ptr(j);

        const auto node_it = mesh_node_id_to_boundary_node_id.find(node->id());
        if (node_it == mesh_node_id_to_boundary_node_id.end())
        {
          mesh_node_id_to_boundary_node_id.insert({node->id(), boundary_node_id});
          boundary_nodes.push_back(
              _boundary_mesh->add_point(libMesh::Point(*node), boundary_node_id));
          bnd_elem_node_ids[j] = boundary_node_id;
          ++boundary_node_id;
        }
        else
          bnd_elem_node_ids[j] = node_it->second;
      }

      std::unique_ptr<libMesh::Elem> new_bnd_elem = libMesh::Elem::build(boundary_elem->type());
      if (new_bnd_elem->dim() >= 3)
        mooseError("Unsupported element type (libMesh elem_type enum): ", boundary_elem->type());

      libMesh::Elem * const bnd_elem = _boundary_mesh->add_elem(std::move(new_bnd_elem));
      _boundary_elem_to_mesh_elem.insert({bnd_elem->id(), side});

      for (const auto j : make_range(boundary_elem->n_nodes()))
        bnd_elem->set_node(j, boundary_nodes[bnd_elem_node_ids[j]]);
    }

  _boundary_mesh->prepare_for_use();
}

void
SideSetPatch::partition()
{
  if (_n_patches == 0)
  {
    // facet mode: one patch per face; the partitioner is not used
    _n_patches = _boundary_mesh->n_elem();
    unsigned int patch_id = 0;
    for (auto & elem : _boundary_mesh->active_element_ptr_range())
      elem->processor_id() = patch_id++;
  }
  else if (_partitioner == "grid")
    partitionWithGrid();
  else
  {
    if (_partitioner == "default" || _partitioner == "metis" || _partitioner == "parmetis")
      ; // the default partitioner of a replicated mesh is METIS
    else if (_partitioner == "linear")
      _boundary_mesh->partitioner().reset(new libMesh::LinearPartitioner);
    else if (_partitioner == "centroid")
    {
      if (_centroid_partitioner_direction.empty())
        mooseError("If using the centroid partitioner you _must_ specify "
                   "centroid_partitioner_direction!");

      if (_centroid_partitioner_direction == "x")
        _boundary_mesh->partitioner().reset(
            new libMesh::CentroidPartitioner(libMesh::CentroidPartitioner::X));
      else if (_centroid_partitioner_direction == "y")
        _boundary_mesh->partitioner().reset(
            new libMesh::CentroidPartitioner(libMesh::CentroidPartitioner::Y));
      else if (_centroid_partitioner_direction == "z")
        _boundary_mesh->partitioner().reset(
            new libMesh::CentroidPartitioner(libMesh::CentroidPartitioner::Z));
      else if (_centroid_partitioner_direction == "radial")
        _boundary_mesh->partitioner().reset(
            new libMesh::CentroidPartitioner(libMesh::CentroidPartitioner::RADIAL));
      else
        mooseError("Unrecognized centroid partitioner direction ",
                   _centroid_partitioner_direction,
                   ".");
    }
    else if (_partitioner == "hilbert_sfc")
      _boundary_mesh->partitioner().reset(new libMesh::HilbertSFCPartitioner);
    else if (_partitioner == "morton_sfc")
      _boundary_mesh->partitioner().reset(new libMesh::MortonSFCPartitioner);
    else
      mooseError("The ", _partitioner, " partitioner is not supported for splitting sidesets.");

    _boundary_mesh->partition(_n_patches);
  }

  // make sure every partition has at least one element; if not rename and adjust _n_patches
  checkPartitionAndCompress();
}

void
SideSetPatch::partitionWithGrid()
{
  auto bounding_box = libMesh::MeshTools::create_bounding_box(*_boundary_mesh);
  const auto & min = bounding_box.min();
  const auto & max = bounding_box.max();
  const auto & delta = max - min;

  std::vector<unsigned int> nelems(3);
  if (_dim == 1)
  {
    // find the largest component in delta
    unsigned int largest_id = 0;
    Real largest = delta(0);
    for (const auto j : make_range(1, 3))
      if (largest < delta(j))
      {
        largest = delta(j);
        largest_id = j;
      }

    nelems = {1, 1, 1};
    nelems[largest_id] = _n_patches;
  }
  else
  {
    // find the smallest component in delta
    unsigned int smallest_id = 0;
    Real smallest = delta(0);
    for (const auto j : make_range(1, 3))
      if (smallest > delta(j))
      {
        smallest = delta(j);
        smallest_id = j;
      }

    // the ids of the two larger dimensions
    unsigned int id1 = 1, id2 = 2;
    if (smallest_id == 1)
      id1 = 0;
    else if (smallest_id == 2)
      id2 = 0;

    nelems[smallest_id] = 1;
    nelems[id1] = std::round(std::sqrt(delta(id1) / delta(id2) * _n_patches));
    nelems[id2] = std::round(std::sqrt(delta(id2) / delta(id1) * _n_patches));
    const std::size_t final_n_patches = nelems[id1] * nelems[id2];
    if (_n_patches != final_n_patches)
    {
      mooseInfo("For creating radiation patches for sideset ",
                _sideset,
                " using the grid partitioner the number of patches was changed from ",
                _n_patches,
                " to ",
                final_n_patches,
                ".");
      _n_patches = final_n_patches;
    }
  }

  const libMesh::Point delem(delta(0) / nelems[0], delta(1) / nelems[1], delta(2) / nelems[2]);
  for (auto & elem : _boundary_mesh->active_element_ptr_range())
  {
    const libMesh::Point centroid = elem->vertex_average();
    std::array<unsigned int, 3> ielem;
    for (const auto i : make_range(LIBMESH_DIM))
    {
      const auto dist = centroid(i) - min(i);
      ielem[i] = (std::abs(dist) < (libMesh::TOLERANCE * libMesh::TOLERANCE))
                     ? 0
                     : std::floor(dist / delem(i));
    }
    elem->processor_id() = ielem[0] + nelems[0] * (ielem[1] + ielem[2] * nelems[1]);
    mooseAssert(elem->processor_id() < _n_patches, "Invalid processor id");
  }
}

void
SideSetPatch::checkPartitionAndCompress()
{
  std::set<libMesh::processor_id_type> processor_ids;
  for (auto & elem : _boundary_mesh->active_element_ptr_range())
    processor_ids.insert(elem->processor_id());

  if (processor_ids.size() == _n_patches)
    return;

  // at least one partition does not have an element assigned to it: adjust _n_patches
  mooseInfo("Some partitions for sideset ",
            _sideset,
            " are empty. Adjusting number of patches from ",
            _n_patches,
            " to ",
            processor_ids.size(),
            ".");
  _n_patches = processor_ids.size();

  std::vector<libMesh::processor_id_type> processor_ids_vec(processor_ids.begin(),
                                                            processor_ids.end());
  std::sort(processor_ids_vec.begin(), processor_ids_vec.end());

  std::map<libMesh::processor_id_type, libMesh::processor_id_type> processor_id_remap;
  for (const auto j : index_range(processor_ids_vec))
    processor_id_remap[processor_ids_vec[j]] = j;

  for (auto & elem : _boundary_mesh->active_element_ptr_range())
  {
    const auto it = processor_id_remap.find(elem->processor_id());
    mooseAssert(it != processor_id_remap.end(), "Processor id not found in remap.");
    elem->processor_id() = it->second;
  }
}

unsigned int
SideSetPatch::patchIndex(const libMesh::dof_id_type elem_id, const unsigned int side) const
{
  const auto it = _patch_index.find({elem_id, side});
  return it == _patch_index.end() ? libMesh::invalid_uint : it->second;
}
