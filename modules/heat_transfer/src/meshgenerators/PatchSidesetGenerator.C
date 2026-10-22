//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PatchSidesetGenerator.h"
#include "SideSetPatch.h"
#include "InputParameters.h"
#include "MooseMesh.h"
#include "MooseTypes.h"
#include "MooseUtils.h"
#include "MooseMeshUtils.h"

#include "libmesh/boundary_info.h"
#include "libmesh/elem.h"
#include "libmesh/utility.h"

#include <algorithm>
#include <limits>

registerMooseObject("HeatTransferApp", PatchSidesetGenerator);

InputParameters
PatchSidesetGenerator::validParams()
{
  InputParameters params = MeshGenerator::validParams();

  params.addRequiredParam<MeshGeneratorName>("input", "The mesh we want to modify");
  params.addRequiredParam<BoundaryName>("boundary",
                                        "The boundary that will be divided into patches");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "n_patches",
      "n_patches>=0",
      "Number of patches. If zero, one patch is created for every face of the sideset.");

  params.addParam<bool>(
      "remove_original_sideset",
      false,
      "If true, remove the original sideset after the patch sidesets are created. The patch "
      "sidesets then fully replace the original sideset on the mesh.");

  MooseEnum partitioning = MooseMesh::partitioning(); // default MOOSE partitioning
  partitioning += "grid";                             // ...but also add our own
  params.addParam<MooseEnum>(
      "partitioner",
      partitioning,
      "Specifies a mesh partitioner to use when splitting the mesh for a parallel computation.");

  MooseEnum direction("x y z radial");
  params.addParam<MooseEnum>("centroid_partitioner_direction",
                             direction,
                             "Specifies the sort direction if using the centroid partitioner. "
                             "Available options: x, y, z, radial");

  params.addParamNamesToGroup("partitioner centroid_partitioner_direction", "Partitioning");

  params.addClassDescription(
      "Divides the given sideset into smaller patches of roughly equal size.");

  return params;
}

PatchSidesetGenerator::PatchSidesetGenerator(const InputParameters & parameters)
  : MeshGenerator(parameters),
    _input(getMesh("input")),
    _n_patches(getParam<unsigned int>("n_patches")),
    _sideset_name(getParam<BoundaryName>("boundary")),
    _remove_original_sideset(getParam<bool>("remove_original_sideset"))
{
}

std::unique_ptr<MeshBase>
PatchSidesetGenerator::generate()
{
  std::unique_ptr<MeshBase> mesh = std::move(_input);

  _mesh->errorIfDistributedMesh("PatchSidesetGenerator");

  if (!MooseMeshUtils::hasBoundaryNameOrID(*mesh, _sideset_name))
    paramError("boundary",
               "The provided boundary name or ID (" + _sideset_name +
                   ") does not exist in the mesh.");

  _sideset = MooseMeshUtils::getBoundaryID(_sideset_name, *mesh);

  // split the faces of the sideset into patches; the patch index of a face is the
  // processor id of its element in the boundary mesh
  const SideSetPatch patches(*mesh,
                             _sideset,
                             _n_patches,
                             getParam<MooseEnum>("partitioner"),
                             getParam<MooseEnum>("centroid_partitioner_direction"));
  _n_patches = patches.nPatches();

  // adding one boundary id per patch must not overflow the boundary id type; this is a
  // pre-flight check because MooseMeshUtils::getBoundaryIDs(..., generate_unknown = true)
  // would silently wrap around (and produce negative, invalid ids) otherwise
  const auto & boundary_ids_in_mesh = mesh->get_boundary_info().get_boundary_ids();
  const auto max_boundary_id =
      boundary_ids_in_mesh.empty()
          ? 0
          : *std::max_element(boundary_ids_in_mesh.begin(), boundary_ids_in_mesh.end());
  const auto max_boundary_id_type =
      static_cast<long long>(std::numeric_limits<boundary_id_type>::max());
  if (static_cast<long long>(max_boundary_id) + static_cast<long long>(_n_patches) >
      max_boundary_id_type)
    paramError("n_patches",
               "Splitting sideset ",
               _sideset_name,
               " into ",
               _n_patches,
               " patches would create boundary ids above ",
               max_boundary_id_type,
               ", the largest boundary id this libMesh build supports (the sideset currently ends "
               "at id ",
               static_cast<long long>(max_boundary_id),
               "). Use the radiation transfer action with create_patch_sidesets = false to index "
               "the patches internally, or rebuild libMesh with --with-boundary-id-bytes=4.");

  BoundaryInfo & boundary_info = mesh->get_boundary_info();

  // prepare sideset names and boundary_ids added to mesh
  std::vector<BoundaryName> sideset_names =
      sidesetNameHelper(boundary_info.get_sideset_name(_sideset));

  std::vector<boundary_id_type> boundary_ids =
      MooseMeshUtils::getBoundaryIDs(*mesh, sideset_names, true);

  mooseAssert(sideset_names.size() == _n_patches,
              "sideset_names must have as many entries as user-requested number of patches.");
  mooseAssert(boundary_ids.size() == _n_patches,
              "boundary_ids must have as many entries as user-requested number of patches.");

  // loop over all elements of the boundary mesh and assign the side of the _original_
  // element to the patch sideset
  for (const auto & elem : patches.boundaryMesh().active_element_ptr_range())
  {
    const auto & side = libmesh_map_find(patches.boundaryElemToMeshElem(), elem->id());
    mooseAssert(elem->processor_id() < boundary_ids.size(),
                "Processor id larger than number of patches.");
    boundary_info.add_side(
        std::get<0>(side), std::get<1>(side), boundary_ids[elem->processor_id()]);
  }

  // make sure new boundary names are set
  for (const auto j : index_range(boundary_ids))
  {
    boundary_info.sideset_name(boundary_ids[j]) = sideset_names[j];
    boundary_info.nodeset_name(boundary_ids[j]) = sideset_names[j];
  }

  if (_remove_original_sideset)
  {
    // every face of the original sideset now belongs to exactly one patch sideset,
    // so the original sideset (and its nodeset) can be removed from the mesh
    boundary_info.remove_id(_sideset);
  }

  return mesh;
}

std::vector<BoundaryName>
PatchSidesetGenerator::sidesetNameHelper(const std::string & base_name) const
{
  std::vector<BoundaryName> rv;
  for (const auto j : make_range(_n_patches))
  {
    std::stringstream ss;
    ss << base_name << "_" << j;
    rv.push_back(ss.str());
  }
  return rv;
}
