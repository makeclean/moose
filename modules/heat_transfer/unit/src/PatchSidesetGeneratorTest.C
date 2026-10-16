//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "gtest/gtest.h"

#include "ActionWarehouse.h"
#include "Factory.h"
#include "MeshGeneratorMesh.h"
#include "MeshGeneratorSystem.h"
#include "Moose.h"
#include "MooseApp.h"
#include "MooseMain.h"
#include "MooseMesh.h"
#include "MooseMeshUtils.h"
#include "PatchSidesetGenerator.h"

#include "libmesh/boundary_info.h"
#include "libmesh/mesh_base.h"
#include "libmesh/point.h"

#include <algorithm>
#include <memory>
#include <set>
#include <string>
#include <utility>

namespace
{

/// Collects the (element id, side) pairs belonging to the sideset with the given id
std::set<std::pair<dof_id_type, unsigned short>>
sidePairs(const MeshBase & mesh, const BoundaryID id)
{
  std::set<std::pair<dof_id_type, unsigned short>> pairs;
  for (const auto & [elem_id, side, side_set_id] : mesh.get_boundary_info().build_active_side_list())
    if (side_set_id == id)
      pairs.insert({elem_id, side});
  return pairs;
}

/// Returns true if a nodeset with the given name exists on the mesh
bool
hasNodesetName(const MeshBase & mesh, const std::string & name)
{
  const auto & nodeset_names = mesh.get_boundary_info().get_nodeset_name_map();
  return std::any_of(nodeset_names.begin(),
                     nodeset_names.end(),
                     [&name](const auto & id_name) { return id_name.second == name; });
}

} // namespace

class PatchSidesetGeneratorTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    const char * argv[2] = {"foo", "\0"};
    _app = Moose::createMooseApp("HeatTransferApp", 1, (char **)argv);
    _factory = &_app->getFactory();

    // the mesh generator system requires a MooseMesh to be present in the action warehouse
    auto mesh_params = _factory->getValidParams("MeshGeneratorMesh");
    mesh_params.set<MooseEnum>("parallel_type") = "REPLICATED";
    mesh_params.set<std::string>("final_generator") = "patch_gen";
    _moose_mesh =
        _factory->create<MeshGeneratorMesh>("MeshGeneratorMesh", "moose_mesh", mesh_params);
    _moose_mesh->setMeshBase(_moose_mesh->buildMeshBaseObject(3));
    _moose_mesh->buildMesh();
    _app->actionWarehouse().mesh() = _moose_mesh;
  }

  /**
   * Runs the generator chain CartesianMeshGenerator -> SideSetsFromBoundingBoxGenerator ->
   * PatchSidesetGenerator and returns the final mesh
   *
   * The base mesh is a 2x2x2 grid of unit cubes spanning [0, 2]^3, and the input sideset
   * "boundary" consists of the four external faces with normal (1, 0, 0)
   */
  std::unique_ptr<MeshBase>
  generate(const unsigned int n_patches,
           const bool remove_original_sideset,
           const std::string & partitioner = "grid")
  {
    auto & mgs = _app->getMeshGeneratorSystem();

    auto cube_params = _factory->getValidParams("CartesianMeshGenerator");
    cube_params.set<MooseEnum>("dim") = "3";
    cube_params.set<std::vector<Real>>("dx") = {1.0, 1.0};
    cube_params.set<std::vector<Real>>("dy") = {1.0, 1.0};
    cube_params.set<std::vector<Real>>("dz") = {1.0, 1.0};
    mgs.addMeshGenerator("CartesianMeshGenerator", "cube", cube_params);

    auto sideset_params = _factory->getValidParams("SideSetsFromBoundingBoxGenerator");
    sideset_params.set<MeshGeneratorName>("input") = "cube";
    sideset_params.set<BoundaryName>("boundary_new") = "boundary";
    sideset_params.set<RealVectorValue>("bottom_left") = RealVectorValue(0.0, 0.0, 0.0);
    sideset_params.set<RealVectorValue>("top_right") = RealVectorValue(2.0, 2.0, 2.0);
    sideset_params.set<Point>("normal") = Point(1.0, 0.0, 0.0);
    sideset_params.set<bool>("include_only_external_sides") = true;
    mgs.addMeshGenerator("SideSetsFromBoundingBoxGenerator", "boundary_mesh", sideset_params);

    auto patch_params = _factory->getValidParams("PatchSidesetGenerator");
    patch_params.set<MeshGeneratorName>("input") = "boundary_mesh";
    patch_params.set<BoundaryName>("boundary") = "boundary";
    patch_params.set<unsigned int>("n_patches") = n_patches;
    patch_params.set<bool>("remove_original_sideset") = remove_original_sideset;
    patch_params.set<MooseEnum>("partitioner") = partitioner;
    mgs.addMeshGenerator("PatchSidesetGenerator", "patch_gen", patch_params);

    _app->actionWarehouse().executeActionsWithAction("create_added_mesh_generators");
    mgs.createAddedMeshGenerators();
    _app->actionWarehouse().executeActionsWithAction("execute_mesh_generators");
    mgs.executeMeshGenerators();

    return mgs.getSavedMesh("main");
  }

  /**
   * Collects the faces of the patch sidesets <base_name>_0 .. <base_name>_(n-1), expecting each
   * of them to exist on the mesh
   */
  std::set<std::pair<dof_id_type, unsigned short>>
  collectPatchSides(const MeshBase & mesh, const std::string & base_name, const unsigned int n)
  {
    std::set<std::pair<dof_id_type, unsigned short>> patch_sides;
    for (const auto j : make_range(n))
    {
      const BoundaryName patch_name = base_name + std::to_string(j);
      EXPECT_TRUE(MooseMeshUtils::hasBoundaryNameOrID(mesh, patch_name))
          << "missing patch sideset " << patch_name;
      const auto sides = sidePairs(mesh, MooseMeshUtils::getBoundaryID(patch_name, mesh));
      patch_sides.insert(sides.begin(), sides.end());
    }
    return patch_sides;
  }

  std::shared_ptr<MooseApp> _app;
  std::shared_ptr<MeshGeneratorMesh> _moose_mesh;
  Factory * _factory;
};

TEST_F(PatchSidesetGeneratorTest, patchesSubdivideSideset)
{
  const auto mesh = generate(4, false);

  // the parent sideset is retained by default, with all four of its faces intact
  EXPECT_TRUE(MooseMeshUtils::hasBoundaryNameOrID(*mesh, "boundary"));
  const auto parent_sides = sidePairs(*mesh, MooseMeshUtils::getBoundaryID("boundary", *mesh));
  ASSERT_EQ(parent_sides.size(), 4);

  // the four patch sidesets each hold one face and cover the parent faces exactly once
  const auto patch_sides = collectPatchSides(*mesh, "boundary_", 4);
  ASSERT_EQ(patch_sides.size(), 4);
  for (const auto j : make_range(4))
    EXPECT_EQ(sidePairs(*mesh, MooseMeshUtils::getBoundaryID("boundary_" + std::to_string(j),
                                                             *mesh))
                  .size(),
              1);
  EXPECT_EQ(patch_sides, parent_sides);
}

TEST_F(PatchSidesetGeneratorTest, removeOriginalSideset)
{
  const auto mesh = generate(4, true);

  // the parent sideset and its nodeset are fully removed from the mesh
  EXPECT_FALSE(MooseMeshUtils::hasBoundaryNameOrID(*mesh, "boundary"));
  EXPECT_FALSE(hasNodesetName(*mesh, "boundary"));

  // the patch sidesets cover the four original faces exactly once
  const auto patch_sides = collectPatchSides(*mesh, "boundary_", 4);
  EXPECT_EQ(patch_sides.size(), 4);
}

TEST_F(PatchSidesetGeneratorTest, facetModeOnePatchPerFace)
{
  // n_patches = 0: one patch per face; the requested partitioner must be ignored
  const auto mesh = generate(0, false, "metis");

  const auto parent_sides = sidePairs(*mesh, MooseMeshUtils::getBoundaryID("boundary", *mesh));
  ASSERT_EQ(parent_sides.size(), 4);

  const auto patch_sides = collectPatchSides(*mesh, "boundary_", 4);
  ASSERT_EQ(patch_sides.size(), 4);
  for (const auto j : make_range(4))
    EXPECT_EQ(sidePairs(*mesh, MooseMeshUtils::getBoundaryID("boundary_" + std::to_string(j),
                                                             *mesh))
                  .size(),
              1);
  EXPECT_EQ(patch_sides, parent_sides);
}

TEST_F(PatchSidesetGeneratorTest, singlePatchKeepsWholeSideset)
{
  const auto mesh = generate(1, false);

  const auto parent_sides = sidePairs(*mesh, MooseMeshUtils::getBoundaryID("boundary", *mesh));
  ASSERT_EQ(parent_sides.size(), 4);

  const BoundaryName patch_name = "boundary_0";
  ASSERT_TRUE(MooseMeshUtils::hasBoundaryNameOrID(*mesh, patch_name));
  EXPECT_EQ(sidePairs(*mesh, MooseMeshUtils::getBoundaryID(patch_name, *mesh)), parent_sides);
}
