# RadiationTransferAction

This action corresponds to the
[GrayDiffuseRadiation](/GrayDiffuseRadiation/index.md) system; see its documentation page for more information.

With [!param](/GrayDiffuseRadiation/create_patch_sidesets) set to `true`, each
entry of [!param](/GrayDiffuseRadiation/n_patches) subdivides the corresponding
[!param](/GrayDiffuseRadiation/boundary) sideset into that many radiation patches by appending a
[PatchSidesetGenerator](/Mesh/PatchSidesetGenerator) to the mesh generator tree. A value of `1`
keeps a single patch per sideset, and a value of `0` creates one patch for every face of the
sideset. When a sideset is subdivided in the vacuum calculation, the original sideset is removed
from the mesh and replaced by its patch sidesets, so subdivided boundaries must be given by name.

By default [!param](/GrayDiffuseRadiation/create_patch_sidesets) is `false`: the patches are then
indexed internally instead of being added to the mesh, which avoids exhausting libMesh's boundary
id space when a sideset is split into many patches. Internal patches are only supported with
[!param](/GrayDiffuseRadiation/view_factor_calculator) = `vacuum_ray_tracing`; the other
calculators require `create_patch_sidesets = true` when
[!param](/GrayDiffuseRadiation/n_patches) is not one.
