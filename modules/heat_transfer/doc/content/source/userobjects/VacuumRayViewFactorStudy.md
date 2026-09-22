# VacuumRayViewFactorStudy

!syntax description /UserObjects/VacuumRayViewFactorStudy

## Description

`VacuumRayViewFactorStudy` fires rays into the space surrounding a set of surfaces to
compute the view factors exchanged between those surfaces and the fraction of
radiation that leaves an open (not fully enclosed) geometry. It is the analog of
[ViewFactorRayStudy.md] for geometries that do not form a closed cavity: the
surfaces radiate into a surrounding void (a vacuum) rather than into an enclosed
cavity filled with a transparent medium, so part of the emitted radiation escapes
the geometry without being intercepted by any other surface.

For MOOSE's ray tracing to work, the void into which the surfaces radiate must be
resolvable by the raytracer. `VacuumRayViewFactorStudy` uses the XDG library for
this purpose: the void is represented as the implicit complement of a surface mesh
and the raytracer sends rays through it.

The study is created automatically by the [GrayDiffuseRadiation](syntax/GrayDiffuseRadiation/index.md)
action when the [!param](/GrayDiffuseRadiation/view_factor_calculator) is set to `vacuum_ray_tracing`.
The resulting view factors and escape fractions are consumed by [VacuumRayViewFactor.md].

## Theory

The view factors are computed by discretizing the angular integral over the
hemisphere above each surface. The angular quadrature is the same half-range
Gauss-Legendre-Chebyshev quadrature adopted from [!citep](WaltersLCQ) that
[ViewFactorRayStudy.md] uses. The polar angle is measured with respect to the
outward normal and is sampled with [!param](/UserObjects/VacuumRayViewFactorStudy/polar_quad_order)
Gauss-Legendre nodes on the interval $(0, 1)$, and the azimuthal angle is sampled
with $4 \times$ [!param](/UserObjects/VacuumRayViewFactorStudy/azimuthal_quad_order)
Chebyshev nodes per hemisphere.

Each surface is spatially sampled with a quadrature controlled by the
[!param](/UserObjects/VacuumRayViewFactorStudy/face_type) and
[!param](/UserObjects/VacuumRayViewFactorStudy/face_order) parameters. The default
`GRID` / `CONSTANT` choice places one quadrature point at the centroid of each
element face.

From each quadrature point, rays are fired along the outward normal directions of
the angular quadrature. The rays traverse the void and are terminated at their
first intersection with another surface, or when they leave the geometry. A ray
that hits surface $j$ increments the view factor weight between the originating
surface and surface $j$; a ray that leaves the geometry increments the escape
weight of the originating surface. Because the rays are fired along the outward
normal, surfaces of a convex geometry never intercept rays from another surface and
their escape fraction is exactly one.

## Important Conventions

- The calculation is three-dimensional and requires a serial replicated mesh.
- Every surface that borders the void must be listed in the
  [!param](/UserObjects/VacuumRayViewFactorStudy/boundary) parameter; an error is
  raised otherwise, because radiation reaching an unlisted surface would otherwise
  be silently lost.

## Example Input Syntax

In this example, `VacuumRayViewFactorStudy` is set up implicitly by the
[GrayDiffuseRadiation](syntax/GrayDiffuseRadiation/index.md) action for an open
(hollow) box whose walls exchange radiation through a void:

!listing modules/heat_transfer/test/tests/vacuum_radiator/open_box.i
block=GrayDiffuseRadiation

!syntax parameters /UserObjects/VacuumRayViewFactorStudy

!syntax inputs /UserObjects/VacuumRayViewFactorStudy

!syntax children /UserObjects/VacuumRayViewFactorStudy

!bibtex bibliography