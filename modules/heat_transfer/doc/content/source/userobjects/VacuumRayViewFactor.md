# VacuumRayViewFactor

!syntax description /UserObjects/VacuumRayViewFactor

## Description

`VacuumRayViewFactor` computes the view factors $F_{i,j}$ exchanged between surfaces
radiating into an open (vacuum) environment and the escape fraction $e_i$ of each
surface, which accounts for the radiation that leaves the geometry without being
intercepted by another surface. It is the analog of [RayTracingViewFactor.md] for
geometries that are not fully enclosed, and it relies on a
[VacuumRayViewFactorStudy.md] to perform the ray tracing.

The physical setup is the same as in the net radiation method for
[gray, diffuse radiative exchange](modules/heat_transfer/index.md#gray_diffuse_radiative_exchange),
except that the surfaces no longer completely enclose the medium: part of the
radiation emitted by a surface leaves the geometry so that the row sum of the view
factors no longer equals one,

\begin{equation}
  \sum_j F_{i,j} = 1 - e_i,
\end{equation}

where $e_i$ is the escape fraction of surface $i$. The irradiation onto surface $i$
then consists of the radiation arriving from the other surfaces plus the radiation
arriving from the environment,

\begin{equation}
  H_i = \sum_j F_{i,j} J_j + e_i B_{\text{env}},
\end{equation}

where $B_{\text{env}}$ is the black body emission of the environment, which is zero
for a vacuum and $\sigma T_{\text{env}}^4$ for a black body environment at
temperature $T_{\text{env}}$.

`VacuumRayViewFactor` also validates the computed view factors: the deviation of
$\sum_j F_{i,j} + e_i$ from one is checked against the
[!param](/UserObjects/VacuumRayViewFactor/view_factor_tol) tolerance. Because the
escape fraction must not be redistributed among the participating surfaces, the
normalization of view factors is not supported by this object and the
[!param](/UserObjects/VacuumRayViewFactor/normalize_view_factor) parameter is
forced to `false`.

The escape fraction is accessible through the
[!param](/UserObjects/VacuumRayViewFactor/ray_study_name) study's public interface
and is used by [GrayLambertSurfaceRadiationBase.md] to close the energy balance of
open geometries.

## Important Conventions

- The calculation is three-dimensional and requires a serial replicated mesh.
- Every surface that borders the void must be a participating surface, i.e. be
  listed in the [!param](/UserObjects/VacuumRayViewFactor/boundary) parameter.

## Example Input Syntax

In this example, `VacuumRayViewFactor` is set up implicitly by the
[GrayDiffuseRadiation](syntax/GrayDiffuseRadiation/index.md) action for a convex
surface radiating into a vacuum:

!listing modules/heat_transfer/test/tests/vacuum_radiator/convex_escape.i
block=GrayDiffuseRadiation

!syntax parameters /UserObjects/VacuumRayViewFactor

!syntax inputs /UserObjects/VacuumRayViewFactor

!syntax children /UserObjects/VacuumRayViewFactor