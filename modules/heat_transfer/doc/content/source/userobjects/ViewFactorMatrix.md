# ViewFactorMatrix

`ViewFactorMatrix` is the sparse storage used by the view factor user objects
([ViewFactorBase](/ViewFactorBase/index.md) and its derivatives) and by
[GrayLambertSurfaceRadiationBase](/GrayLambertSurfaceRadiationBase/index.md).

## Description

View factor matrices of patched (per-facet) boundaries are very sparse: a facet
only sees the facets within its line of sight, so the number of nonzero entries
grows roughly linearly with the number of patches, while a dense storage would
grow quadratically. `ViewFactorMatrix` stores each row as a vector of
(column, value) pairs sorted by column; an entry that is not stored is zero.

The class provides the operations the view factor user objects need:

- `add(row, col, delta)`: adds `delta` to the entry, creating it if absent
  (used while accumulating ray weights or quadrature contributions),
- `merge(other)`: adds all entries of another matrix (thread join),
- `gatherSum(comm)`: sums the entries over the MPI communicator so that every
  processor ends with the identical global matrix; entries whose sum is zero
  are dropped. Every row participates in the collectives, even rows that are
  empty on a given processor,
- `get(row, col)`: returns the value of an entry, or zero if it is not stored,
- `row(i)`: the stored (column, value) pairs of row `i`, sorted by column.

The radiosity solve in
[GrayLambertSurfaceRadiationBase](/GrayLambertSurfaceRadiationBase/index.md)
stamps the matrix entries into a PETSc `Mat` (the pattern is the nonzero
structure of the view factors plus the diagonal) and solves with GMRES and an
ILU(0) preconditioner instead of a dense LU solve, which does not scale to
tens of thousands of patches.
