//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "libmesh/libmesh_common.h"

#include <utility>
#include <vector>

namespace libMesh::Parallel
{
class Communicator;
}

/**
 * Sparse storage of the square view factor matrix.
 *
 * Each row holds its stored (column, value) pairs sorted by column; entries
 * that are not stored are zero. View factor matrices of patched (per-facet)
 * boundaries are very sparse, because a facet only sees the facets that lie
 * within its line of sight: the number of stored entries grows roughly
 * linearly with the number of patches, while a dense storage would grow
 * quadratically.
 */
class ViewFactorMatrix
{
public:
  /// A stored entry: (column, value)
  using Entry = std::pair<unsigned int, libMesh::Real>;

  ViewFactorMatrix() = default;

  explicit ViewFactorMatrix(unsigned int n_rows);

  /// Number of rows (and columns)
  unsigned int nRows() const { return _rows.size(); }

  /// Total number of stored entries
  unsigned int nnz() const;

  /// The stored (column, value) pairs of row i, sorted by column
  const std::vector<Entry> & row(unsigned int i) const { return _rows[i]; }
  std::vector<Entry> & row(unsigned int i) { return _rows[i]; }

  /// The value at (row, col); zero if the entry is not stored
  libMesh::Real get(unsigned int row, unsigned int col) const;

  /// Adds delta to the entry (row, col), creating the entry if it is absent
  void add(unsigned int row, unsigned int col, libMesh::Real delta);

  /// Adds all entries of other to this matrix
  void merge(const ViewFactorMatrix & other);

  /// Removes all entries; the number of rows is unchanged
  void clear();

  /// Sums this matrix over the communicator; entries whose sum is zero are
  /// dropped. Collective: must be called by every processor.
  void gatherSum(const libMesh::Parallel::Communicator & comm);

private:
  std::vector<std::vector<Entry>> _rows;
};
