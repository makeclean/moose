//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ViewFactorMatrix.h"
#include "MooseUtils.h"

#include "libmesh/parallel.h"

#include <algorithm>
#include <map>

ViewFactorMatrix::ViewFactorMatrix(unsigned int n_rows) : _rows(n_rows) {}

unsigned int
ViewFactorMatrix::nnz() const
{
  unsigned int n = 0;
  for (const auto & row : _rows)
    n += row.size();
  return n;
}

Real
ViewFactorMatrix::get(unsigned int row, unsigned int col) const
{
  const auto & r = _rows[row];
  const auto it = std::lower_bound(
      r.begin(), r.end(), col, [](const Entry & e, unsigned int c) { return e.first < c; });
  if (it != r.end() && it->first == col)
    return it->second;
  return 0;
}

void
ViewFactorMatrix::add(unsigned int row, unsigned int col, Real delta)
{
  auto & r = _rows[row];
  const auto it = std::lower_bound(
      r.begin(), r.end(), col, [](const Entry & e, unsigned int c) { return e.first < c; });
  if (it != r.end() && it->first == col)
    it->second += delta;
  else
    r.insert(it, {col, delta});
}

void
ViewFactorMatrix::merge(const ViewFactorMatrix & other)
{
  for (const auto i : make_range(other.nRows()))
    for (const auto & [col, value] : other.row(i))
      add(i, col, value);
}

void
ViewFactorMatrix::clear()
{
  for (auto & row : _rows)
    row.clear();
}

void
ViewFactorMatrix::gatherSum(const libMesh::Parallel::Communicator & comm)
{
  if (comm.size() <= 1)
    return;

  // Pack all nonzeros into a single buffer per processor to avoid an O(N) collective storm
  std::vector<unsigned int> flat_indices;
  std::vector<Real> flat_values;
  for (const auto i : make_range(nRows()))
    for (const auto & [col, value] : row(i))
    {
      flat_indices.push_back(i);
      flat_indices.push_back(col);
      flat_values.push_back(value);
    }

  std::vector<std::vector<unsigned int>> all_indices;
  std::vector<std::vector<Real>> all_values;
  comm.allgather(flat_indices, all_indices);
  comm.allgather(flat_values, all_values);

  clear();
  for (const auto p : index_range(all_indices))
  {
    const auto & indices = all_indices[p];
    const auto & values = all_values[p];
    for (const auto k : index_range(values))
      add(indices[2 * k], indices[2 * k + 1], values[k]);
  }
}
