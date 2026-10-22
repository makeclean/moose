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

#include "timpi/communicator.h"

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
  // Every row participates in the collectives, even rows that are empty on
  // this processor, so that all processors call the same collective sequence
  for (const auto i : make_range(nRows()))
  {
    std::vector<unsigned int> cols;
    std::vector<Real> values;
    for (const auto & [col, value] : row(i))
    {
      cols.push_back(col);
      values.push_back(value);
    }

    // the columns and values of this row on every processor; the sum over
    // the processors is then a local computation
    std::vector<std::vector<unsigned int>> all_cols;
    std::vector<std::vector<Real>> all_values;
    comm.allgather(cols, all_cols);
    comm.allgather(values, all_values);

    // std::map keeps the union of the columns sorted
    std::map<unsigned int, Real> sums;
    for (const auto p : make_range(all_cols.size()))
      for (const auto k : make_range(all_cols[p].size()))
        sums[all_cols[p][k]] += all_values[p][k];

    std::vector<Entry> new_row;
    new_row.reserve(sums.size());
    for (const auto & [col, sum] : sums)
      if (sum != 0)
        new_row.emplace_back(col, sum);
    row(i) = std::move(new_row);
  }
}
