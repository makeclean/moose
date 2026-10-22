//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ViewFactorMatrix.h"
#include "timpi/communicator.h"
#include "gtest/gtest.h"

#include <utility>
#include <vector>

TEST(ViewFactorMatrix, getReturnsZeroForMissingEntry)
{
  ViewFactorMatrix m(3);
  EXPECT_DOUBLE_EQ(m.get(0, 0), 0.0);
  m.add(1, 2, 0.5);
  EXPECT_DOUBLE_EQ(m.get(1, 2), 0.5);
  EXPECT_DOUBLE_EQ(m.get(1, 0), 0.0);
  EXPECT_DOUBLE_EQ(m.get(0, 2), 0.0);
}

TEST(ViewFactorMatrix, addAccumulatesAndKeepsRowsSorted)
{
  ViewFactorMatrix m(2);
  m.add(0, 3, 1.0);
  m.add(0, 1, 2.0);
  m.add(0, 2, 4.0);
  m.add(0, 1, 0.5); // accumulates onto the existing entry

  ASSERT_EQ(m.nnz(), 3);
  const auto & row = m.row(0);
  ASSERT_EQ(row.size(), 3);
  EXPECT_EQ(row[0].first, 1u);
  EXPECT_DOUBLE_EQ(row[0].second, 2.5);
  EXPECT_EQ(row[1].first, 2u);
  EXPECT_DOUBLE_EQ(row[1].second, 4.0);
  EXPECT_EQ(row[2].first, 3u);
  EXPECT_DOUBLE_EQ(row[2].second, 1.0);
}

TEST(ViewFactorMatrix, mergeSumsOverlappingColumns)
{
  ViewFactorMatrix a(3);
  a.add(0, 0, 1.0);
  a.add(0, 2, 3.0);
  a.add(1, 1, 2.0);

  ViewFactorMatrix b(3);
  b.add(0, 2, 1.0); // sums with the existing entry
  b.add(0, 1, 5.0); // new column, must be inserted sorted
  b.add(2, 0, 7.0);

  a.merge(b);
  ASSERT_EQ(a.nnz(), 5);
  const auto & row0 = a.row(0);
  ASSERT_EQ(row0.size(), 3);
  EXPECT_EQ(row0[0].first, 0u);
  EXPECT_DOUBLE_EQ(row0[0].second, 1.0);
  EXPECT_EQ(row0[1].first, 1u);
  EXPECT_DOUBLE_EQ(row0[1].second, 5.0);
  EXPECT_EQ(row0[2].first, 2u);
  EXPECT_DOUBLE_EQ(row0[2].second, 4.0);
  EXPECT_DOUBLE_EQ(a.get(1, 1), 2.0);
  EXPECT_DOUBLE_EQ(a.get(2, 0), 7.0);
}

TEST(ViewFactorMatrix, clearRemovesEntriesButKeepsRows)
{
  ViewFactorMatrix m(4);
  m.add(0, 0, 1.0);
  m.add(3, 1, 2.0);
  EXPECT_EQ(m.nnz(), 2u);

  m.clear();
  EXPECT_EQ(m.nRows(), 4u);
  EXPECT_EQ(m.nnz(), 0u);

  // the matrix is usable after a clear
  m.add(2, 3, 1.5);
  EXPECT_DOUBLE_EQ(m.get(2, 3), 1.5);
}

TEST(ViewFactorMatrix, gatherSumSumsOverProcessors)
{
  // the default-constructed communicator wraps the MPI world when MPI is
  // initialized, so this test is meaningful for any process count
  const libMesh::Parallel::Communicator comm;
  const unsigned int rank = comm.rank();
  const unsigned int size = comm.size();

  ViewFactorMatrix m(3);
  // every processor contributes (rank + 1) to (0, 1)
  m.add(0, 1, static_cast<libMesh::Real>(rank) + 1.0);
  // rank 0 contributes +1 and rank 1 contributes -1 to (1, 0); the entry
  // cancels to zero when both are present and is dropped
  if (rank == 0)
    m.add(1, 0, 1.0);
  if (size > 1 && rank == 1)
    m.add(1, 0, -1.0);
  // a local-only entry
  m.add(2, 2, static_cast<libMesh::Real>(rank) + 0.5);

  m.gatherSum(comm);

  libMesh::Real expected = 0.0;
  for (unsigned int p = 0; p < size; ++p)
    expected += static_cast<libMesh::Real>(p) + 1.0;
  EXPECT_DOUBLE_EQ(m.get(0, 1), expected);

  if (size > 1)
  {
    EXPECT_DOUBLE_EQ(m.get(1, 0), 0.0);
    EXPECT_TRUE(m.row(1).empty());
  }
  else
    EXPECT_DOUBLE_EQ(m.get(1, 0), 1.0);

  expected = 0.0;
  for (unsigned int p = 0; p < size; ++p)
    expected += static_cast<libMesh::Real>(p) + 0.5;
  EXPECT_DOUBLE_EQ(m.get(2, 2), expected);
}
