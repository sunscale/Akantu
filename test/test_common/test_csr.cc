/**
 * @file   test_csr.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Jul 30 2012
 * @date last modification: Sun Dec 03 2017
 *
 * @brief  Test the CSR (compressed sparse row) data structure
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_csr.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
/* -------------------------------------------------------------------------- */

using namespace akantu;

class TestCsrFixture : public ::testing::Test {
protected:
  void SetUp() override {
    csr.resizeRows(N);
    csr.clearRows();

    for (UInt i = 0; i < N; ++i) {
      UInt nb_cols(UInt(rand() * double(N) / (RAND_MAX + 1.)));
      nb_cols_per_row.push_back(nb_cols);
      for (UInt j = 0; j < nb_cols; ++j) {
        ++csr.rowOffset(i);
      }
    }

    csr.countToCSR();
    csr.resizeCols();

    csr.beginInsertions();
    for (UInt i = 0; i < N; ++i) {
      UInt nb_cols = nb_cols_per_row[i];
      for (UInt j = 0; j < nb_cols; ++j) {
        csr.insertInRow(i, nb_cols - j);
      }
    }
    csr.endInsertions();
  }

  std::vector<UInt> nb_cols_per_row;
  CSR<UInt> csr;
  size_t N = 1000;
};

TEST_F(TestCsrFixture, CheckInsertion) { EXPECT_EQ(N, this->csr.getNbRows()); }

TEST_F(TestCsrFixture, Iteration) {
  for (UInt i = 0; i < this->csr.getNbRows(); ++i) {
    auto it = this->csr.begin(i);
    auto end = this->csr.end(i);
    UInt nb_cols = this->nb_cols_per_row[i];
    for (; it != end; ++it) {
      EXPECT_EQ(nb_cols, *it);
      nb_cols--;
    }

    EXPECT_EQ(0, nb_cols);
  }
}

TEST_F(TestCsrFixture, ReverseIteration) {
  for (UInt i = 0; i < csr.getNbRows(); ++i) {
    auto it = csr.rbegin(i);
    auto end = csr.rend(i);

    UInt nb_cols = nb_cols_per_row[i];
    UInt j = nb_cols;

    for (; it != end; --it) {
      EXPECT_EQ((nb_cols - j + 1), *it);
      j--;
    }

    EXPECT_EQ(0, j);
  }
}
