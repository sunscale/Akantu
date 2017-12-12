/**
 * @file   test_arange_iterator.cc
 *
 * @author Nicolas Richart
 *
 * @date creation  Fri Aug 11 2017
 *
 * @brief A Documented file.
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */
/* -------------------------------------------------------------------------- */
#include "aka_iterators.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
/* -------------------------------------------------------------------------- */

using namespace akantu;

TEST(TestArangeIterator, Stop) {
  size_t ref_i = 0;
  for (auto i : arange(10)) {
    EXPECT_EQ(ref_i, i);
    ++ref_i;
  }
}

TEST(TestArangeIterator, StartStop) {
  size_t ref_i = 1;
  for (auto i : arange(1, 10)) {
    EXPECT_EQ(ref_i, i);
    ++ref_i;
  }
}

TEST(TestArangeIterator, StartStopStep) {
  size_t ref_i = 1;
  for (auto i : arange(1, 22, 2)) {
    EXPECT_EQ(ref_i, i);
    ref_i += 2;
  }
}

TEST(TestArangeIterator, StartStopStepZipped) {
  int ref_i1 = -1, ref_i2 = 1;
  for (auto && i : zip(arange(-1, -10, -1), arange(1, 18, 2))) {
    EXPECT_EQ(ref_i1, std::get<0>(i));
    EXPECT_EQ(ref_i2, std::get<1>(i));
    ref_i1 += -1;
    ref_i2 += 2;
  }
}
