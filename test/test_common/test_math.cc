/**
 * @file   test_math.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Feb 21 2014
 * @date last modification: Fri Feb 28 2014
 *
 * @brief  test the static Math class
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_math.hh"

using namespace akantu;

void checkVect(Real * x, Real * xref, UInt n) {
  Real diff[n];
  for(UInt i = 0; i < n; ++i) {
    diff[i] = xref[i]- x[i];
  }

  Real norm = Math::norm(n, diff)/Math::norm(n, xref);
  Real tol = 1e-12;
  AKANTU_DEBUG_ASSERT(norm < tol, "x differs form xref");
}

/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[]) {
  Real A[3*5] = {  0,   1,  2,
		   3,   4,  5,
		   6,   7,  8,
		   9,  10, 11,
		   12, 13, 14, };

  Real x1[5] = { 0, 1, 2, 3, 4 };
  Real x2[3] = { 0, 1, 2 };

  Real y1[3];
  Math::matrix_vector(3, 5, A, x1, y1, 1.);
  Real y1_ref[3] = { 90, 100, 110 };

  checkVect(y1, y1_ref, 3);

  Real y2[5];
  Math::matVectMul<true>(3, 5, 1., A, x2, 0., y2);
  Real y2_ref[5] = { 5, 14, 23, 32, 41 };

  checkVect(y2, y2_ref, 5);

  return 0;
}
