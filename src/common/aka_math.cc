/**
 * @file   aka_math.cc
 *
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 * @author Peter Spijker <peter.spijker@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Wed Aug 04 2010
 * @date last modification: Thu Mar 27 2014
 *
 * @brief  Implementation of the math toolbox
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_array.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
void Math::matrix_vector(UInt m, UInt n,
			 const Array<Real> & A,
			 const Array<Real> & x,
			 Array<Real> & y,
			 Real alpha) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(A.getSize() == x.getSize(),
		      "The vector A(" << A.getID()
		      << ") and the vector x(" << x.getID()
		      << ") must have the same size");

  AKANTU_DEBUG_ASSERT(A.getNbComponent() == m * n,
		      "The vector A(" << A.getID()
		      << ") has the good number of component.");

  AKANTU_DEBUG_ASSERT(x.getNbComponent() == n,
		      "The vector x(" << x.getID()
		      << ") do not the good number of component.");

  AKANTU_DEBUG_ASSERT(y.getNbComponent() == n,
		      "The vector y(" << y.getID()
		      << ") do not the good number of component.");


  UInt nb_element = A.getSize();
  UInt offset_A = A.getNbComponent();
  UInt offset_x = x.getNbComponent();

  y.resize(nb_element);

  Real * A_val = A.storage();
  Real * x_val = x.storage();
  Real * y_val = y.storage();

  for (UInt el = 0; el < nb_element; ++el) {
    matrix_vector(m, n, A_val, x_val, y_val, alpha);

    A_val += offset_A;
    x_val += offset_x;
    y_val += offset_x;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Math::matrix_matrix(UInt m, UInt n, UInt k,
			 const Array<Real> & A,
			 const Array<Real> & B,
			 Array<Real> & C,
			 Real alpha) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(A.getSize() == B.getSize(),
		      "The vector A(" << A.getID()
		      << ") and the vector B(" << B.getID()
		      << ") must have the same size");

  AKANTU_DEBUG_ASSERT(A.getNbComponent() == m * k,
		      "The vector A(" << A.getID()
		      << ") has the good number of component.");

  AKANTU_DEBUG_ASSERT(B.getNbComponent() == k * n ,
		      "The vector B(" << B.getID()
		      << ") do not the good number of component.");

  AKANTU_DEBUG_ASSERT(C.getNbComponent() == m * n,
		      "The vector C(" << C.getID()
		      << ") do not the good number of component.");

  UInt nb_element = A.getSize();
  UInt offset_A = A.getNbComponent();
  UInt offset_B = B.getNbComponent();
  UInt offset_C = C.getNbComponent();

  C.resize(nb_element);

  Real * A_val = A.storage();
  Real * B_val = B.storage();
  Real * C_val = C.storage();

  for (UInt el = 0; el < nb_element; ++el) {
    matrix_matrix(m, n, k, A_val, B_val, C_val, alpha);

    A_val += offset_A;
    B_val += offset_B;
    C_val += offset_C;
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void Math::matrix_matrixt(UInt m, UInt n, UInt k,
			  const Array<Real> & A,
			  const Array<Real> & B,
			  Array<Real> & C,
			  Real alpha) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(A.getSize() == B.getSize(),
		      "The vector A(" << A.getID()
		      << ") and the vector B(" << B.getID()
		      << ") must have the same size");

  AKANTU_DEBUG_ASSERT(A.getNbComponent() == m * k,
		      "The vector A(" << A.getID()
		      << ") has the good number of component.");

  AKANTU_DEBUG_ASSERT(B.getNbComponent() == k * n ,
		      "The vector B(" << B.getID()
		      << ") do not the good number of component.");

  AKANTU_DEBUG_ASSERT(C.getNbComponent() == m * n,
		      "The vector C(" << C.getID()
		      << ") do not the good number of component.");

  UInt nb_element = A.getSize();
  UInt offset_A = A.getNbComponent();
  UInt offset_B = B.getNbComponent();
  UInt offset_C = C.getNbComponent();

  C.resize(nb_element);

  Real * A_val = A.storage();
  Real * B_val = B.storage();
  Real * C_val = C.storage();

  for (UInt el = 0; el < nb_element; ++el) {
    matrix_matrixt(m, n, k, A_val, B_val, C_val, alpha);

    A_val += offset_A;
    B_val += offset_B;
    C_val += offset_C;
  }

  AKANTU_DEBUG_OUT();
}



__END_AKANTU__
