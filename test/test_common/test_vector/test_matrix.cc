/**
 * @file   test_matrix.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Mar 03 2011
 * @date last modification: Tue Aug 19 2014
 *
 * @brief  tests for Matrix
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
#include "aka_common.hh"
#include "aka_types.hh"
#include "aka_array.hh"
#include "aka_math.hh"
/* -------------------------------------------------------------------------- */
#include <cstdlib>
#include <sys/time.h>
#include <iostream>



/* -------------------------------------------------------------------------- */
using namespace akantu;


int main(int argc, char *argv[]) {

#define n 4
  UInt nbm = 10000;

  Real time;


  Matrix<Real> eig_test(3, 3);
  eig_test(0, 0) = 0; eig_test(0, 1) = 1; eig_test(0, 2) = 0;
  eig_test(1, 0) = 1; eig_test(1, 1) = 0; eig_test(1, 2) = 0;
  eig_test(2, 0) = 0; eig_test(2, 1) = 0; eig_test(2, 2) = 1;

  Matrix<Real> eig_vects(3, 3);
  Vector<Real> eigs(3);

  eig_test.eig(eigs, eig_vects);

  std::cout << "A: " << eig_test << std::endl;
  std::cout << "eig(A): " << eigs << " " << Vector<Real>(eig_vects(0)) << " " << Vector<Real>(eig_vects(1)) << " " << Vector<Real>(eig_vects(2)) << std::endl;

  Matrix<Real> check_eig_vects(3, 3);
  check_eig_vects(0, 0) = 1./std::sqrt(2); check_eig_vects(0, 1) = 0; check_eig_vects(0, 2) = -1./std::sqrt(2);
  check_eig_vects(1, 0) = 1./std::sqrt(2); check_eig_vects(1, 1) = 0; check_eig_vects(1, 2) = 1./std::sqrt(2);
  check_eig_vects(2, 0) = 0;               check_eig_vects(2, 1) = 1; check_eig_vects(2, 2) = 0;

  for (UInt i = 0; i < 3; ++i) {
    Vector<Real> eig_v = eig_vects(i);
    Vector<Real> check_eig_v = check_eig_vects(i);

    if(! check_eig_v.equal(eig_v, 1e-14)) {
      AKANTU_DEBUG_ERROR("The " << i << "th eigen vector is not correct: " << eig_v << " should be " << check_eig_v);
    }
  }



  Array<Real> A(nbm, n*n);
  Array<Real> B(nbm, n*n);
  Array<Real> C1(nbm, n*n);
  Array<Real> C2(nbm, n*n);
  Array<Real> C3(nbm, n*n);
  Array<Real> C4(nbm, n*n);

  for (UInt i = 0; i < n*n; ++i) {
    A.storage()[i] = drand48();
    B.storage()[i] = drand48();
  }

  for (UInt i = 1; i < nbm; ++i) {
    memcpy(A.storage() + i * n * n, A.storage(), n*n*sizeof(Real));
    memcpy(B.storage() + i * n * n, B.storage(), n*n*sizeof(Real));
  }

  struct timeval begin, end;

  Array<Real>::matrix_iterator itA = A.begin(n,n);
  Array<Real>::matrix_iterator itB = B.begin(n,n);
  itA = A.begin(n,n);
  itB = B.begin(n,n);
  std::cerr << *itA << std::endl;
  std::cerr << *itB << std::endl;

  /* ------------------------------------------------------------------------ */
  gettimeofday(&begin, NULL);
  Math::matrix_matrix(n, n, n, A, B, C1);
  gettimeofday(&end, NULL);

  //time =  (end.tv_sec * 1e3 + end.tv_usec * 1e-3) - (begin.tv_sec * 1e3 + begin.tv_usec * 1e-3);
  time =  (end.tv_sec * 1e6 + end.tv_usec) - (begin.tv_sec * 1e6 + begin.tv_usec);
  std::cout << "matrix_matrix : " << std::fixed << time/nbm << "us" << std::endl;

  /* ------------------------------------------------------------------------ */
  Array<Real>::matrix_iterator itC = C2.begin(n,n);

  gettimeofday(&begin, NULL);
  for (UInt i = 0; i < nbm; ++i) {
    Matrix<Real> & C = *itC;
    C = *itA * *itB;
    ++itA; ++itB;++itC;
  }
  gettimeofday(&end, NULL);

  itC = C2.begin(n,n);
  std::cerr << *itC << std::endl;
  time =  (end.tv_sec * 1e6 + end.tv_usec) - (begin.tv_sec * 1e6 + begin.tv_usec);
  std::cout << "it Mc() = it Ma() * it Mb() : " << std::fixed << time/nbm << "us" << std::endl;

  /* ------------------------------------------------------------------------ */
  Array<Real>::matrix_iterator muitA = A.begin(n,n);
  Array<Real>::matrix_iterator muitB = B.begin(n,n);
  Array<Real>::matrix_iterator muitC = C4.begin(n,n);
  gettimeofday(&begin, NULL);
  for (UInt i = 0; i < nbm; ++i) {
    muitC->mul<false, false>(*muitA, *muitB);
    ++muitA; ++muitB;++muitC;
  }
  gettimeofday(&end, NULL);

  muitC = C4.begin(n,n);
  std::cerr << *muitC << std::endl;

  time =  (end.tv_sec * 1e6 + end.tv_usec) - (begin.tv_sec * 1e6 + begin.tv_usec);
  std::cout << "it Mc.mul( it Ma(),  it Mb()) : " << std::fixed << time/nbm << "us" << std::endl;

  return 0;
}
