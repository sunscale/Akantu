/**
 * @file   test_parser.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Nov 13 2013
 * @date last modification: Sun Jul 09 2017
 *
 * @brief  test the input file parser
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_random_generator.hh"
#include "parser.hh"

#include <iostream>

using namespace akantu;

int main(int argc, char * argv[]) {
  initialize("input_file.dat", argc, argv);

  const Parser & p = getStaticParser();

  std::cout << RandomGenerator<UInt>::seed() << "==123456" << std::endl;

  std::cout << p << std::endl;

  Real toto = p.getParameter("toto");
  std::cout << toto;
  Real ref = 2 * M_PI + std::max(2., 50.);
  if (std::abs(toto - ref) > std::numeric_limits<Real>::epsilon()) {
    std::cout << "!=" << ref << std::endl;
    return 1;
  }

  std::cout << "==" << ref << std::endl;

  Vector<Real> vect = p.getParameter("vect");
  std::cout << vect << std::endl;

  Matrix<Real> mat = p.getParameter("mat");
  std::cout << mat << std::endl;

  RandomParameter<Real> rand1 = p.getParameter("rand1");
  std::cout << rand1 << std::endl;

  RandomParameter<Real> rand2 = p.getParameter("rand2");
  std::cout << rand2 << std::endl;

  RandomParameter<Real> rand3 = p.getParameter("rand3");
  std::cout << rand3 << std::endl;

  finalize();
  return 0;
}
