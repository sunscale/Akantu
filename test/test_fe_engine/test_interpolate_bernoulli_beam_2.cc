/**
 * @file   test_interpolate_bernoulli_beam_2.cc
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 *
 * @date creation: Fri Jul 15 2011
 * @date last modification: Sat Jan 23 2016
 *
 * @brief  Test of the interpolation on the type _bernoulli_beam_2
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
#include <cstdlib>
#include <fstream>
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "fe_engine.hh"
#include "fe_engine_template.hh"
#include "integrator_gauss.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "shape_linked.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main() {

  Mesh beams(2);

  /* --------------------------------------------------------------------------
   */
  // Defining the mesh

  Array<Real> & nodes = const_cast<Array<Real> &>(beams.getNodes());
  nodes.resize(4);

  beams.addConnectivityType(_bernoulli_beam_2);
  Array<UInt> & connectivity =
      const_cast<Array<UInt> &>(beams.getConnectivity(_bernoulli_beam_2));
  connectivity.resize(3);

  for (UInt i = 0; i < 4; ++i) {

    nodes(i, 0) = (i + 1) * 2;
    nodes(i, 1) = 1;
  }
  for (UInt i = 0; i < 3; ++i) {

    connectivity(i, 0) = i;
    connectivity(i, 1) = i + 1;
  }
  akantu::MeshIOMSH mesh_io;
  mesh_io.write("b_beam_2.msh", beams);

  /* --------------------------------------------------------------------------
   */
  // Interpolation

  FEEngineTemplate<IntegratorGauss, ShapeLinked> * fem =
      new FEEngineTemplate<IntegratorGauss, ShapeLinked>(beams, 2);

  fem->initShapeFunctions();

  Array<Real> displ_on_nodes(4, 3);
  Array<Real> displ_on_quad(0, 3);

  for (UInt i = 0; i < 4; ++i) {

    displ_on_nodes(i, 0) = (i + 1) * 2; // Definition of the displacement
    displ_on_nodes(i, 1) = 0;
    displ_on_nodes(i, 2) = 0;
  }

  fem->getShapeFunctions().interpolateOnControlPoints<_bernoulli_beam_2>(
      displ_on_nodes, displ_on_quad, 3, _not_ghost, NULL, false, 0, 0, 0);

  fem->getShapeFunctions().interpolateOnControlPoints<_bernoulli_beam_2>(
      displ_on_nodes, displ_on_quad, 3, _not_ghost, NULL, false, 1, 1, 1);

  fem->getShapeFunctions().interpolateOnControlPoints<_bernoulli_beam_2>(
      displ_on_nodes, displ_on_quad, 3, _not_ghost, NULL, true, 2, 2, 1);

  fem->getShapeFunctions().interpolateOnControlPoints<_bernoulli_beam_2>(
      displ_on_nodes, displ_on_quad, 3, _not_ghost, NULL, false, 3, 2, 3);

  fem->getShapeFunctions().interpolateOnControlPoints<_bernoulli_beam_2>(
      displ_on_nodes, displ_on_quad, 3, _not_ghost, NULL, true, 4, 3, 3);

  Real * don = displ_on_nodes.storage();
  Real * doq = displ_on_quad.storage();

  std::ofstream my_file("out.txt");
  my_file << don << std::endl;
  my_file << doq << std::endl;

  return EXIT_SUCCESS;
}
