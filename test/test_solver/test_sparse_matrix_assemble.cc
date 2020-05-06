/**
 * @file   test_sparse_matrix_assemble.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Dec 13 2010
 * @date last modification: Fri May 19 2017
 *
 * @brief  test the assembling method of the SparseMatrix class
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
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"

#include "dof_synchronizer.hh"
#include "sparse_matrix.hh"

/* -------------------------------------------------------------------------- */
using namespace akantu;

int main(int argc, char * argv[]) {
  initialize(argc, argv);

  UInt spatial_dimension = 2;
  Mesh mesh(spatial_dimension);
  mesh.read("triangle.msh");

  UInt nb_nodes = mesh.getNbNodes();

  DOFManagerDefault dof_manager(mesh, "test_dof_manager");

  Array<Real> test_synchronize(nb_nodes, spatial_dimension, "Test vector");
  dof_manager.registerDOFs("test_synchronize", test_synchronize, _dst_nodal);

  auto & A = dof_manager.getNewMatrix("A", _symmetric);

  // const akantu::Mesh::ConnectivityTypeList & type_list =
  // mesh.getConnectivityTypeList();
  // akantu::Mesh::ConnectivityTypeList::const_iterator it;

  // for(it = type_list.begin(); it != type_list.end(); ++it) {
  //   if(mesh.getSpatialDimension(*it) != spatial_dimension) continue;
  //   akantu::UInt nb_element           = mesh.getNbElement(*it);
  //   akantu::UInt nb_nodes_per_element = mesh.getNbNodesPerElement(*it);
  //   akantu::Element element(*it);

  //   akantu::UInt m = nb_nodes_per_element * spatial_dimension;
  //   akantu::Array<akantu::Real> local_mat(m, m, 1, "local_mat");

  //   for(akantu::UInt e = 0; e < nb_element; ++e) {
  //     element.element = e;
  //     sparse_matrix.addToMatrix(local_mat.storage(), element,
  //     nb_nodes_per_element);
  //   }
  // }

  A.saveMatrix("matrix.mtx");

  finalize();

  return EXIT_SUCCESS;
}
