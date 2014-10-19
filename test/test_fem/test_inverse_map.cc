/**
 * @file   test_inverse_map.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Fri May 25 2012
 * @date last modification: Tue Sep 02 2014
 *
 * @brief  test of the fem class
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
#include "fe_engine.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "shape_lagrange.hh"
#include "integrator_gauss.hh"
/* -------------------------------------------------------------------------- */
#include <cstdlib>
#include <fstream>
#include <iostream>
/* -------------------------------------------------------------------------- */
using namespace akantu;

int main(int argc, char *argv[]) {
  akantu::initialize(argc, argv);

  debug::setDebugLevel(dblTest);
  const ElementType type = TYPE;
  UInt dim = ElementClass<type>::getSpatialDimension();

  MeshIOMSH mesh_io;
  Mesh my_mesh(dim);

  my_mesh.computeBoundingBox();
  const Vector<Real> & lower = my_mesh.getLowerBounds();
  const Vector<Real> & upper = my_mesh.getUpperBounds();

  std::stringstream meshfilename; meshfilename << type << ".msh";
  mesh_io.read(meshfilename.str(), my_mesh);

  UInt nb_elements = my_mesh.getNbElement(type);
  ///
  FEEngineTemplate<IntegratorGauss,ShapeLagrange> *fem =
    new FEEngineTemplate<IntegratorGauss,ShapeLagrange>(my_mesh, dim, "my_fem");

  fem->initShapeFunctions();

  UInt nb_quad_points = fem->getNbQuadraturePoints(type);

  /// get the quadrature points coordinates
  Array<Real> coord_on_quad(nb_quad_points*nb_elements,
			     my_mesh.getSpatialDimension(),
			     "coord_on_quad");

  fem->interpolateOnQuadraturePoints(my_mesh.getNodes(),
				     coord_on_quad,
				     my_mesh.getSpatialDimension(),
				     type);


  /// loop over the quadrature points
  Array<Real>::iterator< Vector<Real> > it = coord_on_quad.begin(dim);
  Vector<Real> natural_coords(dim);

  Matrix<Real> quad = GaussIntegrationElement<type>::getQuadraturePoints();

   for(UInt el = 0 ; el < nb_elements ; ++el){
    for(UInt q = 0 ; q < nb_quad_points ; ++q){
      fem->inverseMap(*it, el, type, natural_coords);
      for (UInt i = 0; i < dim; ++i) {
	__attribute__ ((unused)) const Real eps = 1e-13;
	AKANTU_DEBUG_ASSERT(std::abs((natural_coords(i) - quad(i,q))/(upper(i)-lower(i))) < eps,
			    "real coordinates inversion test failed:"
			    << natural_coords(i) << " - " << quad(i, q)
			    << " = " << (natural_coords(i) - quad(i, q))/(upper(i)-lower(i)));
      }
      ++it;
    }
   }

   std::cout << "inverse completed over " << nb_elements << " elements" << std::endl;


  delete fem;
  finalize();

  return EXIT_SUCCESS;
}
