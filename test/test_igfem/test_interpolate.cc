/**
 * @file   test_interpolate.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  test the interpolation function of the igfem elements
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "shape_igfem.hh"
#include "integrator_gauss_igfem.hh"
#include "fe_engine.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"


/* -------------------------------------------------------------------------- */
#include <cstdlib>
#include <fstream>
#include <iostream>
/* -------------------------------------------------------------------------- */
using namespace akantu;

int main(int argc, char *argv[]) {
  initialize(argc, argv);

  const ElementType type = _igfem_triangle_4;
  ///  debug::setDebugLevel(dblTest);
  UInt dim = ElementClass<type>::getSpatialDimension();

  // create a mesh that containes one IGFEM triangle 4
  UInt nb_elements = 1;
  UInt nb_nodes = 4;
  Mesh mesh(dim);
  std::cout << "Generating mesh..." << std::endl;
  Array<Real> & nodes = const_cast<Array<Real> &>(mesh.getNodes());
  nodes.resize(nb_nodes);
  mesh.addConnectivityType(type);
  Array<UInt> & connectivity = const_cast<Array<UInt> &>(mesh.getConnectivity(type));
  connectivity.resize(nb_elements);
  // set the nodal coordinates
  nodes.storage()[0] = 0.;
  nodes.storage()[1] = 0.;
  nodes.storage()[2] = 1.;
  nodes.storage()[3] = 0.;
  nodes.storage()[4] = 0.;
  nodes.storage()[5] = 1.;
  nodes.storage()[6] = 0.5;
  nodes.storage()[7] = 0.5;
  // set the element connectivity 
  connectivity[0] = 0;
  connectivity[1] = 1;
  connectivity[2] = 2;
  connectivity[3] = 3;

  FEEngine *fem = new FEEngineTemplate<IntegratorGauss,ShapeLagrange,_ek_igfem>(mesh, dim, "my_fem");

  std::stringstream outfilename; outfilename << "out_" << type << ".txt";
  std::ofstream my_file(outfilename.str().c_str());

  fem->initShapeFunctions();

  std::cout << *fem << std::endl;

  Array<Real> const_val(fem->getMesh().getNbNodes(), 2, "const_val");
 
  UInt nb_quadrature_points = fem->getNbQuadraturePoints(type) * nb_elements;

  Array<Real> val_on_quad(nb_quadrature_points, 2, "val_on_quad");

  // for (UInt i = 0; i < const_val.getSize(); ++i) {
  //   const_val.storage()[i * 2 + 0] = 1.;
  //   const_val.storage()[i * 2 + 1] = 2.;
  // }
  const_val.storage()[0] = 1.;
  const_val.storage()[1] = 2.;
  const_val.storage()[2] = 1.;
  const_val.storage()[3] = 2.;
  const_val.storage()[4] = 1.;
  const_val.storage()[5] = 2.;
  const_val.storage()[6] = 0.;
  const_val.storage()[7] = 0.;
  fem->interpolateOnQuadraturePoints(const_val, val_on_quad, 2, type);

  my_file << const_val << std::endl;
  my_file << val_on_quad << std::endl;


  // // interpolate coordinates
  // Array<Real> coord_on_quad(nb_quadrature_points, mesh.getSpatialDimension(), "coord_on_quad");

  // fem->interpolateOnQuadraturePoints(mesh.getNodes(),
  // 				     coord_on_quad,
  // 				     mesh.getSpatialDimension(),
  // 				     type);
  // my_file << mesh.getNodes() << std::endl;
  // my_file << coord_on_quad << std::endl;

  delete fem;

  finalize();
  return EXIT_SUCCESS;
}
