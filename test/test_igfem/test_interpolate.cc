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
#include "dumper_paraview.hh"
#include "dumper_nodal_field.hh" 
/* -------------------------------------------------------------------------- */
#include <cstdlib>
#include <fstream>
#include <iostream>
/* -------------------------------------------------------------------------- */
using namespace akantu;

int main(int argc, char *argv[]) {
  initialize(argc, argv);

  const ElementType type_igfem = _igfem_triangle_4;
  const ElementType type_regular = _triangle_3;
  ///  debug::setDebugLevel(dblTest);
  UInt dim = ElementClass<type_igfem>::getSpatialDimension();

  // create a mesh that containes one IGFEM triangle 4
  UInt nb_elements_igfem = 1;
  UInt nb_elements_regular = 1;
  UInt nb_nodes = 5;
  Mesh mesh(dim);
  std::cout << "Generating mesh..." << std::endl;
  Array<Real> & nodes = const_cast<Array<Real> &>(mesh.getNodes());
  nodes.resize(nb_nodes);
  mesh.addConnectivityType(type_regular);
  mesh.addConnectivityType(type_igfem);
  Array<UInt> & connectivity_igfem = const_cast<Array<UInt> &>(mesh.getConnectivity(type_igfem));
  connectivity_igfem.resize(nb_elements_igfem);
  Array<UInt> & connectivity_regular = const_cast<Array<UInt> &>(mesh.getConnectivity(type_regular));
  connectivity_regular.resize(nb_elements_regular);
  // set the nodal coordinates
  nodes.storage()[0] = 0.;
  nodes.storage()[1] = 0.;
  nodes.storage()[2] = 1.;
  nodes.storage()[3] = 0.;
  nodes.storage()[4] = 0.;
  nodes.storage()[5] = 1.;
  nodes.storage()[6] = 0.5;
  nodes.storage()[7] = 0.5;
  nodes.storage()[8] = -1.;
  // set the element connectivity 
  // first element
  connectivity_igfem[0] = 0;
  connectivity_igfem[1] = 1;
  connectivity_igfem[2] = 2;
  connectivity_igfem[3] = 3;
  // second element
  connectivity_regular[0] = 0;
  connectivity_regular[1] = 2;
  connectivity_regular[2] = 4;

  FEEngine *fem = new FEEngineTemplate<IntegratorGauss,ShapeLagrange,_ek_igfem>(mesh, dim, "my_fem");

  std::stringstream outfilename; outfilename << "out_" << type_igfem << ".txt";
  std::ofstream my_file(outfilename.str().c_str());

  fem->initShapeFunctions();

  std::cout << *fem << std::endl;

  Array<Real> const_val(fem->getMesh().getNbNodes(), 2, "const_val");
 
  UInt nb_quadrature_points = fem->getNbQuadraturePoints(type_igfem) * nb_elements_igfem;

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
  const_val.storage()[8] = 1.;
  const_val.storage()[9] = 2.;
  fem->interpolateOnQuadraturePoints(const_val, val_on_quad, 2, type_igfem);

  my_file << const_val << std::endl;
  my_file << val_on_quad << std::endl;

  //  DumperParaview dumper_regular("mesh");
  DumperParaview dumper_igfem("mesh_igfem");
  dumper_igfem.registerMesh(mesh, dim, _not_ghost, _ek_igfem);
  DumperParaview dumper_regular("mesh_regular");
  dumper_regular.registerMesh(mesh, dim, _not_ghost, _ek_regular);

  //  dumper_regular.registerField("displacement", new dumper::NodalField<Real,false>(const_val, 0, 0));

  dumper_igfem.dump();
  dumper_regular.dump();
  //  dumper_regular.dump();
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
