/**
 * @file   test_interpolate.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  test the interpolation function of the igfem elements
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "fe_engine.hh"
#include "integrator_gauss_igfem.hh"
#include "shape_igfem.hh"
/* -------------------------------------------------------------------------- */
#include <cstdlib>
#include <fstream>
#include <iostream>
/* -------------------------------------------------------------------------- */
using namespace akantu;

void interpolate(const ElementType type);
void generateIGFEMMesh(const ElementType type, Mesh & mesh,
                       const std::string & filename);

int main(int argc, char * argv[]) {
  akantu::initialize(argc, argv);
  debug::setDebugLevel(dblTest);

  /// interpolate test of _igfem_triangle_4
  AKANTU_DEBUG_INFO("integrating _igfem_triangle_4...");
  interpolate(_igfem_triangle_4);

  /// interpolate test of _igfem_triangle_5
  AKANTU_DEBUG_INFO("integrating _igfem_triangle_5...");
  interpolate(_igfem_triangle_5);

  finalize();
  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
void interpolate(const ElementType type) {

  UInt dim = 2;
  std::stringstream mesh_info;
  mesh_info << "mesh_info" << type << ".txt";
  Mesh my_mesh(dim);
  generateIGFEMMesh(type, my_mesh, mesh_info.str());

  Real eps = 3e-13;
  std::cout << "Epsilon : " << eps << std::endl;

  FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_igfem> * fem =
      new FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_igfem>(
          my_mesh, dim, "my_fem");

  fem->initShapeFunctions();

  Array<Real> const_val(fem->getMesh().getNbNodes(), 2, "const_val");

  UInt nb_element = my_mesh.getNbElement(type);
  UInt nb_quadrature_points = fem->getNbIntegrationPoints(type) * nb_element;

  Array<Real> val_on_quad(nb_quadrature_points, 2, "val_on_quad");

  /// the number of standard nodes in the mesh, i.e. not interface nodes
  UInt nb_standard_nodes = 9;

  /// impose constant value at standard nodes
  for (UInt i = 0; i < nb_standard_nodes; ++i) {
    const_val.storage()[i * 2 + 0] = 1.;
    const_val.storage()[i * 2 + 1] = 2.;
  }

  /// for field to be constant the enriched values need to be zero,
  /// because enrichment is not needed since there is no kink in the
  /// applied field
  for (UInt i = nb_standard_nodes; i < const_val.getSize(); ++i) {
    const_val.storage()[i * 2 + 0] = 0.;
    const_val.storage()[i * 2 + 1] = 0.;
  }

  fem->interpolateOnIntegrationPoints(const_val, val_on_quad, 2, type);

  std::cout << "Interpolation of array : " << const_val << std::endl;
  std::cout << "Gives on quads : " << val_on_quad << std::endl;

  /// interpolate coordinates
  Array<Real> coord_on_quad(nb_quadrature_points, my_mesh.getSpatialDimension(),
                            "coord_on_quad");
  /// create an array with the nodal coordinates that need to be
  /// interpolated. The nodal coordinates of the enriched nodes need
  /// to be set to zero, because they represent the enrichment of the
  /// position field, and the enrichments for this field are all zero!
  /// There is no gap in the mesh!
  Array<Real> igfem_nodes(my_mesh.getNbNodes(), dim);
  const ShapeLagrange<_ek_igfem> & shapes = fem->getShapeFunctions();
  shapes.extractValuesAtStandardNodes(my_mesh.getNodes(), igfem_nodes,
                                      _not_ghost);
  fem->interpolateOnIntegrationPoints(igfem_nodes, coord_on_quad,
                                      my_mesh.getSpatialDimension(), type);

  std::cout << "Interpolations of node coordinates : " << my_mesh.getNodes()
            << std::endl;
  std::cout << "Gives : " << coord_on_quad << std::endl;

  delete fem;
}
