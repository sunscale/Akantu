/**
 * @file   test_gradient.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  test the gradient computation for the igfem elements
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

bool gradient(const ElementType type);
void generateIGFEMMesh(const ElementType type, Mesh & mesh,
                       const std::string & filename);

int main(int argc, char * argv[]) {
  akantu::initialize(argc, argv);
  debug::setDebugLevel(dblTest);

  bool test_passed = false;

  /// integrate test of _igfem_triangle_4
  AKANTU_DEBUG_INFO("integrating _igfem_triangle_4...");
  test_passed = gradient(_igfem_triangle_4);
  if (!test_passed) {
    finalize();
    return EXIT_FAILURE;
  }

  /// integrate test of _igfem_triangle_5
  AKANTU_DEBUG_INFO("integrating _igfem_triangle_5...");
  test_passed = gradient(_igfem_triangle_5);
  if (!test_passed) {
    finalize();
    return EXIT_FAILURE;
  }

  finalize();
  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
bool gradient(const ElementType type) {

  bool correct_result = true;
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
  Real alpha[2][3] = {{13, 23, 31}, {11, 7, 5}};

  /// create the 2 component field
  /// for field to be constant the enriched values need to be zero!! Therefore
  /// non-zero values
  /// are only imposed on standard nodes
  UInt nb_standard_nodes = 9;
  const Array<Real> & position = fem->getMesh().getNodes();
  Array<Real> const_val(fem->getMesh().getNbNodes(), 2, "const_val");

  UInt nb_element = my_mesh.getNbElement(type);
  UInt nb_quadrature_points = fem->getNbIntegrationPoints(type) * nb_element;

  Array<Real> grad_on_quad(nb_quadrature_points, 2 * dim, "grad_on_quad");
  for (UInt i = 0; i < nb_standard_nodes; ++i) {
    const_val(i, 0) = 0;
    const_val(i, 1) = 0;
    for (UInt d = 0; d < dim; ++d) {
      const_val(i, 0) += alpha[0][d] * position(i, d);
      const_val(i, 1) += alpha[1][d] * position(i, d);
    }
  }
  /// impose zero at enriched nodes
  for (UInt i = nb_standard_nodes; i < const_val.getSize(); ++i) {
    const_val(i, 0) = 0;
    const_val(i, 1) = 0;
  }

  /// compute the gradient
  fem->gradientOnIntegrationPoints(const_val, grad_on_quad, 2, type);

  std::cout << "Linear array on nodes : " << const_val << std::endl;
  std::cout << "Gradient on quad : " << grad_on_quad << std::endl;

  /// check the results
  Array<Real>::matrix_iterator it = grad_on_quad.begin(2, dim);
  Array<Real>::matrix_iterator it_end = grad_on_quad.end(2, dim);
  for (; it != it_end; ++it) {
    for (UInt d = 0; d < dim; ++d) {
      Matrix<Real> & grad = *it;
      if (!(std::abs(grad(0, d) - alpha[0][d]) < eps) ||
          !(std::abs(grad(1, d) - alpha[1][d]) < eps)) {
        std::cout << "Error gradient is not correct " << (*it)(0, d) << " "
                  << alpha[0][d] << " (" << std::abs((*it)(0, d) - alpha[0][d])
                  << ")"
                  << " - " << (*it)(1, d) << " " << alpha[1][d] << " ("
                  << std::abs((*it)(1, d) - alpha[1][d]) << ")"
                  << " - " << d << std::endl;
        std::cout << *it << std::endl;
        correct_result = false;
        return correct_result;
      }
    }
  }

  // compute gradient of coordinates
  Array<Real> grad_coord_on_quad(nb_quadrature_points, dim * dim,
                                 "grad_coord_on_quad");
  /// create an array with the nodal coordinates that need to be
  /// interpolated. The nodal coordinates of the enriched nodes need
  /// to be set to zero, because they represent the enrichment of the
  /// position field, and the enrichments for this field are all zero!
  /// There is no gap in the mesh!
  Array<Real> igfem_nodes(my_mesh.getNbNodes(), dim);
  const ShapeLagrange<_ek_igfem> & shapes = fem->getShapeFunctions();
  shapes.extractValuesAtStandardNodes(my_mesh.getNodes(), igfem_nodes,
                                      _not_ghost);

  fem->gradientOnIntegrationPoints(igfem_nodes, grad_coord_on_quad,
                                   my_mesh.getSpatialDimension(), type);

  std::cout << "Node positions : " << my_mesh.getNodes() << std::endl;
  std::cout << "Gradient of nodes : " << grad_coord_on_quad << std::endl;

  Array<Real>::matrix_iterator itp = grad_coord_on_quad.begin(dim, dim);
  Array<Real>::matrix_iterator itp_end = grad_coord_on_quad.end(dim, dim);

  for (; itp != itp_end; ++itp) {
    for (UInt i = 0; i < dim; ++i) {
      for (UInt j = 0; j < dim; ++j) {
        if (!(std::abs((*itp)(i, j) - (i == j)) < eps)) {
          std::cout << *itp << std::endl;
          correct_result = false;
          return correct_result;
        }
      }
    }
  }

  delete fem;
  return correct_result;
}
