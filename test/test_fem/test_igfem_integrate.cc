#include "aka_common.hh"
#include "shape_igfem.hh"
#include "integrator_gauss.hh"
#include "shape_igfem.hh"
/* -------------------------------------------------------------------------- */
#include <cstdlib>
#include <fstream>
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[]) {
  akantu::initialize(argc, argv);
  debug::setDebugLevel(dblTest);
  const ElementType type = _igfem_triangle_3;
  UInt dim = ElementClass<type>::getSpatialDimension();

  Real eps = 3e-13;
  std::cout << "Epsilon : " << eps << std::endl;

  Mesh my_mesh(dim);

  std::stringstream meshfilename; meshfilename << "_triangle_3.msh";
  my_mesh.read(meshfilename.str());

  FEEngine *fem = new FEEngineTemplate<IntegratorGauss,ShapeLagrange,_ek_igfem>(my_mesh, dim, "my_fem");

  std::stringstream outfilename; outfilename << "out_" << type << ".txt";
  std::ofstream my_file(outfilename.str().c_str());

  fem->initShapeFunctions();

  std::cout << *fem << std::endl;

  UInt nb_element = my_mesh.getNbElement(type);
  UInt nb_quadrature_points = fem->getNbQuadraturePoints(type) * nb_element;

  Array<Real> const_val(fem->getMesh().getNbNodes(), 2, "const_val");
  Array<Real> val_on_quad(nb_quadrature_points, 2 , "val_on_quad");

  for (UInt i = 0; i < const_val.getSize(); ++i) {
    const_val.storage()[i * 2 + 0] = 1.;
    const_val.storage()[i * 2 + 1] = 2.;
  }

  //interpolate function on quadrature points
  fem->interpolateOnQuadraturePoints(const_val, val_on_quad, 2, type);

  //integrate function on elements
  akantu::Array<akantu::Real> int_val_on_elem(nb_element, 2, "int_val_on_elem");
  fem->integrate(val_on_quad, int_val_on_elem, 2, type);

  // get global integration value
  Real value[2] = {0,0};
  my_file << val_on_quad << std::endl << int_val_on_elem << std::endl;
  for (UInt i = 0; i < fem->getMesh().getNbElement(type); ++i) {
    value[0] += int_val_on_elem.storage()[2*i];
    value[1] += int_val_on_elem.storage()[2*i+1];
  }

  my_file << "integral on the mesh of 1 is " << value[0] << " and of 2 is " << value[1] << std::endl;


  delete fem;
  finalize();

  if(!(std::abs(value[0] - 1.) < eps && std::abs(value[1] - 2.) < eps)) {
    std::cout << "|1 - " << value[0] << "| = " << std::abs(value[0] - 1.) << std::endl
	      << "|2 - " << value[1] << "| = " << std::abs(value[1] - 2.) << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
