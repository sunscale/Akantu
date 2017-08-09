/**
 * @file   test_fe_engine_precomputation.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Jun 14 2010
 * @date last modification: Mon Jul 13 2015
 *
 * @brief test integration on elements, this test consider that mesh is a cube
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "fe_engine.hh"
#include "shape_lagrange.hh"
#include "integrator_gauss.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */
using namespace akantu;

typedef FEEngineTemplate<IntegratorGauss, ShapeLagrange> FEM;
const ElementType type = TYPE;

//                        cst           x           x^2           x^3
//                        x^4          x^5
Real alpha[3][6] = {{0.40062394, 0.13703225, 0.51731446, 0.87830084, 0.5410543,
                     0.71842292}, // x
                    {0.41861835, 0.11080576, 0.49874043, 0.49077504, 0.85073835,
                     0.66259755}, // y
                    {0.92620845, 0.7503478, 0.62962232, 0.31662719, 0.64069644,
                     0.30878135}}; // z

static Vector<Real> eval_poly(UInt degree, const Vector<Real> & x) {
  Vector<Real> res(x.size());
  for (UInt i = 0; i < degree + 1; ++i) {
    for (UInt d = 0; d < x.size(); ++d) {
      res(d) += std::pow(x(d), i) * alpha[d][i];
    }
  }
  return res;
}

static Vector<Real> eval_int(UInt degree, const Vector<Real> & a,
                             const Vector<Real> & b) {
  Vector<Real> res(a.size());
  for (UInt i = 0; i < degree + 1; ++i) {
    for (UInt d = 0; d < a.size(); ++d) {
      res(d) += (std::pow(b(d), i + 1) - std::pow(a(d), i + 1)) * alpha[d][i] / Real(i + 1);
    }
  }

  if (a.size() == 3) {
    res(_x) *= std::abs(b(_y) - a(_y)) * std::abs(b(_z) - a(_z));
    res(_y) *= std::abs(b(_x) - a(_x)) * std::abs(b(_z) - a(_z));
    res(_z) *= std::abs(b(_y) - a(_y)) * std::abs(b(_x) - a(_x));
  } else if (a.size() == 2) {
    res(_x) *= std::abs(b(_y) - a(_y));
    res(_y) *= std::abs(b(_x) - a(_x));
  }

  return res;
}

template <UInt degree>
static Vector<Real> integrate_poly(UInt poly_degree, FEM & fem) {
  Mesh & mesh = fem.getMesh();
  UInt dim = mesh.getSpatialDimension();

  Matrix<Real> integration_points =
      fem.getIntegrator().getIntegrationPoints<type, degree>();

  // Vector<Real> integration_weights =
  //     fem.getIntegrator().getIntegrationWeights<type, degree>();

  // for (UInt i = 0; i < integration_points.cols(); ++i) {
  //   std::cout << "q(" << i << ") = " << Vector<Real>(integration_points(i))
  //   << " - w(" << i << ") = "<< integration_weights[i] << std::endl;
  // }

  UInt nb_integration_points = integration_points.cols();
  UInt nb_element = mesh.getNbElement(type);

  UInt shapes_size = ElementClass<type>::getShapeSize();
  Array<Real> shapes(0, shapes_size);
  fem.getShapeFunctions().computeShapesOnIntegrationPoints<type>(
      mesh.getNodes(), integration_points, shapes, _not_ghost);

  UInt vect_size = nb_integration_points * nb_element;
  Array<Real> integration_points_pos(vect_size, dim);
  fem.getShapeFunctions().interpolateOnIntegrationPoints<type>(
      mesh.getNodes(), integration_points_pos, dim, shapes);

  Array<Real> polynomial(vect_size, dim);

  Array<Real>::vector_iterator P_it = polynomial.begin(dim);
  Array<Real>::vector_iterator P_end = polynomial.end(dim);
  Array<Real>::const_vector_iterator x_it = integration_points_pos.begin(dim);

  for (; P_it != P_end; ++P_it, ++x_it) {
    *P_it = eval_poly(poly_degree, *x_it);
     //   std::cout << "Q = " << *x_it << std::endl;
     //   std::cout << "P(Q) = " << *P_it << std::endl;
  }

  Vector<Real> res(dim);

  Array<Real> polynomial_1d(vect_size, 1);
  for (UInt d = 0; d < dim; ++d) {
    Array<Real>::const_vector_iterator P_it = polynomial.begin(dim);
    Array<Real>::const_vector_iterator P_end = polynomial.end(dim);
    Array<Real>::scalar_iterator P1_it = polynomial_1d.begin();
    for (; P_it != P_end; ++P_it, ++P1_it) {
      *P1_it = (*P_it)(d);
      // std::cout << "P(Q, d) = " << *P1_it << std::endl;
    }
    res(d) = fem.getIntegrator().integrate<type, degree>(polynomial_1d);
  }

  return res;
}

int main(int argc, char * argv[]) {
  akantu::initialize(argc, argv);

  const UInt dim = ElementClass<type>::getSpatialDimension();
  Mesh mesh(dim);

  std::stringstream meshfilename;
  meshfilename << type << ".msh";
  mesh.read(meshfilename.str());

  FEM fem(mesh, dim, "my_fem");

  const Vector<Real> & lower = mesh.getLowerBounds();
  const Vector<Real> & upper = mesh.getUpperBounds();

  bool ok = true;

  for (UInt d = 0; d < 6; ++d) {
    Vector<Real> res(dim);
    switch (d) {
    case 0:
      res = integrate_poly<1>(d, fem);
      break;
    case 1:
      res = integrate_poly<1>(d, fem);
      break;
    case 2:
      res = integrate_poly<2>(d, fem);
      break;
    case 3:
      res = integrate_poly<3>(d, fem);
      break;
    case 4:
      res = integrate_poly<4>(d, fem);
      break;
    case 5:
      res = integrate_poly<5>(d, fem);
      break;
    }

    Vector<Real> exact(dim);
    exact = eval_int(d, lower, upper);

    Vector<Real> error(dim);
    error = exact - res;

    Real error_n = error.norm<L_inf>();
    if (error_n > 5e-14) {
      std::cout << d << " -> Resultat " << res << " - ";
      std::cout << "Exact" << exact << " -- ";
      std::cout << error << " {" << error_n << "}" << std::endl;
      ok = false;
    }
  }

  finalize();

  if (ok)
    return 0;
  else
    return 1;
}
