/**
 * @file   test_material_elasto_plastic_linear_isotropic_hardening.cc
 *
 * @author Jaehyun Cho <jaehyun.cho@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Thu Dec 03 2015
 *
 * @brief  test for material type elasto plastic linear isotropic hardening using
 *         tension-compression test
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
#include "solid_mechanics_model.hh"
#include "non_linear_solver.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;

/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[]) {
  initialize("test_material_elasto_plastic_linear_isotropic_hardening.dat", argc, argv);

  const UInt spatial_dimension = 2;
  const Real u_increment = 0.1;
  const UInt steps = 20;

  Mesh mesh(spatial_dimension);
  mesh.read("test_material_elasto_plastic_linear_isotropic_hardening.msh");

  SolidMechanicsModel model(mesh);
  model.initFull(_analysis_method = _static);

  model.getNonLinearSolver("static").set("max_iterations", 300);
  model.getNonLinearSolver("static").set("threshold", 1e-6);

  model.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "left");
  model.applyBC(BC::Dirichlet::FixedValue(0.0, _y), "bottom");

  std::cout.precision(4);
  for (UInt i = 0 ; i < steps ; ++i) {
    model.applyBC(BC::Dirichlet::FixedValue(i * u_increment, _x), "right");
    model.solveStep();

    Real strainxx = i * u_increment / 10.;

    const Array<UInt> & edge_nodes = mesh.getElementGroup("right").getNodes();
    Array<Real> & residual = model.getInternalForce();
    Real reaction = 0;

    for (UInt n = 0 ; n < edge_nodes.size() ; n++) {
      reaction -= residual(edge_nodes(n), 0);
    }

    std::cout << strainxx << "," << reaction << std::endl;
  }

  finalize();
  return EXIT_SUCCESS;
}
