/**
 * @file   test_solid_mechanics_model_material_eigenstrain.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sat Apr 16 2011
 * @date last modification: Thu Feb 01 2018
 *
 * @brief  test the internal field prestrain
 *
 * @section LICENSE
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
#include "mesh_utils.hh"
#include "non_linear_solver.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {
  initialize("material_elastic.dat", argc, argv);

  UInt dim = 3;

  /// load mesh
  Mesh mesh(dim);
  mesh.read("cube_3d_tet_4.msh");

  /// declaration of model
  SolidMechanicsModel model(mesh);

  /// model initialization
  // model.initFull(_analysis_method=akantu._explicit_lumped_mass)
  model.initFull(_analysis_method = _implicit_dynamic);
  // model.initFull(_analysis_method = akantu._implicit_dynamic)

  auto & solver = model.getNonLinearSolver();
  solver.set("threshold", 1e-4);
  solver.set("max_iterations", 100);
  solver.set("convergence_type", SolveConvergenceCriteria::_residual);

  const Array<Real> & coordinates = mesh.getNodes();
  Array<Real> & displacement = model.getDisplacement();
  Array<bool> & boundary = model.getBlockedDOFs();
  MeshUtils::buildFacets(mesh);

  /* ------------------------------------------------------------------------ */
  /* Dynamic eolution                                                         */
  /* ------------------------------------------------------------------------ */
  model.solveStep();
  std::cout << "Converged in " << Int(solver.get("nb_iterations")) << " ("
            << Real(solver.get("error")) << ")" << std::endl;

  finalize();

  return EXIT_SUCCESS;
}
