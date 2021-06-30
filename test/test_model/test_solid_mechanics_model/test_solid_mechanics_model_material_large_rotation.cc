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
#include "sparse_matrix_aij.hh"
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

  model.setBaseName("waves");
  model.addDumpFieldVector("displacement");
  model.addDumpFieldVector("acceleration");
  model.addDumpFieldVector("velocity");
  model.addDumpFieldVector("internal_force");
  model.addDumpFieldVector("external_force");
  model.addDumpField("strain");
  model.addDumpField("stress");
  model.addDumpField("blocked_dofs");

  /* ------------------------------------------------------------------------ */
  // get mass center
  /* ------------------------------------------------------------------------ */

  model.assembleMass();
  auto & M = model.getDOFManager().getMatrix("M");
  Array<Real> _mass(M.size(), 1);
  _mass.zero();
  std::cout << "AAAA " << M.size() << std::endl;
  std::cout << "AAAA " << _mass.size() << std::endl;

  for (UInt i = 0; i < M.size(); ++i) {
    for (UInt j = 0; j < M.size(); ++j) {
      std::cout << i << ", " << j <<std::endl;
      _mass[i] += M(i, j);
    }
  }
  std::array<Real, 3> mass_center{0., 0., 0.};
  std::cout << "AAAA " << _mass.size() << std::endl;
  Real total_mass = 0.;
  for (UInt i = 0; i < _mass.size(); ++i) {
    for (UInt j = 0; j < 3; ++j) {
      mass_center[j] += _mass(i * 3 + j);
      total_mass += _mass(i * 3 + j);
    }
  }
  mass_center[0] /= total_mass / 3.;
  mass_center[1] /= total_mass / 3.;
  mass_center[2] /= total_mass / 3.;
  std::cout << "total mass" << total_mass << std::endl;
  std::cout << mass_center[0] << " " << mass_center[1] << " " << mass_center[2]
            << std::endl;

  /* ---------------------------------------------------------------------- */
  /* Dynamic evolution */
  /* ---------------------------------------------------------------------- */
  model.dump();

  model.solveStep();
  model.dump();

  std::cout << "Converged in " << Int(solver.get("nb_iterations")) << " ("
            << Real(solver.get("error")) << ")" << std::endl;

  finalize();

  return EXIT_SUCCESS;
}
