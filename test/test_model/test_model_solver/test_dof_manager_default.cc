/**
 * @file   test_dof_manager_default.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Feb 24 12:28:44 2016
 *
 * @brief  Test default dof manager
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "dof_manager_default.hh"
#include "solver_callback.hh"
#include "time_step_solver.hh"
#include "sparse_matrix_aij.hh"

using namespace akantu;

/**
 *   =\o-----o-----o-> F
 *     |           |
 *     |---- L ----|
 */
class MySolverCallback : public SolverCallback {
public:
  MySolverCallback(Real F, DOFManagerDefault & dof_manager, UInt nb_dofs = 3)
      : dof_manager(dof_manager),
        K(dynamic_cast<SparseMatrixAIJ &>(
            dof_manager.getNewMatrix("K", _symmetric))),
        dispacement(nb_dofs, 1, "disp"), blocked(nb_dofs, 1),
        forces(nb_dofs, 1), nb_dofs(nb_dofs) {
    dof_manager.registerDOFs("disp", dispacement, _dst_generic);
    dof_manager.registerBlockedDOFs("disp", blocked);

    dispacement.set(0.);
    forces.set(0.);
    blocked.set(false);

    forces(nb_dofs - 1, _x) = F;
    blocked(0, _x) = true;

    for (UInt i = 0; i < nb_dofs; ++i)
      K.addToProfile(i, i);
    for (UInt i = 0; i < nb_dofs - 1; ++i)
      K.addToProfile(i, i + 1);
  }

  void assembleJacobian() {
    for (UInt i = 1; i < nb_dofs - 1; ++i)
      K(i, i) = 2;
    for (UInt i = 0; i < nb_dofs - 1; ++i)
      K(i, i + 1) = -1;

    K(0, 0) = K(nb_dofs - 1, nb_dofs - 1) = 1;

    // K *= 1 / L_{el}
    K *= nb_dofs - 1;
  }

  void assembleResidual() { dof_manager.assembleToResidual("disp", forces); }

  void predictor() {}
  void corrector() {}

  DOFManagerDefault & dof_manager;
  SparseMatrixAIJ & K;
  Array<Real> dispacement;
  Array<bool> blocked;
  Array<Real> forces;

  UInt nb_dofs;
};

int main(int argc, char * argv[]) {
  initialize(argc, argv);

  DOFManagerDefault dof_manager("test_dof_manager");
  MySolverCallback callback(10., dof_manager, 11);

  dof_manager.getNewMatrix("J", "K");

  NonLinearSolver & nls =
      dof_manager.getNewNonLinearSolver("my_nls", _nls_linear);
  TimeStepSolver & tss =
      dof_manager.getNewTimeStepSolver("my_tss", _tsst_static, nls);
  tss.setIntegrationScheme("disp", _ist_pseudo_time);
  tss.solveStep(callback);

  dof_manager.getMatrix("K").saveMatrix("K_dof_manager_default.mtx");

  Array<Real>::const_scalar_iterator disp_it = callback.dispacement.begin();
  Array<Real>::const_scalar_iterator force_it = callback.forces.begin();
  Array<bool>::const_scalar_iterator blocked_it = callback.blocked.begin();
  std::cout << std::setw(8) << "disp"
            << " " << std::setw(8) << "force"
            << " " << std::setw(8) << "blocked" << std::endl;

  for (; disp_it != callback.dispacement.end();
       ++disp_it, ++force_it, ++blocked_it) {
    std::cout << std::setw(8) << *disp_it << " " << std::setw(8) << *force_it
              << " " << std::setw(8) << std::boolalpha << *blocked_it << std::endl;
  }

  finalize();
  return EXIT_SUCCESS;
}
