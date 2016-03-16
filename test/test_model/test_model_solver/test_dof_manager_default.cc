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
#include "sparse_matrix_aij.hh"

using namespace akantu;

/**
 *   =\o---o---o-> F
 */
class MySolverCallback : public SolverCallback {
public:
  MySolverCallback(Real F, DOFManagerDefault & dof_manager)
      : dof_manager(dof_manager),
        K(dynamic_cast<SparseMatrixAIJ &>(
            dof_manager.getNewMatrix("K", _symmetric))),
        dispacement(3, 1, "disp"), blocked(3, 1), forces(3, 1) {
    dof_manager.registerDOFs("disp", dispacement, _dst_generic);
    dof_manager.registerBlockedDOFs("disp", blocked);

    dispacement.set(0.);
    forces.set(0.);
    blocked.set(false);

    forces(2, _x) = F;
    blocked(0, _x) = true;

    K.addToProfile(0, 0);
    K.addToProfile(0, 1);
    K.addToProfile(1, 1);
    K.addToProfile(1, 2);
    K.addToProfile(2, 2);
  }

  void assembleJacobian() {
    K(0, 0) = 1.;
    K(1, 1) = 2.;
    K(2, 2) = 1.;
    K(0, 1) = -1.;
    K(1, 2) = -1.;
  }
  void assembleResidual() { dof_manager.assembleToResidual("disp", forces); }

  void predictor() {}
  void corrector() {}

  DOFManagerDefault & dof_manager;
  SparseMatrixAIJ & K;
  Array<Real> dispacement;
  Array<bool> blocked;
  Array<Real> forces;
};

int main(int argc, char * argv[]) {
  initialize(argc, argv);

  DOFManagerDefault dof_manager("test_dof_manager");
  MySolverCallback callback(10., dof_manager);

  dof_manager.getNewMatrix("J", "K");

  NonLinearSolver & nls =
      dof_manager.getNewNonLinearSolver("my_nls", _nls_linear);
  TimeStepSolver & tss =
      dof_manager.getNewTimeStepSolver("my_tss", _tsst_static, nls);
  tss.setIntegrationScheme("disp", _ist_pseudo_time);
  tss.solveStep(callback);

  Array<Real>::const_scalar_iterator disp_it = callback.dispacement.begin();
  Array<Real>::const_scalar_iterator force_it = callback.forces.begin();
  Array<bool>::const_scalar_iterator blocked_it = callback.blocked.begin();
  for (; disp_it != callback.dispacement.end(); ++disp_it, ++force_it, ++blocked_it) {
    std::cout << *disp_it << " " << *force_it << " " <<  *blocked_it << std::endl;
  }

  finalize();
  return EXIT_SUCCESS;
}
