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
#include "dof_manager.hh"
#include "mesh.hh"
#include "mesh_accessor.hh"
#include "model_solver.hh"
#include "sparse_matrix.hh"
/* -------------------------------------------------------------------------- */
#include <fstream>
/* -------------------------------------------------------------------------- */

using namespace akantu;

class MyModel;
static void genMesh(Mesh & mesh, UInt nb_nodes);
static void printResults(MyModel & model);

/**
 *   =\o-----o-----o-> F
 *     |           |
 *     |---- L ----|
 */
class MyModel : public ModelSolver {
public:
  MyModel(Real F, Mesh & mesh)
      : ModelSolver(mesh, "model_solver", 0),
        dispacement(mesh.getNbNodes(), 1, "disp"),
        blocked(mesh.getNbNodes(), 1), forces(mesh.getNbNodes(), 1), mesh(mesh),
        nb_dofs(mesh.getNbNodes()), EA(1.) {

    this->initDOFManager("mumps");

    this->getDOFManager().registerDOFs("disp", dispacement, _dst_nodal);
    this->getDOFManager().registerBlockedDOFs("disp", blocked);
    this->getDOFManager().getNewMatrix("K", _symmetric);
    this->getDOFManager().getNewMatrix("J", "K");

    dispacement.set(0.);
    forces.set(0.);
    blocked.set(false);

    forces(nb_dofs - 1, _x) = F;
    blocked(0, _x) = true;
  }

  void assembleJacobian() {
    SparseMatrix & K = this->getDOFManager().getMatrix("K");
    K.clear();

    Matrix<Real> k(2, 2);
    k(0, 0) = k(1, 1) = 1;
    k(0, 1) = k(1, 0) = -1;

    Array<Real> k_all_el(this->nb_dofs - 1, 4);
    Array<Real>::matrix_iterator k_it = k_all_el.begin(2, 2);

    Array<UInt>::const_vector_iterator cit =
      this->mesh.getConnectivity(_segment_2).begin(2);
    Array<UInt>::const_vector_iterator cend =
        this->mesh.getConnectivity(_segment_2).end(2);

    for (; cit != cend; ++cit, ++k_it) {
      const Vector<UInt> & conn = *cit;
      UInt n1 = conn(0);
      UInt n2 = conn(1);

      Real p1 = this->mesh.getNodes()(n1, _x);
      Real p2 = this->mesh.getNodes()(n2, _x);

      Real L = std::abs(p2 - p1);

      Matrix<Real> & k_el = *k_it;
      k_el = k;
      k_el *= EA / L;
    }
    this->getDOFManager().assembleElementalMatricesToMatrix(
        "K", "disp", k_all_el, _segment_2);
  }

  void assembleResidual() {
    this->getDOFManager().assembleToResidual("disp", forces);

    Array<Real> forces_internal_el(this->nb_dofs - 1, 2);
    Array<Real>::vector_iterator f_it = forces_internal_el.begin(2);

    Array<UInt>::const_vector_iterator cit =
        this->mesh.getConnectivity(_segment_2).begin(2);
    Array<UInt>::const_vector_iterator cend =
        this->mesh.getConnectivity(_segment_2).end(2);

    for (; cit != cend; ++cit, ++f_it) {
      const Vector<UInt> & conn = *cit;
      UInt n1 = conn(0);
      UInt n2 = conn(1);
      Real p1 = this->mesh.getNodes()(n1, _x);
      Real p2 = this->mesh.getNodes()(n2, _x);

      Real L = std::abs(p2 - p1);

      Real u1 = this->dispacement(n1, _x);
      Real u2 = this->dispacement(n2, _x);

      Real f_n = EA / L * (u1 - u2);

      Vector<Real> & f = *f_it;

      f(0) = -f_n;
      f(1) = f_n;
    }

    this->getDOFManager().assembleElementalArrayToResidual(
        "disp", forces_internal_el, _segment_2, _not_ghost, -1.);
  }

  void predictor() {}
  void corrector() {}

  Array<Real> dispacement;
  Array<bool> blocked;
  Array<Real> forces;

private:
  Mesh & mesh;
  UInt nb_dofs;
  Real EA;
};

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  initialize(argc, argv);

  std::cout << std::setprecision(7);

  UInt nb_nodes = 11;
  Mesh mesh(1);

  genMesh(mesh, nb_nodes);

  MyModel model(10., mesh);

  model.getNewSolver("static", _tsst_static, _nls_linear);
  model.setIntegrationScheme("static", "disp", _ist_pseudo_time);

  model.solveStep();

  printResults(model);

  finalize();
  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
void genMesh(Mesh & mesh, UInt nb_nodes) {
  MeshAccessor mesh_accessor(mesh);
  Array<Real> & nodes = mesh_accessor.getNodes();
  Array<UInt> & conn = mesh_accessor.getConnectivity(_segment_2);

  nodes.resize(nb_nodes);

  for (UInt n = 0; n < nb_nodes; ++n) {
    nodes(n, _x) = n * (1. / (nb_nodes - 1));
  }

  conn.resize(nb_nodes - 1);
  for (UInt n = 0; n < nb_nodes - 1; ++n) {
    conn(n, 0) = n;
    conn(n, 1) = n + 1;
  }
}

/* -------------------------------------------------------------------------- */
void printResults(MyModel & model) {
  Array<Real>::const_scalar_iterator disp_it = model.dispacement.begin();
  Array<Real>::const_scalar_iterator force_it = model.forces.begin();
  Array<bool>::const_scalar_iterator blocked_it = model.blocked.begin();

  std::ofstream output("disp.csv");
  output << "node"
         << ", " << std::setw(8) << "disp"
         << ", " << std::setw(8) << "force"
         << ", " << std::setw(8) << "blocked" << std::endl;

  UInt node = 0;
  for (; disp_it != model.dispacement.end();
       ++disp_it, ++force_it, ++blocked_it, ++node) {
    output << node << ", "
           << std::setw(8) << *disp_it << ", "
           << std::setw(8) << *force_it << ", "
           << std::setw(8) << *blocked_it
           << std::endl;
  }
}
