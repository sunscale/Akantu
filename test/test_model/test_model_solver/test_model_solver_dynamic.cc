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
#include "model_solver.hh"
#include "mesh.hh"
#include "mesh_accessor.hh"
#include "dof_manager.hh"
#include "sparse_matrix.hh"

using namespace akantu;

class MyModel;
static void genMesh(Mesh & mesh, UInt nb_nodes);

/**
 *   =\o-----o-----o-> F
 *     |           |
 *     |---- L ----|
 */
class MyModel : public ModelSolver {
public:
  MyModel(Real F, Mesh & mesh)
      : ModelSolver(mesh, "model_solver", 0),
        displacement(mesh.getNbNodes(), 1, "disp"),
        velocity(mesh.getNbNodes(), 1, "velo"),
        acceleration(mesh.getNbNodes(), 1, "accel"),
        blocked(mesh.getNbNodes(), 1, "blocked"),
        forces(mesh.getNbNodes(), 1, "force_ext"),
        internal_forces(mesh.getNbNodes(), 1, "force_int"), mesh(mesh),
        nb_dofs(mesh.getNbNodes()), E(1.), A(1.), rho(1.) {

    this->initDOFManager("mumps");

    this->getDOFManager().registerDOFs("disp", displacement, _dst_nodal);
    this->getDOFManager().registerDOFsDerivative("disp", 1, velocity);
    this->getDOFManager().registerDOFsDerivative("disp", 2, acceleration);

    this->getDOFManager().registerBlockedDOFs("disp", blocked);

    this->getDOFManager().getNewMatrix("K", _symmetric);
    this->getDOFManager().getNewMatrix("M", "K");
    this->getDOFManager().getNewMatrix("J", "K");

    this->getDOFManager().getNewLumpedMatrix("M");

    this->assembleLumpedMass();
    this->assembleMass();
    this->assembleJacobian();

    displacement.set(0.);
    velocity.set(0.);
    acceleration.set(0.);

    forces.set(0.);
    blocked.set(false);

    //    forces(nb_dofs - 1, _x) = F;
    displacement(nb_dofs - 1, _x) = 1;
    blocked(0, _x) = true;
  }

  void assembleLumpedMass() {
    Array<Real> & M = this->getDOFManager().getLumpedMatrix("M");
    M.clear();

    Array<Real> m_all_el(this->nb_dofs - 1, 2);
    Array<Real>::vector_iterator m_it = m_all_el.begin(2);

    Array<UInt>::const_vector_iterator cit =
        this->mesh.getConnectivity(_segment_2).begin(2);
    Array<UInt>::const_vector_iterator cend =
        this->mesh.getConnectivity(_segment_2).end(2);

    for (; cit != cend; ++cit, ++m_it) {
      const Vector<UInt> & conn = *cit;
      UInt n1 = conn(0);
      UInt n2 = conn(1);

      Real p1 = this->mesh.getNodes()(n1, _x);
      Real p2 = this->mesh.getNodes()(n2, _x);

      Real L = std::abs(p2 - p1);

      Real M_n = rho * A * L / 2.;
      (*m_it)(0) = (*m_it)(1) = M_n;
    }

    this->getDOFManager().assembleElementalArrayLocalArray(
        m_all_el, M, _segment_2, _not_ghost);
  }


  void assembleMass() {
    SparseMatrix & M = this->getDOFManager().getMatrix("M");
    M.clear();

    Array<Real> m_all_el(this->nb_dofs - 1, 4);
    Array<Real>::matrix_iterator m_it = m_all_el.begin(2, 2);

    Array<UInt>::const_vector_iterator cit =
        this->mesh.getConnectivity(_segment_2).begin(2);
    Array<UInt>::const_vector_iterator cend =
        this->mesh.getConnectivity(_segment_2).end(2);

    for (; cit != cend; ++cit, ++m_it) {
      const Vector<UInt> & conn = *cit;
      UInt n1 = conn(0);
      UInt n2 = conn(1);

      Real p1 = this->mesh.getNodes()(n1, _x);
      Real p2 = this->mesh.getNodes()(n2, _x);

      Real L = std::abs(p2 - p1);

      Matrix<Real> & m_el = *m_it;
      m_el.set(rho * A * L / 4.);
    }
    this->getDOFManager().assembleElementalMatricesToMatrix(
        "M", "disp", m_all_el, _segment_2);
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
      k_el *= E * A / L;
    }
    this->getDOFManager().assembleElementalMatricesToMatrix(
        "K", "disp", k_all_el, _segment_2);
  }

  void assembleResidual() {
    this->getDOFManager().assembleToResidual("disp", forces);
    this->getDOFManager().assembleMatMulVectToResidual("disp", "K",
                                                       this->displacement,
                                                       -1);
  }



  Real getPotentialEnergy() {
    // return (1./2. * this->mulVectMatVect(this->displacement, this->getDOFManager().getMatrix("K"),
    //                                      this->displacement));

    Real res = 0;
    Array<UInt>::const_vector_iterator cit =
        this->mesh.getConnectivity(_segment_2).begin(2);
    Array<UInt>::const_vector_iterator cend =
        this->mesh.getConnectivity(_segment_2).end(2);

    for (; cit != cend; ++cit) {
      const Vector<UInt> & conn = *cit;
      UInt n1 = conn(0);
      UInt n2 = conn(1);
      Real p1 = this->mesh.getNodes()(n1, _x);
      Real p2 = this->mesh.getNodes()(n2, _x);

      Real L = std::abs(p2 - p1);

      Real u1 = this->displacement(n1, _x);
      Real u2 = this->displacement(n2, _x);

      Real strain = (u2 - u1) / L;

      res += strain * E * strain;
    }

    return res / 2.;
  }

  Real getKineticEnergy() {
    // return (1./2. * this->mulVectMatVect(this->velocity, this->getDOFManager().getMatrix("K"),
    //                                      this->velocity));
    Real res = 0;
    Array<Real> &m = this->getDOFManager().getLumpedMatrix("M");
    Array<Real>::const_scalar_iterator it  = velocity.begin();
    Array<Real>::const_scalar_iterator end = velocity.end();
    Array<Real>::const_scalar_iterator m_it = m.begin();

    for (; it != end; ++it, ++m_it) {
      res += *m_it * *it * *it;
    }
    return res / 2.;
  }

  Real mulVectMatVect(const Array<Real> & x, const SparseMatrix & A, const Array<Real> & y) {
    Array<Real> Ay(this->nb_dofs);
    A.matVecMul(y, Ay);
    Real res = 0.;

    Array<Real>::const_scalar_iterator it  = Ay.begin();
    Array<Real>::const_scalar_iterator end = Ay.end();
    Array<Real>::const_scalar_iterator x_it = x.begin();

    for (; it != end; ++it, ++x_it) {
      res += *x_it * *it;
    }
    return res;
  }

  void predictor() {}
  void corrector() {}

  Array<Real> displacement;
  Array<Real> velocity;
  Array<Real> acceleration;
  Array<bool> blocked;
  Array<Real> forces;
  Array<Real> internal_forces;

private:
  Mesh & mesh;
  UInt nb_dofs;
  Real E, A, rho;
};

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  initialize(argc, argv);

  UInt nb_nodes = 21;
  UInt max_steps = 500;
  Real time_step = 0.01;
  Mesh mesh(1);

  genMesh(mesh, nb_nodes);

  MyModel model(10., mesh);

  //  model.getNewSolver("dynamic", _tsst_dynamic, _nls_newton_raphson);
  //  model.setIntegrationScheme("dynamic", "disp", _ist_trapezoidal_rule_2,
  //                             IntegrationScheme::_displacement);
  model.getNewSolver("dynamic", _tsst_dynamic, _nls_lumped);
  model.setIntegrationScheme("dynamic", "disp", _ist_central_difference,
                             IntegrationScheme::_acceleration);

  model.setTimeStep(time_step);

  const Array<Real> & disp = model.displacement;
  std::cout << std::setw(8) << "time" << ", "
            << std::setw(8) << "disp" << ", "
            << std::setw(8) << "epot" << ", "
            << std::setw(8) << "ekin" << std::endl;

  for (UInt i = 1; i < max_steps + 1; ++i) {
    model.solveStep();


    Real epot = model.getPotentialEnergy();
    Real ekin = model.getKineticEnergy();
    std::cout << std::setw(8) << time_step * i << ", "
              << std::setw(8) << disp(nb_nodes - 1, _x) << ", "
              << std::setw(8) << epot << ", "
              << std::setw(8) << ekin << std::endl;
  }

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
