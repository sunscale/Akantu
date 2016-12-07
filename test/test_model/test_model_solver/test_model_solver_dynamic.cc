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

#ifndef EXPLICIT
#define EXPLICIT true
#endif

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
  MyModel(Real F, Mesh & mesh, bool lumped)
      : ModelSolver(mesh, "model_solver", 0),
        displacement(mesh.getNbNodes(), 1, "disp"),
        velocity(mesh.getNbNodes(), 1, "velo"),
        acceleration(mesh.getNbNodes(), 1, "accel"),
        blocked(mesh.getNbNodes(), 1, "blocked"),
        forces(mesh.getNbNodes(), 1, "force_ext"),
        internal_forces(mesh.getNbNodes(), 1, "force_int"), mesh(mesh),
        nb_dofs(mesh.getNbNodes()), E(1.), A(1.), rho(1.), lumped(lumped) {

    this->initDOFManager();

    this->getDOFManager().registerDOFs("disp", displacement, _dst_nodal);
    this->getDOFManager().registerDOFsDerivative("disp", 1, velocity);
    this->getDOFManager().registerDOFsDerivative("disp", 2, acceleration);

    this->getDOFManager().registerBlockedDOFs("disp", blocked);

    this->getDOFManager().getNewMatrix("K", _symmetric);
    this->getDOFManager().getNewMatrix("M", "K");
    this->getDOFManager().getNewMatrix("J", "K");

    this->getDOFManager().getNewLumpedMatrix("M");

    if (lumped) {
      this->assembleLumpedMass();
    } else {
      this->assembleMass();
      this->assembleStiffness();
    }
    this->assembleJacobian();

    displacement.set(0.);
    velocity.set(0.);
    acceleration.set(0.);

    forces.set(0.);
    blocked.set(false);

    forces(nb_dofs - 1, _x) = F;
    // displacement(nb_dofs - 1, _x) = .1;
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

      Real M_n = rho * A * L / 2;
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

    Matrix<Real> m(2, 2);
    m(0, 0) = m(1, 1) = 2;
    m(0, 1) = m(1, 0) = 1;

    // under integrated
    // m(0, 0) = m(1, 1) = 3./2.;
    // m(0, 1) = m(1, 0) = 3./2.;

    // lumping the mass matrix
    // m(0, 0) += m(0, 1);
    // m(1, 1) += m(1, 0);
    // m(0, 1) = m(1, 0) = 0;

    for (; cit != cend; ++cit, ++m_it) {
      const Vector<UInt> & conn = *cit;
      UInt n1 = conn(0);
      UInt n2 = conn(1);

      Real p1 = this->mesh.getNodes()(n1, _x);
      Real p2 = this->mesh.getNodes()(n2, _x);

      Real L = std::abs(p2 - p1);

      Matrix<Real> & m_el = *m_it;
      m_el = m;
      m_el *= rho * A * L / 6.;
    }
    this->getDOFManager().assembleElementalMatricesToMatrix(
        "M", "disp", m_all_el, _segment_2);
  }

  void assembleJacobian() {}
  void assembleStiffness() {
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

      Real u1 = this->displacement(n1, _x);
      Real u2 = this->displacement(n2, _x);

      Real f_n = E * A / L * (u2 - u1);

      Vector<Real> & f = *f_it;

      f(0) = -f_n;
      f(1) = f_n;
    }

    internal_forces.clear();
    this->getDOFManager().assembleElementalArrayLocalArray(
        forces_internal_el, internal_forces, _segment_2, _not_ghost);
    this->getDOFManager().assembleToResidual("disp", internal_forces, -1.);
  }

  Real getPotentialEnergy() {
    Real res = 0;

    if (!lumped) {
      res = this->mulVectMatVect(this->displacement,
                                 this->getDOFManager().getMatrix("K"),
                                 this->displacement);
    } else {
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

        res += strain * E * strain * A * L;
      }
    }
    return res / 2.;
  }

  Real getKineticEnergy() {
    Real res = 0;
    if (!lumped) {
      res = this->mulVectMatVect(
          this->velocity, this->getDOFManager().getMatrix("M"), this->velocity);
    } else {
      Array<Real> & m = this->getDOFManager().getLumpedMatrix("M");
      Array<Real>::const_scalar_iterator it = velocity.begin();
      Array<Real>::const_scalar_iterator end = velocity.end();
      Array<Real>::const_scalar_iterator m_it = m.begin();

      for (; it != end; ++it, ++m_it) {
        res += *m_it * *it * *it;
      }
    }
    return res / 2.;
  }

  Real mulVectMatVect(const Array<Real> & x, const SparseMatrix & A,
                      const Array<Real> & y) {
    Array<Real> Ay(this->nb_dofs);
    A.matVecMul(y, Ay);
    Real res = 0.;

    Array<Real>::const_scalar_iterator it = Ay.begin();
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

  bool lumped;
};

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  initialize(argc, argv);

  UInt nb_nodes = 201;
  UInt max_steps = 2000;
  Real time_step = 0.001;
  Mesh mesh(1);
  Real F = -9.81;
  bool _explicit = EXPLICIT;

  genMesh(mesh, nb_nodes);

  MyModel model(F, mesh, _explicit);

  if (!_explicit) {
    model.getNewSolver("dynamic", _tsst_dynamic, _nls_newton_raphson);
    model.setIntegrationScheme("dynamic", "disp", _ist_trapezoidal_rule_2,
                               IntegrationScheme::_displacement);
  } else {
    model.getNewSolver("dynamic", _tsst_dynamic_lumped, _nls_lumped);
    model.setIntegrationScheme("dynamic", "disp", _ist_central_difference,
                               IntegrationScheme::_acceleration);
  }

  model.setTimeStep(time_step);

  const Array<Real> & disp = model.displacement;
  const Array<Real> & velo = model.velocity;
  const Array<Real> & reac = model.internal_forces;

// #if EXPLICIT == true
//   std::ofstream output("output_dynamic_explicit.csv");
// #else
//   std::ofstream output("output_dynamic_implicit.csv");
// #endif

  std::cout << std::setw(8) << "time"
            << "," << std::setw(8) << "disp"
            << "," << std::setw(8) << "velo"
            << "," << std::setw(8) << "reac"
            << "," << std::setw(8) << "wext"
            << "," << std::setw(8) << "epot"
            << "," << std::setw(8) << "ekin"
            << "," << std::setw(8) << "total"
            << std::endl;

  Real wext = 0;

  Real epot = model.getPotentialEnergy();
  Real ekin = model.getKineticEnergy();
  Real einit = ekin +  epot;
  std::cout << std::setw(8) << 0.
            << "," << std::setw(8) << disp(nb_nodes - 1, _x)
            << "," << std::setw(8) << velo(nb_nodes - 1, _x)
            << "," << std::setw(8) << (-reac(0, _x))
            << "," << std::setw(8) << wext
            << "," << std::setw(8) << epot
            << "," << std::setw(8) << ekin
            << "," << std::setw(8) << (ekin + epot - wext - einit)
            << std::endl;

#if EXPLICIT == false
  // NonLinearSolver & solver =
  //   model.getDOFManager().getNonLinearSolver("dynamic");
#endif

  for (UInt i = 1; i < max_steps + 1; ++i) {
    model.solveStep();

// #if EXPLICIT == false
    //     UInt nb_iter = solver.get("nb_iterations");
//     Real error = solver.get("error");
//     bool converged = solver.get("converged");
//     std::cout << error << " " << nb_iter << " -> " << converged << std::endl;
// #endif

    wext += F * velo(nb_nodes - 1, 0) * time_step;
    epot = model.getPotentialEnergy();
    ekin = model.getKineticEnergy();
    Real etot = ekin + epot - wext - einit;

    std::cout << std::setw(8) << time_step * i
              << "," << std::setw(8) << disp(nb_nodes - 1, _x)
              << "," << std::setw(8) << velo(nb_nodes - 1, _x)
              << "," << std::setw(8) << (-reac(0, _x))
              << "," << std::setw(8) << wext
              << "," << std::setw(8) << epot
              << "," << std::setw(8) << ekin
              << "," << std::setw(8) << (etot)
              << std::endl;

    if (std::abs(etot) > 1e-7) {
      AKANTU_DEBUG_ERROR("The total energy of the system is not conserved!");
    }
  }

  // output.close();
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
