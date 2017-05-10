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
#include "non_linear_solver.hh"
#include "sparse_matrix.hh"
#include "static_communicator.hh"
#include "synchronizer_registry.hh"
#include "data_accessor.hh"
#include "element_synchronizer.hh"
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
class MyModel : public ModelSolver, public DataAccessor<Element> {
public:
  MyModel(Real F, Mesh & mesh, bool lumped)
      : ModelSolver(mesh, "model_solver", 0),
        mesh(mesh), nb_dofs(mesh.getNbNodes()),
        nb_elements(mesh.getNbElement()),
        E(1.), A(1.), rho(1.), lumped(lumped),
        displacement(nb_dofs, 1, "disp"),
        velocity(nb_dofs, 1, "velo"),
        acceleration(nb_dofs, 1, "accel"),
        blocked(nb_dofs, 1, "blocked"),
        forces(nb_dofs, 1, "force_ext"),
        internal_forces(nb_dofs, 1, "force_int"),
        stresses(nb_elements, 1, "stress"),
        strains(nb_elements, 1, "strain"),
        initial_lengths(nb_elements, 1, "L0") {
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

    UInt global_nb_nodes = mesh.getNbGlobalNodes();
    for(UInt n = 0; n < nb_dofs; ++n) {
      if (mesh.getGlobalNodesIds()(n) == (global_nb_nodes -1))
        forces(n, _x) = F;

      if (mesh.getGlobalNodesIds()(n) == 0)
        blocked(n, _x) = true;
    }

    auto cit  = this->mesh.getConnectivity(_segment_2).begin(2);
    auto cend = this->mesh.getConnectivity(_segment_2).end(2);
    auto L_it = this->initial_lengths.begin();

    for (; cit != cend; ++cit, ++L_it) {
      const Vector<UInt> & conn = *cit;
      UInt n1 = conn(0);
      UInt n2 = conn(1);
      Real p1 = this->mesh.getNodes()(n1, _x);
      Real p2 = this->mesh.getNodes()(n2, _x);

      *L_it = std::abs(p2 - p1);
    }

    this->registerDataAccessor(*this);
    this->registerSynchronizer(const_cast<ElementSynchronizer &>(this->mesh.getElementSynchronizer()),
                               _gst_user_1);
  }

  void assembleLumpedMass() {
    Array<Real> & M = this->getDOFManager().getLumpedMatrix("M");
    M.clear();

    this->assembleLumpedMass(_not_ghost);
    if (this->mesh.getNbElement(_segment_2, _ghost) > 0)
      this->assembleLumpedMass(_ghost);
  }

  void assembleLumpedMass(const GhostType & ghost_type) {
    Array<Real> & M = this->getDOFManager().getLumpedMatrix("M");

    Array<Real> m_all_el(this->mesh.getNbElement(_segment_2, ghost_type), 2);

    Array<Real>::vector_iterator m_it = m_all_el.begin(2);

    auto cit  = this->mesh.getConnectivity(_segment_2, ghost_type).begin(2);
    auto cend = this->mesh.getConnectivity(_segment_2, ghost_type).end(2);

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
        m_all_el, M, _segment_2, ghost_type);
  }

  void assembleMass() {
    SparseMatrix & M = this->getDOFManager().getMatrix("M");
    M.clear();

    Array<Real> m_all_el(this->nb_elements, 4);
    Array<Real>::matrix_iterator m_it = m_all_el.begin(2, 2);

    auto cit  = this->mesh.getConnectivity(_segment_2).begin(2);
    auto cend = this->mesh.getConnectivity(_segment_2).end(2);

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
    k(0, 0) = k(1, 1) =  1;
    k(0, 1) = k(1, 0) = -1;

    Array<Real> k_all_el(this->nb_elements, 4);

    auto k_it = k_all_el.begin(2, 2);
    auto cit  = this->mesh.getConnectivity(_segment_2).begin(2);
    auto cend = this->mesh.getConnectivity(_segment_2).end(2);

    for (; cit != cend; ++cit, ++k_it) {
      const auto & conn = *cit;
      UInt n1 = conn(0);
      UInt n2 = conn(1);

      Real p1 = this->mesh.getNodes()(n1, _x);
      Real p2 = this->mesh.getNodes()(n2, _x);

      Real L = std::abs(p2 - p1);

      auto & k_el = *k_it;
      k_el = k;
      k_el *= E * A / L;
    }

    this->getDOFManager().assembleElementalMatricesToMatrix(
        "K", "disp", k_all_el, _segment_2);
  }

  void assembleResidual() {
    this->getDOFManager().assembleToResidual("disp", forces);
    internal_forces.clear();

    this->assembleResidual(_not_ghost);

    this->synchronize(_gst_user_1);

    this->getDOFManager().assembleToResidual("disp", internal_forces, -1.);
  }

  void assembleResidual(const GhostType & ghost_type) {
    Array<Real> forces_internal_el(
        this->mesh.getNbElement(_segment_2, ghost_type), 2);

    auto cit  = this->mesh.getConnectivity(_segment_2, ghost_type).begin(2);
    auto cend = this->mesh.getConnectivity(_segment_2, ghost_type).end(2);

    auto f_it = forces_internal_el.begin(2);

    auto strain_it  = this->strains.begin();
    auto stress_it  = this->stresses.begin();
    auto L_it = this->initial_lengths.begin();

    for (; cit != cend; ++cit, ++f_it, ++strain_it, ++stress_it, ++L_it) {
      const auto & conn = *cit;
      UInt n1 = conn(0);
      UInt n2 = conn(1);

      Real u1 = this->displacement(n1, _x);
      Real u2 = this->displacement(n2, _x);

      *strain_it = (u2 - u1) / *L_it;
      *stress_it = E * *strain_it;
      Real f_n = A * *stress_it;
      Vector<Real> & f = *f_it;

      f(0) = -f_n;
      f(1) = f_n;
    }

    this->getDOFManager().assembleElementalArrayLocalArray(
        forces_internal_el, internal_forces, _segment_2, ghost_type);
  }

  Real getPotentialEnergy() {
    Real res = 0;

    if (!lumped) {
      res = this->mulVectMatVect(this->displacement,
                                 this->getDOFManager().getMatrix("K"),
                                 this->displacement);
    } else {
      auto strain_it  = this->strains.begin();
      auto stress_it  = this->stresses.begin();
      auto strain_end = this->strains.end();
      auto L_it = this->initial_lengths.begin();

      for (; strain_it != strain_end; ++strain_it, ++stress_it, ++L_it) {
        res += *strain_it * *stress_it * A * *L_it;
      }

      StaticCommunicator::getStaticCommunicator().allReduce(res, _so_sum);
    }

    return res / 2.;
  }

  Real getKineticEnergy() {
    Real res = 0;
    if (!lumped) {
      res = this->mulVectMatVect(
          this->velocity, this->getDOFManager().getMatrix("M"), this->velocity);
    } else {
      auto & m = this->getDOFManager().getLumpedMatrix("M");
      auto it = velocity.begin();
      auto end = velocity.end();
      auto m_it = m.begin();

      for (UInt node = 0; it != end; ++it, ++m_it, ++node) {
        if (mesh.isLocalOrMasterNode(node))
          res += *m_it * *it * *it;
      }

      StaticCommunicator::getStaticCommunicator().allReduce(res, _so_sum);
    }

    return res / 2.;
  }

  Real getExternalWorkIncrement() {
    Real res = 0;

    auto it = velocity.begin();
    auto end = velocity.end();
    auto if_it = internal_forces.begin();
    auto ef_it = forces.begin();
    auto b_it = blocked.begin();

    for (UInt node = 0; it != end; ++it, ++if_it, ++ef_it, ++b_it, ++node) {
      if (mesh.isLocalOrMasterNode(node))
        res += (*b_it ? - *if_it : *ef_it) * *it;
    }

    StaticCommunicator::getStaticCommunicator().allReduce(res, _so_sum);

    return res * this->getTimeStep();
  }

  Real mulVectMatVect(const Array<Real> & x, const SparseMatrix & A,
                      const Array<Real> & y) {
    Array<Real> Ay(this->nb_dofs, 1, 0.);
    A.matVecMul(y, Ay);
    Real res = 0.;

    Array<Real>::const_scalar_iterator it = Ay.begin();
    Array<Real>::const_scalar_iterator end = Ay.end();
    Array<Real>::const_scalar_iterator x_it = x.begin();

    for (; it != end; ++it, ++x_it) {
      res += *x_it * *it;
    }

    StaticCommunicator::getStaticCommunicator().allReduce(res, _so_sum);
    return res;
  }

  void predictor() {}
  void corrector() {}

  /* ------------------------------------------------------------------------ */
  UInt getNbData(const Array<Element> & elements, const SynchronizationTag &) const {
    return elements.getSize() * sizeof(Real);
  }

  void packData(CommunicationBuffer & buffer, const Array<Element> & elements,
                const SynchronizationTag & tag) const {
    if(tag == _gst_user_1) {
      for(const auto & el : elements) {
        buffer << this->stresses(el.element);
      }
    }
  }


  void unpackData(CommunicationBuffer & buffer, const Array<Element> & elements,
                  const SynchronizationTag & tag) {
    if(tag == _gst_user_1) {
      auto cit  = this->mesh.getConnectivity(_segment_2, _ghost).begin(2);

      for(const auto & el : elements) {
        Real stress;
        buffer >> stress;

        Real f = A * stress;

        Vector<UInt> conn = cit[el.element];
        this->internal_forces(conn(0), _x) += -f;
        this->internal_forces(conn(1), _x) +=  f;
      }
    }
  }

private:
  Mesh & mesh;
  UInt nb_dofs;
  UInt nb_elements;
  Real E, A, rho;

  bool lumped;
public:
  Array<Real> displacement;
  Array<Real> velocity;
  Array<Real> acceleration;
  Array<bool> blocked;
  Array<Real> forces;
  Array<Real> internal_forces;

  Array<Real> stresses;
  Array<Real> strains;

  Array<Real> initial_lengths;
};

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  initialize(argc, argv);

  UInt prank = StaticCommunicator::getStaticCommunicator().whoAmI();
  UInt global_nb_nodes = 201;
  UInt max_steps = 200;
  Real time_step = 0.001;
  Mesh mesh(1);
  Real F = -9.81;
  bool _explicit = EXPLICIT;

  genMesh(mesh, global_nb_nodes);

  mesh.distribute();

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

  // #if EXPLICIT == true
  //   std::ofstream output("output_dynamic_explicit.csv");
  // #else
  //   std::ofstream output("output_dynamic_implicit.csv");
  // #endif

  if (prank == 0) {
    std::cout << std::scientific;
    std::cout << std::setw(14) << "time"
              << "," << std::setw(14) << "wext"
              << "," << std::setw(14) << "epot"
              << "," << std::setw(14) << "ekin"
              << "," << std::setw(14) << "total" << std::endl;
  }
  Real wext = 0.;

  model.assembleResidual();

  Real epot = model.getPotentialEnergy();
  Real ekin = model.getKineticEnergy();
  Real einit = ekin + epot;
  Real etot = ekin + epot - wext - einit;
  if (prank == 0) {
    std::cout << std::setw(14) << 0.
              << "," << std::setw(14) << wext
              << "," << std::setw(14) << epot
              << "," << std::setw(14) << ekin
              << "," << std::setw(14) << etot
              << std::endl;
  }

#if EXPLICIT == false
  NonLinearSolver & solver =
    model.getDOFManager().getNonLinearSolver("dynamic");
#endif

  for (UInt i = 1; i < max_steps + 1; ++i) {
    model.solveStep();

#if EXPLICIT == false
    UInt nb_iter = solver.get("nb_iterations");
    Real error = solver.get("error");
    bool converged = solver.get("converged");
    if(prank == 0)
      std::cerr << error << " " << nb_iter << " -> " << converged << std::endl;
#endif

    epot = model.getPotentialEnergy();
    ekin = model.getKineticEnergy();
    wext += model.getExternalWorkIncrement();
    etot = ekin + epot - wext - einit;
    if (prank == 0) {
      std::cout << std::setw(14) << time_step * i
                << "," << std::setw(14) << wext
                << "," << std::setw(14) << epot
                << "," << std::setw(14) << ekin
                << "," << std::setw(14) << etot
                << std::endl;
    }

    if (std::abs(etot) > 1e1) {
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
