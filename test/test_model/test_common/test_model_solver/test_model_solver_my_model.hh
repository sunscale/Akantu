/**
 * @file   test_model_solver_my_model.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Apr 13 2016
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Test default dof manager
 *
 *
 * Copyright (©) 2016-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_iterators.hh"
#include "boundary_condition.hh"
#include "communicator.hh"
#include "data_accessor.hh"
#include "dof_manager_default.hh"
#include "element_synchronizer.hh"
#include "mesh.hh"
#include "model_solver.hh"
#include "periodic_node_synchronizer.hh"
#include "solver_vector_default.hh"
#include "sparse_matrix.hh"

/* -------------------------------------------------------------------------- */

namespace akantu {

#ifndef __AKANTU_TEST_MODEL_SOLVER_MY_MODEL_HH__
#define __AKANTU_TEST_MODEL_SOLVER_MY_MODEL_HH__

/**
 *   =\o-----o-----o-> F
 *     |           |
 *     |---- L ----|
 */
class MyModel : public ModelSolver,
                public BoundaryCondition<MyModel>,
                public DataAccessor<Element> {
public:
  MyModel(Real F, Mesh & mesh, bool lumped,
          const ID & dof_manager_type = "default")
      : ModelSolver(mesh, ModelType::_model, "model_solver", 0),
        nb_dofs(mesh.getNbNodes()), nb_elements(mesh.getNbElement(_segment_2)),
        lumped(lumped), E(1.), A(1.), rho(1.), mesh(mesh),
        displacement(nb_dofs, 1, "disp"), velocity(nb_dofs, 1, "velo"),
        acceleration(nb_dofs, 1, "accel"), blocked(nb_dofs, 1, "blocked"),
        forces(nb_dofs, 1, "force_ext"),
        internal_forces(nb_dofs, 1, "force_int"),
        stresses(nb_elements, 1, "stress"), strains(nb_elements, 1, "strain"),
        initial_lengths(nb_elements, 1, "L0") {
    this->initBC(*this, displacement, forces);
    this->initDOFManager(dof_manager_type);
    this->getDOFManager().registerDOFs("disp", displacement, _dst_nodal);
    this->getDOFManager().registerDOFsDerivative("disp", 1, velocity);
    this->getDOFManager().registerDOFsDerivative("disp", 2, acceleration);
    this->getDOFManager().registerBlockedDOFs("disp", blocked);

    displacement.set(0.);
    velocity.set(0.);
    acceleration.set(0.);

    forces.set(0.);
    blocked.set(false);

    UInt global_nb_nodes = mesh.getNbGlobalNodes();
    for (auto && n : arange(nb_dofs)) {
      auto global_id = mesh.getNodeGlobalId(n);
      if (global_id == (global_nb_nodes - 1))
        forces(n, _x) = F;

      if (global_id == 0)
        blocked(n, _x) = true;
    }

    auto cit = this->mesh.getConnectivity(_segment_2).begin(2);
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
    this->registerSynchronizer(
        const_cast<ElementSynchronizer &>(this->mesh.getElementSynchronizer()),
        SynchronizationTag::_user_1);
  }

  void assembleLumpedMass() {
    auto & M = this->getDOFManager().getLumpedMatrix("M");
    M.clear();

    this->assembleLumpedMass(_not_ghost);
    if (this->mesh.getNbElement(_segment_2, _ghost) > 0)
      this->assembleLumpedMass(_ghost);

    is_lumped_mass_assembled = true;
  }

  void assembleLumpedMass(const GhostType & ghost_type) {
    Array<Real> M(nb_dofs, 1, 0.);

    Array<Real> m_all_el(this->mesh.getNbElement(_segment_2, ghost_type), 2);

    for (auto && data :
         zip(make_view(this->mesh.getConnectivity(_segment_2), 2),
             make_view(m_all_el, 2))) {
      const auto & conn = std::get<0>(data);
      auto & m_el = std::get<1>(data);

      UInt n1 = conn(0);
      UInt n2 = conn(1);

      Real p1 = this->mesh.getNodes()(n1, _x);
      Real p2 = this->mesh.getNodes()(n2, _x);

      Real L = std::abs(p2 - p1);

      Real M_n = rho * A * L / 2;
      m_el(0) = m_el(1) = M_n;
    }

    this->getDOFManager().assembleElementalArrayLocalArray(
        m_all_el, M, _segment_2, ghost_type);

    this->getDOFManager().assembleToLumpedMatrix("disp", M, "M");
  }

  void assembleMass() {
    SparseMatrix & M = this->getDOFManager().getMatrix("M");
    M.clear();

    Array<Real> m_all_el(this->nb_elements, 4);

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

    for (auto && data :
         zip(make_view(this->mesh.getConnectivity(_segment_2), 2),
             make_view(m_all_el, 2, 2))) {
      const auto & conn = std::get<0>(data);
      auto & m_el = std::get<1>(data);
      UInt n1 = conn(0);
      UInt n2 = conn(1);

      Real p1 = this->mesh.getNodes()(n1, _x);
      Real p2 = this->mesh.getNodes()(n2, _x);

      Real L = std::abs(p2 - p1);

      m_el = m;
      m_el *= rho * A * L / 6.;
    }

    this->getDOFManager().assembleElementalMatricesToMatrix(
        "M", "disp", m_all_el, _segment_2);

    is_mass_assembled = true;
  }

  MatrixType getMatrixType(const ID &) override { return _symmetric; }

  void assembleMatrix(const ID & matrix_id) override {
    if (matrix_id == "K") {
      if (not is_stiffness_assembled)
        this->assembleStiffness();
    } else if (matrix_id == "M") {
      if (not is_mass_assembled)
        this->assembleMass();
    } else if (matrix_id == "C") {
      // pass, no damping matrix
    } else {
      AKANTU_EXCEPTION("This solver does not know what to do with a matrix "
                       << matrix_id);
    }
  }

  void assembleLumpedMatrix(const ID & matrix_id) override {
    if (matrix_id == "M") {
      if (not is_lumped_mass_assembled)
        this->assembleLumpedMass();
    } else {
      AKANTU_EXCEPTION("This solver does not know what to do with a matrix "
                       << matrix_id);
    }
  }

  void assembleStiffness() {
    SparseMatrix & K = this->getDOFManager().getMatrix("K");
    K.clear();

    Matrix<Real> k(2, 2);
    k(0, 0) = k(1, 1) = 1;
    k(0, 1) = k(1, 0) = -1;

    Array<Real> k_all_el(this->nb_elements, 4);

    auto k_it = k_all_el.begin(2, 2);
    auto cit = this->mesh.getConnectivity(_segment_2).begin(2);
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

    is_stiffness_assembled = true;
  }

  void assembleResidual() override {
    this->getDOFManager().assembleToResidual("disp", forces);
    internal_forces.clear();

    this->assembleResidualInternal(_not_ghost);

    this->synchronize(SynchronizationTag::_user_1);

    this->getDOFManager().assembleToResidual("disp", internal_forces, -1.);
  }

  void assembleResidualInternal(const GhostType & ghost_type) {
    Array<Real> forces_internal_el(
        this->mesh.getNbElement(_segment_2, ghost_type), 2);

    auto cit = this->mesh.getConnectivity(_segment_2, ghost_type).begin(2);
    auto cend = this->mesh.getConnectivity(_segment_2, ghost_type).end(2);

    auto f_it = forces_internal_el.begin(2);

    auto strain_it = this->strains.begin();
    auto stress_it = this->stresses.begin();
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

    if (not lumped) {
      res = this->mulVectMatVect(this->displacement, "K", this->displacement);
    } else {
      auto strain_it = this->strains.begin();
      auto stress_it = this->stresses.begin();
      auto strain_end = this->strains.end();
      auto L_it = this->initial_lengths.begin();

      for (; strain_it != strain_end; ++strain_it, ++stress_it, ++L_it) {
        res += *strain_it * *stress_it * A * *L_it;
      }

      mesh.getCommunicator().allReduce(res, SynchronizerOperation::_sum);
    }

    return res / 2.;
  }

  Real getKineticEnergy() {
    Real res = 0;
    if (not lumped) {
      res = this->mulVectMatVect(this->velocity, "M", this->velocity);
    } else {
      Array<Real> & m = dynamic_cast<SolverVectorDefault &>(
          this->getDOFManager().getLumpedMatrix("M"));
      auto it = velocity.begin();
      auto end = velocity.end();
      auto m_it = m.begin();

      for (UInt node = 0; it != end; ++it, ++m_it, ++node) {
        if (mesh.isLocalOrMasterNode(node))
          res += *m_it * *it * *it;
      }

      mesh.getCommunicator().allReduce(res, SynchronizerOperation::_sum);
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
        res += (*b_it ? -*if_it : *ef_it) * *it;
    }

    mesh.getCommunicator().allReduce(res, SynchronizerOperation::_sum);

    return res * this->getTimeStep();
  }

  Real mulVectMatVect(const Array<Real> & x, const ID & A_id,
                      const Array<Real> & y) {
    Array<Real> Ay(nb_dofs);
    this->getDOFManager().assembleMatMulVectToArray("disp", A_id, y, Ay);

    Real res = 0.;
    for (auto && data : zip(arange(nb_dofs), make_view(Ay), make_view(x))) {
      res += std::get<2>(data) * std::get<1>(data) *
             mesh.isLocalOrMasterNode(std::get<0>(data));
    }

    mesh.getCommunicator().allReduce(res, SynchronizerOperation::_sum);
    return res;
  }

  /* ------------------------------------------------------------------------ */
  UInt getNbData(const Array<Element> & elements,
                 const SynchronizationTag &) const override {
    return elements.size() * sizeof(Real);
  }

  void packData(CommunicationBuffer & buffer, const Array<Element> & elements,
                const SynchronizationTag & tag) const override {
    if (tag == SynchronizationTag::_user_1) {
      for (const auto & el : elements) {
        buffer << this->stresses(el.element);
      }
    }
  }

  void unpackData(CommunicationBuffer & buffer, const Array<Element> & elements,
                  const SynchronizationTag & tag) override {
    if (tag == SynchronizationTag::_user_1) {
      auto cit = this->mesh.getConnectivity(_segment_2, _ghost).begin(2);

      for (const auto & el : elements) {
        Real stress;
        buffer >> stress;

        Real f = A * stress;

        Vector<UInt> conn = cit[el.element];
        this->internal_forces(conn(0), _x) += -f;
        this->internal_forces(conn(1), _x) += f;
      }
    }
  }

  const Mesh & getMesh() const { return mesh; }

  UInt getSpatialDimension() const { return 1; }

  auto & getBlockedDOFs() { return blocked; }

private:
  UInt nb_dofs;
  UInt nb_elements;

  bool lumped;

  bool is_stiffness_assembled{false};
  bool is_mass_assembled{false};
  bool is_lumped_mass_assembled{false};

public:
  Real E, A, rho;

  Mesh & mesh;
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

#endif /* __AKANTU_TEST_MODEL_SOLVER_MY_MODEL_HH__ */

} // namespace akantu
