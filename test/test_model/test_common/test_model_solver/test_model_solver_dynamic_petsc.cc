/**
 * @file   test_model_solver_dynamic.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Apr 13 2016
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Test default dof manager
 *
 * @section LICENSE
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
#include "communicator.hh"
#include "element_group.hh"
#include "mesh.hh"
#include "mesh_accessor.hh"
#include "non_linear_solver.hh"
/* -------------------------------------------------------------------------- */
#include "boundary_condition_functor.hh"
#include "mpi_communicator_data.hh"
/* -------------------------------------------------------------------------- */
#include "dumpable_inline_impl.hh"
#include "dumper_element_partition.hh"
#include "dumper_iohelper_paraview.hh"
/* -------------------------------------------------------------------------- */
#include <fstream>
/* -------------------------------------------------------------------------- */
#include <petscmat.h>
#include <petscsnes.h>
#include <petscvec.h>
/* -------------------------------------------------------------------------- */

#ifndef EXPLICIT
#define EXPLICIT true
#endif

template <typename func>
void CHECK_ERR_CXX(func && func_, PetscErrorCode ierr) {
  if (PetscUnlikely(ierr != 0)) {
    const char * desc;
    PetscErrorMessage(ierr, &desc, nullptr);
    AKANTU_EXCEPTION("Error in PETSc call to \'" << func_ << "\': " << desc);
  }
}

using namespace akantu;

static void genMesh(Mesh & mesh, UInt nb_nodes);

class MyModel {
public:
  MyModel(Real F, Mesh & mesh, bool lumped)
      : nb_dofs(mesh.getNbNodes()), nb_elements(mesh.getNbElement(_segment_2)),
        lumped(lumped), E(1.), A(1.), rho(1.), mesh(mesh),
        displacement(nb_dofs, 1, "disp"), velocity(nb_dofs, 1, "velo"),
        acceleration(nb_dofs, 1, "accel"), blocked(nb_dofs, 1, "blocked"),
        forces(nb_dofs, 1, "force_ext"),
        internal_forces(nb_dofs, 1, "force_int"),
        stresses(nb_elements, 1, "stress"), strains(nb_elements, 1, "strain"),
        initial_lengths(nb_elements, 1, "L0") {

    auto n_global = mesh.getNbGlobalNodes();
    int n_local = 0;

    std::vector<PetscInt> nodes_global_ids(nb_dofs);
    for (auto && data : enumerate(nodes_global_ids)) {
      auto n = std::get<0>(data);
      n_local += mesh.isLocalOrMasterNode(n);
      std::get<1>(data) = mesh.getNodeGlobalId(n);
    }

    mpi_comm = dynamic_cast<MPICommunicatorData &>(
                   mesh.getCommunicator().getCommunicatorData())
                   .getMPICommunicator();

    MeshAccessor mesh_accessor(mesh);

    ierr = ISLocalToGlobalMappingCreate(
        mpi_comm, 1, mesh.getNbNodes(), nodes_global_ids.data(),
        PETSC_COPY_VALUES, &petsc_local_to_global);
    CHECK_ERR_CXX("ISLocalToGlobalMappingCreate", ierr);

    auto setName = [](auto && Obj, auto && name) {
      PetscObjectSetName(reinterpret_cast<PetscObject>(Obj), name);
    };

    ierr = VecCreate(mpi_comm, &rhs);
    ierr = VecSetSizes(rhs, n_local, n_global);
    ierr = VecSetFromOptions(rhs);
    ierr = VecSetLocalToGlobalMapping(rhs, petsc_local_to_global);
    setName(rhs, "rhs");

    ierr = VecDuplicate(rhs, &x);
    ierr = VecDuplicate(rhs, &x_save);
    ierr = VecDuplicate(rhs, &dx);
    ierr = VecDuplicate(rhs, &f_int);
    ierr = VecDuplicate(rhs, &f_dirichlet);
    setName(x, "x");
    setName(x_save, "x save");
    setName(dx, "dx");
    setName(f_int, "f_int");
    setName(f_dirichlet, "f_dirichlet");

    ierr = MatCreate(mpi_comm, &M);
    ierr = MatSetSizes(M, n_local, n_local, n_global, n_global);
    ierr = MatSetFromOptions(M);
    ierr = MatSetOption(M, MAT_SYMMETRIC, PETSC_TRUE);
    ierr = MatSetOption(M, MAT_ROW_ORIENTED, PETSC_TRUE);
    ierr = MatSetUp(M);
    ierr = MatSetLocalToGlobalMapping(M, petsc_local_to_global,
                                      petsc_local_to_global);
    setName(M, "M");

    assembleMass();

    ierr = MatDuplicate(M, MAT_DO_NOT_COPY_VALUES, &K);
    setName(K, "K");
    ierr = MatDuplicate(M, MAT_DO_NOT_COPY_VALUES, &J);
    setName(J, "J");

    ierr = SNESCreate(mpi_comm, &snes);
    ierr = SNESSetFromOptions(snes);
    ierr = SNESSetFunction(snes, rhs, MyModel::FormFunction, this);
    ierr = SNESSetJacobian(snes, J, J, MyModel::FormJacobian, this);

    PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_INDEX);

    displacement.set(0.);
    velocity.set(0.);
    acceleration.set(0.);

    forces.set(0.);
    blocked.set(false);
    blocked(0, 0) = true;
    blocked(nb_dofs - 1, 0) = true;
    displacement(0, 0) = 0;
    displacement(nb_dofs - 1, 0) = 1;

    for (auto && data :
         zip(make_view(this->mesh.getConnectivity(_segment_2), 2),
             make_view(this->initial_lengths))) {
      const auto & conn = std::get<0>(data);
      auto & L = std::get<1>(data);

      auto p1 = this->mesh.getNodes()(conn(0), _x);
      auto p2 = this->mesh.getNodes()(conn(1), _x);

      L = std::abs(p2 - p1);
    }
  }

  // static PetscErrorCode SNESMonitor(SNES snes,PetscInt its,PetscReal
  // fnorm,void *ctx) {
  //   auto & _this = *reinterpret_cast<MyModel *>(ctx);
  //   //SNESMonitorDefault(snes, its, fnorm, PETSC_VIEWER_STDOUT_WORLD);
  // }

  static PetscErrorCode FormFunction(SNES /*snes*/, Vec /*dx*/, Vec /*f*/,
                                     void * ctx) {
    auto & _this = *reinterpret_cast<MyModel *>(ctx);
    _this.assembleResidual();
    return 0;
  }

  static PetscErrorCode FormJacobian(SNES /*snes*/, Vec /*dx*/, Mat /*J*/,
                                     Mat /*P*/, void * ctx) {
    auto & _this = *reinterpret_cast<MyModel *>(ctx);
    _this.assembleJacobian();
    return 0;
  }

  ~MyModel() {
    ierr = MatDestroy(&M);
    ierr = MatDestroy(&K);
    ierr = MatDestroy(&J);

    ierr = VecDestroy(&rhs);
    ierr = VecDestroy(&x);
    ierr = VecDestroy(&dx);
    ierr = VecDestroy(&x_save);
    ierr = VecDestroy(&f_int);

    PetscFinalize();
  }

  void solveStep() {
    std::cout << "solveStep" << std::endl;
    copy(x_save, displacement);

    ierr = SNESSolve(snes, NULL, dx);
    CHECK_ERR_CXX("SNESSolve", ierr);

    setSolutionToDisplacement();
    assembleResidual();
  }

  void applyBC() {
    std::vector<PetscInt> rows;
    for (auto && data : enumerate(blocked)) {
      if (std::get<1>(data)) {
        rows.push_back(std::get<0>(data));
      }
    }

    copy(x, displacement);
    ierr = MatZeroRowsColumnsLocal(J, rows.size(), rows.data(), 1., x,
                                   f_dirichlet);
    VecView(f_dirichlet, PETSC_VIEWER_STDOUT_WORLD);
    CHECK_ERR_CXX("MatZeroRowsColumnsLocal", ierr);
  }

  void setSolutionToDisplacement() {
    std::cout << "setSolutionToDisplacement" << std::endl;
    ierr = VecWAXPY(x, 1, x_save, dx);
    copy(displacement, x);
  }

  void assembleJacobian() {
    std::cout << "assembleJacobian" << std::endl;
    setSolutionToDisplacement();

    assembleStiffness();

    ierr = MatZeroEntries(J);
    CHECK_ERR_CXX("MatZeroEntries", ierr);

    ierr = MatAXPY(J, 1., K, SAME_NONZERO_PATTERN);
    CHECK_ERR_CXX("MatAXPY", ierr);

    MatView(J, PETSC_VIEWER_STDOUT_WORLD);
    applyBC();
    MatView(J, PETSC_VIEWER_STDOUT_WORLD);
  }

  void assembleMass() {
    std::cout << "assembleMass" << std::endl;
    ierr = MatZeroEntries(M);
    CHECK_ERR_CXX("MatZeroEntries", ierr);

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

      Vector<Int> conn_int(conn.size());
      for (auto && data : zip(conn_int, conn)) {
        std::get<0>(data) = std::get<1>(data);
      }

      ierr = MatSetValuesLocal(M, conn_int.size(), conn_int.storage(),
                               conn_int.size(), conn_int.storage(), m.storage(),
                               ADD_VALUES);
    }

    ierr = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
    ierr = MatSetOption(M, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);

    PetscViewer viewer;
    ierr = PetscViewerASCIIOpen(mpi_comm, "M.mtx", &viewer);
    PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATRIXMARKET);
    ierr = MatView(M, viewer);
    PetscViewerPopFormat(viewer);
    ierr = PetscViewerDestroy(&viewer);
    // this->getDOFManager().assembleElementalMatricesToMatrix(
    //   "M", "disp", m_all_el, _segment_2);

    is_mass_assembled = true;
  }

  // MatrixType getMatrixType(const ID &) { return _symmetric; }

  // void assembleMatrix(const ID & matrix_id) {
  //   if (matrix_id == "K") {
  //     if (not is_stiffness_assembled)
  //       this->assembleStiffness();
  //   } else if (matrix_id == "M") {
  //     if (not is_mass_assembled)
  //       this->assembleMass();
  //   } else if (matrix_id == "C") {
  //     // pass, no damping matrix
  //   } else {
  //     AKANTU_EXCEPTION("This solver does not know what to do with a matrix "
  //                      << matrix_id);
  //   }
  // }

  void assembleLumpedMatrix(const ID & matrix_id) {
    std::cout << "assembleLumpedMatrix" << std::endl;
    AKANTU_EXCEPTION("This solver does not know what to do with a matrix "
                     << matrix_id);
  }

  void assembleStiffness() {
    std::cout << "assembleStiffness" << std::endl;
    // SparseMatrix & K = this->getDOFManager().getMatrix("K");
    // K.clear();
    ierr = MatZeroEntries(K);
    CHECK_ERR_CXX("MatZeroEntries", ierr);

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

      Vector<Int> conn_int(conn.size());
      for (auto && data : zip(conn_int, conn)) {
        std::get<0>(data) = std::get<1>(data);
      }

      ierr = MatSetValuesLocal(K, conn_int.size(), conn_int.storage(),
                               conn_int.size(), conn_int.storage(),
                               k_el.storage(), ADD_VALUES);
    }

    ierr = MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
    CHECK_ERR_CXX("MatAssemblyBegin", ierr);

    ierr = MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
    CHECK_ERR_CXX("MatAssemblyEnd", ierr);

    ierr = MatSetOption(K, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
    CHECK_ERR_CXX("MatSetOption", ierr);

    PetscViewer viewer;
    ierr = PetscViewerASCIIOpen(mpi_comm, "K.mtx", &viewer);
    CHECK_ERR_CXX("PetscViewerASCIIOpen", ierr);

    PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATRIXMARKET);
    ierr = MatView(K, viewer);
    CHECK_ERR_CXX("MatView", ierr);

    PetscViewerPopFormat(viewer);
    ierr = PetscViewerDestroy(&viewer);
    CHECK_ERR_CXX("PetscViewerDestroy", ierr);

    // this->getDOFManager().assembleElementalMatricesToMatrix(
    //     "K", "disp", k_all_el, _segment_2);

    is_stiffness_assembled = true;
  }

  void copy(Array<Real> & y, Vec x) {
    std::cout << "copy <-" << std::endl;
    const PetscScalar * x_local;
    ierr = VecGetArrayRead(x, &x_local);

    for (auto && data : zip(y, range(x_local + 0, x_local + y.size()))) {
      std::get<0>(data) = std::get<1>(data);
    }
    ierr = VecRestoreArrayRead(x, &x_local);

    // VecView(x, PETSC_VIEWER_STDOUT_WORLD);
    // std::cout << y.getID() << " " << Vector<Real>(y.storage(), y.size())
    //           << std::endl;
  }

  void print(const Array<Real> & x) const {
    std::cout << x.getID() << " " << Vector<Real>(x.storage(), x.size())
              << std::endl;
  }

  void copy(Vec x, const Array<Real> & y) {
    std::cout << "copy ->" << std::endl;
    PetscScalar * x_local;
    ierr = VecGetArray(x, &x_local);

    for (auto && data : zip(y, range(x_local + 0, x_local + y.size()))) {
      std::get<1>(data) = std::get<0>(data);
    }
    ierr = VecRestoreArray(x, &x_local);

    // std::cout << y.getID() << " " << Vector<Real>(y.storage(), y.size())
    //           << std::endl;
    // VecView(x, PETSC_VIEWER_STDOUT_WORLD);
  }

  void assembleResidual() {
    std::cout << "assembleResidual" << std::endl;
    //   this->getDOFManager().assembleToResidual("disp", forces);
    setSolutionToDisplacement();
    copy(rhs, forces);
    // VecAXPY(rhs, -1., f_dirichlet);

    print(displacement);
    this->assembleResidual(_not_ghost);
    // this->synchronize(SynchronizationTag::_user_1);

    // this->getDOFManager().assembleToResidual("disp", internal_forces, -1.);
    VecAXPY(rhs, 1., f_int);

    for (auto && data : enumerate(blocked)) {
      if (std::get<1>(data)) {
        VecSetValueLocal(rhs, std::get<0>(data), 0., INSERT_VALUES);
      }
    }
    VecAssemblyBegin(rhs);
    VecAssemblyEnd(rhs);

    VecView(rhs, PETSC_VIEWER_STDOUT_WORLD);
  }

  void assembleResidual(const GhostType & ghost_type) {
    std::cout << "assembleResidual" << std::endl;
    VecZeroEntries(f_int);

    auto cit = this->mesh.getConnectivity(_segment_2, ghost_type).begin(2);
    auto cend = this->mesh.getConnectivity(_segment_2, ghost_type).end(2);

    auto strain_it = this->strains.begin();
    auto stress_it = this->stresses.begin();
    auto L_it = this->initial_lengths.begin();

    for (; cit != cend; ++cit, ++strain_it, ++stress_it, ++L_it) {
      const auto & conn = *cit;
      UInt n1 = conn(0);
      UInt n2 = conn(1);

      Real u1 = this->displacement(n1, _x);
      Real u2 = this->displacement(n2, _x);

      *strain_it = (u2 - u1) / *L_it;
      *stress_it = E * *strain_it;
      Real f_n = A * *stress_it;

      std::cout << n1 << "[" << u1 << "]"
                << " <-> " << n2 << "[" << u2 << "]"
                << " : " << f_n << std::endl;

      ierr = VecSetValueLocal(f_int, n1, -f_n, ADD_VALUES);
      ierr = VecSetValueLocal(f_int, n2, f_n, ADD_VALUES);
    }

    ierr = VecAssemblyBegin(f_int);
    ierr = VecAssemblyEnd(f_int);
    // this->getDOFManager().assembleElementalArrayLocalArray(
    //     forces_internal_el, internal_forces, _segment_2, ghost_type);
  }

  Real getPotentialEnergy() {
    std::cout << "getPotentialEnergy" << std::endl;
    copy(x, displacement);
    Vec Ax;

    ierr = VecDuplicate(x, &Ax);
    ierr = MatMult(K, x, Ax);
    PetscScalar res;
    ierr = VecDot(x, Ax, &res);

    return res / 2.;
  }

  Real getKineticEnergy() {
    std::cout << "getKineticEnergy" << std::endl;
    return 0;
  }

  // Real getExternalWorkIncrement() {
  //   Real res = 0;

  //   auto it = velocity.begin();
  //   auto end = velocity.end();
  //   auto if_it = internal_forces.begin();
  //   auto ef_it = forces.begin();
  //   auto b_it = blocked.begin();

  //   for (UInt node = 0; it != end; ++it, ++if_it, ++ef_it, ++b_it, ++node) {
  //     if (mesh.isLocalOrMasterNode(node))
  //       res += (*b_it ? -*if_it : *ef_it) * *it;
  //   }

  //   mesh.getCommunicator().allReduce(res, SynchronizerOperation::_sum);

  //   return res * this->getTimeStep();
  // }

  // void predictor() {}
  // void corrector() {}

  // /* ------------------------------------------------------------------------
  // */ UInt getNbData(const Array<Element> & elements,
  //                const SynchronizationTag &) const {
  //   return elements.size() * sizeof(Real);
  // }

  // void packData(CommunicationBuffer & buffer, const Array<Element> &
  // elements,
  //               const SynchronizationTag & tag) const {
  //   if (tag == SynchronizationTag::_user_1) {
  //     for (const auto & el : elements) {
  //       buffer << this->stresses(el.element);
  //     }
  //   }
  // }

  // void unpackData(CommunicationBuffer & buffer, const Array<Element> &
  // elements,
  //                 const SynchronizationTag & tag) {
  //   if (tag == SynchronizationTag::_user_1) {
  //     auto cit = this->mesh.getConnectivity(_segment_2, _ghost).begin(2);

  //     for (const auto & el : elements) {
  //       Real stress;
  //       buffer >> stress;

  //       Real f = A * stress;

  //       Vector<UInt> conn = cit[el.element];
  //       this->internal_forces(conn(0), _x) += -f;
  //       this->internal_forces(conn(1), _x) += f;
  //     }
  //   }
  // }

  Real getExternalWorkIncrement() {
    std::cout << "getExternalWorkIncrement" << std::endl;
    return 0.;
  }

  template <class Functor> void applyBC(Functor && func, const ID & group_id) {
    auto & group = mesh.getElementGroup(group_id).getNodeGroup().getNodes();

    auto blocked_dofs = make_view(blocked, 1).begin();
    auto disps = make_view(displacement, 1).begin();
    auto poss = make_view(mesh.getNodes(), 1).begin();
    for (auto && node : group) {
      auto disp = Vector<Real>(disps[node]);
      auto pos = Vector<Real>(poss[node]);
      auto flags = Vector<bool>(blocked_dofs[node]);
      func(node, flags, disp, pos);
    }
  }

  const Mesh & getMesh() const { return mesh; }

  UInt getSpatialDimension() const { return 1; }

  auto & getBlockedDOFs() { return blocked; }

  void setTimeStep(Real dt) {
    std::cout << "setTimeStep" << std::endl;
    this->dt = dt;
  }

private:
  PetscErrorCode ierr{0};
  MPI_Comm mpi_comm;
  ISLocalToGlobalMapping petsc_local_to_global;

  UInt nb_dofs;
  UInt nb_elements;

  bool lumped;

  bool is_stiffness_assembled{false};
  bool is_mass_assembled{false};
  bool is_lumped_mass_assembled{false};

  Mat K{nullptr}, J{nullptr}, M{nullptr};
  Vec rhs{nullptr}, x{nullptr}, x_save{nullptr}, dx{nullptr}, f_int{nullptr},
      f_dirichlet{nullptr};

  SNES snes;

  Real dt{0};
  Array<Real> save_displacement;

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

/* -------------------------------------------------------------------------- */
class Sinusoidal : public BC::Dirichlet::DirichletFunctor {
public:
  Sinusoidal(MyModel & model, Real amplitude, Real pulse_width, Real t)
      : model(model), A(amplitude), k(2 * M_PI / pulse_width),
        t(t), v{std::sqrt(model.E / model.rho)} {}

  void operator()(UInt n, Vector<bool> & /*flags*/, Vector<Real> & disp,
                  const Vector<Real> & coord) const {
    auto x = coord(_x);
    model.velocity(n, _x) = k * v * A * sin(k * (x - v * t));
    disp(_x) = A * cos(k * (x - v * t));
  }

private:
  MyModel & model;
  Real A{1.};
  Real k{2 * M_PI};
  Real t{1.};
  Real v{1.};
};

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  initialize(argc, argv);

  PetscInitialize(&argc, &argv, nullptr, nullptr);

  UInt prank = Communicator::getStaticCommunicator().whoAmI();
  UInt global_nb_nodes = 3;
  UInt max_steps = 400;
  Real time_step = 0.001;
  Mesh mesh(1);
  Real F = -9.81;
  bool _explicit = EXPLICIT;
  // const Real pulse_width = 0.2;
  const Real A = 0.01;

  if (prank == 0)
    genMesh(mesh, global_nb_nodes);

  mesh.distribute();

  // mesh.makePeriodic(_x);

  MyModel model(F, mesh, _explicit);

  //  model.forces.clear();
  //  model.blocked.clear();

  // model.applyBC(Sinusoidal(model, A, pulse_width, 0.), "all");
  // model.applyBC(BC::Dirichlet::FlagOnly(_x), "border");

  // if (!_explicit) {
  //   model.getNewSolver("dynamic", TimeStepSolverType::_dynamic,
  //                      NonLinearSolverType::_newton_raphson);
  //   model.setIntegrationScheme("dynamic", "disp",
  //                              IntegrationSchemeType::_trapezoidal_rule_2,
  //                              IntegrationScheme::_displacement);
  // } else {
  //   model.getNewSolver("dynamic", TimeStepSolverType::_dynamic_lumped,
  //                      NonLinearSolverType::_lumped);
  //   model.setIntegrationScheme("dynamic", "disp",
  //                              IntegrationSchemeType::_central_difference,
  //                              IntegrationScheme::_acceleration);
  // }

  model.setTimeStep(time_step);

  if (prank == 0) {
    std::cout << std::scientific;
    std::cout << std::setw(14) << "time"
              << "," << std::setw(14) << "wext"
              << "," << std::setw(14) << "epot"
              << "," << std::setw(14) << "ekin"
              << "," << std::setw(14) << "total"
              << "," << std::setw(14) << "max_disp"
              << "," << std::setw(14) << "min_disp" << std::endl;
  }
  Real wext = 0.;

  // model.getDOFManager().clearResidual();
  // model.assembleResidual();

  Real epot = 0; // model.getPotentialEnergy();
  Real ekin = 0; // model.getKineticEnergy();
  Real einit = ekin + epot;
  Real etot = ekin + epot - wext - einit;

  Real max_disp = 0., min_disp = 0.;
  for (auto && disp : model.displacement) {
    max_disp = std::max(max_disp, disp);
    min_disp = std::min(min_disp, disp);
  }

  if (prank == 0) {
    std::cout << std::setw(14) << 0. << "," << std::setw(14) << wext << ","
              << std::setw(14) << epot << "," << std::setw(14) << ekin << ","
              << std::setw(14) << etot << "," << std::setw(14) << max_disp
              << "," << std::setw(14) << min_disp << std::endl;
  }

  // #if EXPLICIT == false
  //   NonLinearSolver & solver =
  //       model.getDOFManager().getNonLinearSolver("dynamic");

  //   solver.set("max_iterations", 20);
  // #endif

  auto * dumper = new DumperParaview("dynamic", "./paraview");
  mesh.registerExternalDumper(*dumper, "dynamic", true);
  mesh.addDumpMesh(mesh);

  mesh.addDumpFieldExternalToDumper("dynamic", "displacement",
                                    model.displacement);
  mesh.addDumpFieldExternalToDumper("dynamic", "velocity", model.velocity);
  mesh.addDumpFieldExternalToDumper("dynamic", "forces", model.forces);
  mesh.addDumpFieldExternalToDumper("dynamic", "acceleration",
                                    model.acceleration);

  mesh.dump();
  max_steps = 1;
  for (UInt i = 1; i < max_steps + 1; ++i) {
    // model.applyBC(Sinusoidal(model, A, pulse_width, time_step * (i - 1)),
    //             "border");

    model.solveStep();
    mesh.dump();

    epot = model.getPotentialEnergy();
    ekin = model.getKineticEnergy();
    wext += model.getExternalWorkIncrement();
    etot = ekin + epot - wext - einit;

    Real max_disp = 0., min_disp = 0.;
    for (auto && disp : model.displacement) {
      max_disp = std::max(max_disp, disp);
      min_disp = std::min(min_disp, disp);
    }

    if (prank == 0) {
      std::cout << std::setw(14) << time_step * i << "," << std::setw(14)
                << wext << "," << std::setw(14) << epot << "," << std::setw(14)
                << ekin << "," << std::setw(14) << etot << "," << std::setw(14)
                << max_disp << "," << std::setw(14) << min_disp << std::endl;
    }
  }

  // output.close();
  //  finalize();
  // PetscFinalize();

  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
void genMesh(Mesh & mesh, UInt nb_nodes) {
  MeshAccessor mesh_accessor(mesh);
  Array<Real> & nodes = mesh_accessor.getNodes();
  Array<UInt> & conn = mesh_accessor.getConnectivity(_segment_2);

  nodes.resize(nb_nodes);

  // auto & all = mesh.createNodeGroup("all_nodes");

  for (UInt n = 0; n < nb_nodes; ++n) {
    nodes(n, _x) = n * (1. / (nb_nodes - 1));
    // all.add(n);
  }

  // mesh.createElementGroupFromNodeGroup("all", "all_nodes");

  conn.resize(nb_nodes - 1);
  for (UInt n = 0; n < nb_nodes - 1; ++n) {
    conn(n, 0) = n;
    conn(n, 1) = n + 1;
  }

  // Array<UInt> & conn_points = mesh_accessor.getConnectivity(_point_1);
  // conn_points.resize(2);

  // conn_points(0, 0) = 0;
  // conn_points(1, 0) = nb_nodes - 1;

  // auto & border = mesh.createElementGroup("border", 0);
  // border.add({_point_1, 0, _not_ghost}, true);
  // border.add({_point_1, 1, _not_ghost}, true);

  mesh_accessor.makeReady();
}
