/**
 * @file   solid_phase_coupler.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Fri Sep 28 2018
 * @date last modification: Fri Sep 28 2018
 *
 * @brief  class for coupling of solid mechancis and phasefield model
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
#include "boundary_condition.hh"
#include "data_accessor.hh"
#include "fe_engine.hh"
#include "material.hh"
#include "material_phasefield.hh"
#include "model.hh"
#include "phase_field_model.hh"
#include "solid_mechanics_model.hh"
#include "sparse_matrix.hh"
#include "time_step_solver.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_COUPLER_SOLID_PHASEFIELD_HH__
#define __AKANTU_COUPLER_SOLID_PHASEFIELD_HH__

/* ------------------------------------------------------------------------ */
/* Coupling : Solid Mechanics / PhaseField                                  */
/* ------------------------------------------------------------------------ */
namespace akantu {
template <ElementKind kind, class IntegrationOrderFunctor>
class IntegratorGauss;
template <ElementKind kind> class ShapeLagrange;
class DOFManager;
} // namespace akantu

namespace akantu {

class CouplerSolidPhaseField
    : public Model,
      public DataAccessor<Element>,
      public DataAccessor<UInt>,
      public BoundaryCondition<CouplerSolidPhaseField> {

  /* ------------------------------------------------------------------------ */
  /*  Constructor/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

  using MyFEEngineType = FEEngineTemplate<IntegratorGauss, ShapeLagrange>;

public:
  CouplerSolidPhaseField(
      Mesh & mesh, UInt spatial_dimension = _all_dimensions,
      const ID & id = "coupler_solid_phasefield",
      const MemoryID & memory_id = 0,
      const ModelType model_type = ModelType::_coupler_solid_phasefield);

  ~CouplerSolidPhaseField() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  /// initialize the complete model
  void initFullImpl(const ModelOptions & options) override;

  /// initialize the modelType
  void initModel() override;

  /// get some default values for derived classes
  std::tuple<ID, TimeStepSolverType>
  getDefaultSolverID(const AnalysisMethod & method) override;

  /* ------------------------------------------------------------------------ */
  /* Solver Interface                                                         */
  /* ------------------------------------------------------------------------ */
public:
  /// assembles the contact stiffness matrix
  virtual void assembleStiffnessMatrix();

  /// assembles the contant internal forces
  virtual void assembleInternalForces();

public:
  /// computes damage on quad points for solid mechanics model from
  /// damage array from phasefield model
  void computeDamageOnQuadPoints(const GhostType &);

  /// computes strain on quadrature points for phasefield model from
  /// displacement gradient from solid mechanics model
  void computeStrainOnQuadPoints(const GhostType & ghost_type);

  /// solve the coupled model
  void solve();

private:
  /// computes small strain from displacement gradient
  void gradUToEpsilon(const Matrix<Real> & grad_u, Matrix<Real> & epsilon);

  /// test the convergence criteria
  bool checkConvergence(Array<Real> &, Array<Real> &, Array<Real> &,
                        Array<Real> &);

protected:
  /// callback for the solver, this adds f_{ext} - f_{int} to the residual
  void assembleResidual() override;

  /// callback for the solver, this adds f_{ext} or  f_{int} to the residual
  void assembleResidual(const ID & residual_part) override;
  bool canSplitResidual() override { return true; }

  /// get the type of matrix needed
  MatrixType getMatrixType(const ID & matrix_id) override;

  /// callback for the solver, this assembles different matrices
  void assembleMatrix(const ID & matrix_id) override;

  /// callback for the solver, this assembles the stiffness matrix
  void assembleLumpedMatrix(const ID & matrix_id) override;

  /// callback for the model to instantiate the matricess when needed
  void initSolver(TimeStepSolverType, NonLinearSolverType) override;

  /// solve the coupled model
  //void solveStep(const ID & solver_id = "") override;

  /// solve a step using a given pre instantiated time step solver and
  /// non linear solver with a user defined callback instead of the
  /// model itself /!\ This can mess up everything
  //void solveStep(SolverCallback & callback, const ID & solver_id = "") override;

  /* ------------------------------------------------------------------------ */
  /* Mass matrix for solid mechanics model                                    */
  /* ------------------------------------------------------------------------ */
public:
  /// assemble the lumped mass matrix
  void assembleMassLumped();

  /// assemble the mass matrix for consistent mass resolutions
  void assembleMass();

protected:
  /// assemble the lumped mass matrix for local and ghost elements
  void assembleMassLumped(GhostType ghost_type);

  /// assemble the mass matrix for either _ghost or _not_ghost elements
  void assembleMass(GhostType ghost_type);

protected:
  /* --------------------------------------------------------------------------
   */
  TimeStepSolverType getDefaultSolverType() const override;
  /* --------------------------------------------------------------------------
   */
  ModelSolverOptions
  getDefaultSolverOptions(const TimeStepSolverType & type) const;

public:
  bool isDefaultSolverExplicit() { return method == _explicit_lumped_mass; }

  /* ------------------------------------------------------------------------ */
public:
  // DataAccessor<Element>
  UInt getNbData(const Array<Element> &,
                 const SynchronizationTag &) const override;
  void packData(CommunicationBuffer &, const Array<Element> &,
                const SynchronizationTag &) const override;
  void unpackData(CommunicationBuffer &, const Array<Element> &,
                  const SynchronizationTag &) override;

  UInt getNbData(const Array<UInt> & indexes,
                 const SynchronizationTag & tag) const override {}

  void packData(CommunicationBuffer & buffer, const Array<UInt> & dofs,
                const SynchronizationTag & tag) const override{}

  void unpackData(CommunicationBuffer & buffer, const Array<UInt> & dofs,
                  const SynchronizationTag & tag) override {}


  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  FEEngine & getFEEngineBoundary(const ID & name = "") override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// return the dimension of the system space
  AKANTU_GET_MACRO(SpatialDimension, Model::spatial_dimension, UInt);

  /// get the solid mechanics model
  AKANTU_GET_MACRO(SolidMechanicsModel, *solid, SolidMechanicsModel &);

  /// get the contact mechanics model
  AKANTU_GET_MACRO(PhaseFieldModel, *phasefield, PhaseFieldModel &);

  /* ------------------------------------------------------------------------ */
  /* Dumpable interface                                                       */
  /* ------------------------------------------------------------------------ */
public:
  std::shared_ptr<dumper::Field>
  createNodalFieldReal(const std::string & field_name,
                       const std::string & group_name,
                       bool padding_flag) override;

  std::shared_ptr<dumper::Field>
  createNodalFieldBool(const std::string & field_name,
                       const std::string & group_name,
                       bool padding_flag) override;

  std::shared_ptr<dumper::Field>
  createElementalField(const std::string & field_name,
                       const std::string & group_name, bool padding_flag,
                       const UInt & spatial_dimension,
                       const ElementKind & kind) override;

  virtual void dump(const std::string & dumper_name);

  virtual void dump(const std::string & dumper_name, UInt step);

  virtual void dump(const std::string & dumper_name, Real time, UInt step);

  void dump() override;

  virtual void dump(UInt step);

  virtual void dump(Real time, UInt step);

  /* ------------------------------------------------------------------------ */
  /* Members                                                                  */
  /* ------------------------------------------------------------------------ */
private:
  /// solid mechanics model
  SolidMechanicsModel * solid{nullptr};

  /// phasefield model
  PhaseFieldModel * phasefield{nullptr};

  Array<Real> * displacement{nullptr};

  ///
  Array<Real> * displacement_increment{nullptr};

  /// external forces array
  Array<Real> * external_force{nullptr};
};

} // namespace akantu

#endif /* __AKANTU_COUPLER_SOLID_PHASEFIELD_HH__ */
