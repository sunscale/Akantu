/**
 * @file   phase_field_model.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Sun Jul 30 2018
 * @date last modification: Mon Feb 05 2018
 *
 * @brief  Model class for Phase Field problem
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
#include "model.hh"
/* -------------------------------------------------------------------------- */
#include <array>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_PHASE_FIELD_MODEL_HH__
#define __AKANTU_PHASE_FIELD_MODEL_HH__

namespace akantu {
class PhaseField;
class PhaseFieldSelector;
template <ElementKind kind, class IntegrationOrderFuntor>
class IntegratorGauss;
template <ElementKind kind> class ShapeLagrange;
} // namespace akantu

/* -------------------------------------------------------------------------- */
namespace akantu {

/* -------------------------------------------------------------------------- */
class PhaseFieldModel
  : public Model,
    public DataAccessor<Element>,
    public DataAccessor<UInt>,
    public BoundaryCondition<PhaseFieldModel> {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using FEEngineType = FEEngineTemplate<IntegratorGauss, ShapeLagrange>;

  PhaseFieldModel(Mesh & mesh, UInt spatial_dimension = _all_dimensions,
                  const ID & id = "phase_field_model",
		  const ModelType model_type = ModelType::_phase_field_model);

  ~PhaseFieldModel() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  /// generic function to initialize everything ready for explicit dynamics
  void initFullImpl(const ModelOptions & options) override;

  /// initialize all internal array for phasefields
  void initPhaseFields();

  /// allocate all vectors
  void initSolver(TimeStepSolverType, NonLinearSolverType) override;

  /// initialize the model
  void initModel() override;

  /// predictor
  void predictor() override;

  /// corrector
  void corrector() override;

  /// compute the heat flux
  void assembleResidual() override;

  /// get the type of matrix needed
  MatrixType getMatrixType(const ID &) override;

  /// callback to assemble a Matrix
  void assembleMatrix(const ID &) override;

  /// callback to assemble a lumped Matrix
  void assembleLumpedMatrix(const ID &) override;

  std::tuple<ID, TimeStepSolverType>
  getDefaultSolverID(const AnalysisMethod & method) override;

  ModelSolverOptions
  getDefaultSolverOptions(const TimeStepSolverType & type) const override;

  /// function to print the containt of the class
  void printself(std::ostream & stream, int indent = 0) const override;

 /* ------------------------------------------------------------------------ */
  /* Materials (phase_field_model.cc)                            */
  /* ------------------------------------------------------------------------ */
public:
  /// register an empty phasefield of a given type
  PhaseField & registerNewPhaseField(const ID & mat_name, const ID & mat_type,
				   const ID & opt_param);

  /// reassigns phasefields depending on the phasefield selector
  void reassignPhaseField();

protected:
  /// register a phasefield in the dynamic database
  PhaseField & registerNewPhaseField(const ParserSection & phase_section);

  /// read the phasefield files to instantiate all the phasefields
  void instantiatePhaseFields();

  /// set the element_id_by_phasefield and add the elements to the good phasefields
  void
  assignPhaseFieldToElements(const ElementTypeMapArray<UInt> * filter = nullptr);

 
  /* ------------------------------------------------------------------------ */
  /* Methods for static                                                       */
  /* ------------------------------------------------------------------------ */
public:
  /// assembles the phasefield stiffness matrix
  virtual void assembleStiffnessMatrix();

  /// compute the internal forces
  virtual void assembleInternalForces();

  // compute the internal forces
  void assembleInternalForces(const GhostType & ghost_type);
  
  /* ------------------------------------------------------------------------ */
  /* Methods for dynamic                                                      */
  /* ------------------------------------------------------------------------ */
public:
  /// set the stable timestep
  void setTimeStep(Real time_step, const ID & solver_id = "") override;

protected:
  /// callback for the solver, this is called at beginning of solve
  void beforeSolveStep() override;
  /// callback for the solver, this is called at end of solve
  void afterSolveStep(bool converged = true) override;

  /* ------------------------------------------------------------------------ */
  /* Data Accessor inherited members                                          */
  /* ------------------------------------------------------------------------ */
public:
  UInt getNbData(const Array<Element> & elements,
                 const SynchronizationTag & tag) const override;

  void packData(CommunicationBuffer & buffer, const Array<Element> & elements,
                const SynchronizationTag & tag) const override;

  void unpackData(CommunicationBuffer & buffer, const Array<Element> & elements,
                  const SynchronizationTag & tag) override;

  UInt getNbData(const Array<UInt> & indexes,
                 const SynchronizationTag & tag) const override;

  void packData(CommunicationBuffer & buffer, const Array<UInt> & dofs,
                const SynchronizationTag & tag) const override;

  void unpackData(CommunicationBuffer & buffer, const Array<UInt> & dofs,
                  const SynchronizationTag & tag) override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// return the dimension of the system space
  AKANTU_GET_MACRO(SpatialDimension, Model::spatial_dimension, UInt);

  /// return the damage array
  AKANTU_GET_MACRO_DEREF_PTR(Damage, damage);

  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(Damage, damage);
  
  /// get the PhaseFieldModel::internal_force vector (internal forces)
  AKANTU_GET_MACRO_DEREF_PTR(InternalForce, internal_force);

  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(InternalForce, internal_force);
  
  /// get the PhaseFieldModel::external_force vector (external forces)
  AKANTU_GET_MACRO_DEREF_PTR(ExternalForce, external_force);

  AKANTU_GET_MACRO_DEREF_PTR_NOT_CONST(ExternalForce, external_force);
  
  /// get the PhaseFieldModel::force vector (external forces)
  Array<Real> & getForce() {
    AKANTU_DEBUG_WARNING("getForce was maintained for backward compatibility, "
                         "use getExternalForce instead");
    return *external_force;
  }

  /// get the PhaseFieldModel::blocked_dofs vector
  AKANTU_GET_MACRO_DEREF_PTR(BlockedDOFs, blocked_dofs);

  
  /// get an iterable on the phasefields
  inline decltype(auto) getPhaseFields();

  /// get an iterable on the phasefields
  inline decltype(auto) getPhaseFields() const;

  /// get a particular phasefield (by phasefield index)
  inline PhaseField & getPhaseField(UInt mat_index);

  /// get a particular phasefield (by phasefield index)
  inline const PhaseField & getPhaseField(UInt mat_index) const;

  /// get a particular phasefield (by phasefield name)
  inline PhaseField & getPhaseField(const std::string & name);

  /// get a particular phasefield (by phasefield name)
  inline const PhaseField & getPhaseField(const std::string & name) const;

  /// get a particular phasefield id from is name
  inline UInt getPhaseFieldIndex(const std::string & name) const;

  /// give the number of phasefields
  inline UInt getNbPhaseFields() const { return phasefields.size(); }

  /// give the phasefield internal index from its id
  Int getInternalIndexFromID(const ID & id) const;

  AKANTU_GET_MACRO(PhaseFieldByElement, phasefield_index,
                   const ElementTypeMapArray<UInt> &);
  AKANTU_GET_MACRO(PhaseFieldLocalNumbering, phasefield_local_numbering,
                   const ElementTypeMapArray<UInt> &);

  /// vectors containing local material element index for each global element
  /// index
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(PhaseFieldByElement, phasefield_index,
                                         UInt);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(PhaseFieldByElement, phasefield_index, UInt);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(PhaseFieldLocalNumbering,
                                         phasefield_local_numbering, UInt);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(PhaseFieldLocalNumbering,
                                   phasefield_local_numbering, UInt);

  AKANTU_GET_MACRO_NOT_CONST(PhaseFieldSelector, *phasefield_selector,
                             PhaseFieldSelector &);

  AKANTU_SET_MACRO(PhaseFieldSelector, phasefield_selector,
                   std::shared_ptr<PhaseFieldSelector>);

  FEEngine & getFEEngineBoundary(const ID & name = "") override;


  /* ------------------------------------------------------------------------ */
  /* Dumpable Interface                                                       */
  /* ------------------------------------------------------------------------ */
public:
  std::shared_ptr<dumpers::Field>
  createNodalFieldReal(const std::string & field_name,
                       const std::string & group_name,
                       bool padding_flag) override;

  std::shared_ptr<dumpers::Field>
  createNodalFieldBool(const std::string & field_name,
                       const std::string & group_name,
                       bool padding_flag) override;

  std::shared_ptr<dumpers::Field>
  createElementalField(const std::string & field_name,
                       const std::string & group_name, bool padding_flag,
                       UInt spatial_dimension,
                       ElementKind kind) override;

  virtual void dump(const std::string & dumper_name) override;

  virtual void dump(const std::string & dumper_name, UInt step) override;

  virtual void dump(const std::string & dumper_name, Real time, UInt step) override;

  void dump() override;

  virtual void dump(UInt step) override;

  virtual void dump(Real time, UInt step) override;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// number of iterations
  UInt n_iter;

  /// damage array
  std::unique_ptr<Array<Real>> damage;

  /// damage array at the previous time step
  std::unique_ptr<Array<Real>> previous_damage;

  /// boundary vector 
  std::unique_ptr<Array<bool>> blocked_dofs;

  /// external force vector
  std::unique_ptr<Array<Real>> external_force;

  /// residuals array
  std::unique_ptr<Array<Real>> internal_force;

  /// Arrays containing the phasefield index for each element
  ElementTypeMapArray<UInt> phasefield_index;

  /// Arrays containing the position in the element filter of the phasefield
  /// (phasefield's local numbering)
  ElementTypeMapArray<UInt> phasefield_local_numbering;

  /// class defining of to choose a phasefield
  std::shared_ptr<PhaseFieldSelector> phasefield_selector;
  
  /// mapping between phasefield name and phasefield internal id
  std::map<std::string, UInt> phasefields_names_to_id;
  
  /// list of used phasefields
  std::vector<std::unique_ptr<PhaseField>> phasefields;

  /// tells if the phasefield are instantiated
  bool are_phasefields_instantiated{false};
};

} // namespace akantu


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "phasefield.hh"
#include "parser.hh"

#include "phase_field_model_inline_impl.cc"
/* -------------------------------------------------------------------------- */

#endif
