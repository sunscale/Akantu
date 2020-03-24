/**
 * @file   solid_mechanics_model.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Jul 27 2010
 * @date last modification: Wed Feb 21 2018
 *
 * @brief  Model of Solid Mechanics
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
#include "non_local_manager_callback.hh"
#include "solid_mechanics_model_event_handler.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SOLID_MECHANICS_MODEL_HH__
#define __AKANTU_SOLID_MECHANICS_MODEL_HH__

namespace akantu {
class Material;
class MaterialSelector;
class DumperIOHelper;
class NonLocalManager;
template <ElementKind kind, class IntegrationOrderFunctor>
class IntegratorGauss;
template <ElementKind kind> class ShapeLagrange;
} // namespace akantu

/* -------------------------------------------------------------------------- */
namespace akantu {

/* -------------------------------------------------------------------------- */
class SolidMechanicsModel
    : public Model,
      public DataAccessor<Element>,
      public DataAccessor<UInt>,
      public BoundaryCondition<SolidMechanicsModel>,
      public NonLocalManagerCallback,
      public EventHandlerManager<SolidMechanicsModelEventHandler> {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  class NewMaterialElementsEvent : public NewElementsEvent {
  public:
    AKANTU_GET_MACRO_NOT_CONST(MaterialList, material, Array<UInt> &);
    AKANTU_GET_MACRO(MaterialList, material, const Array<UInt> &);

  protected:
    Array<UInt> material;
  };

  using MyFEEngineType = FEEngineTemplate<IntegratorGauss, ShapeLagrange>;

protected:
  using EventManager = EventHandlerManager<SolidMechanicsModelEventHandler>;

public:
  SolidMechanicsModel(
      Mesh & mesh, UInt spatial_dimension = _all_dimensions,
      const ID & id = "solid_mechanics_model", const MemoryID & memory_id = 0,
      const ModelType model_type = ModelType::_solid_mechanics_model);

  ~SolidMechanicsModel() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  /// initialize completely the model
  void initFullImpl(
      const ModelOptions & options = SolidMechanicsModelOptions()) override;

public:
  /// initialize all internal arrays for materials
  virtual void initMaterials();

protected:
  /// initialize the model
  void initModel() override;

  /// function to print the containt of the class
  void printself(std::ostream & stream, int indent = 0) const override;

  /// get some default values for derived classes
  std::tuple<ID, TimeStepSolverType>
  getDefaultSolverID(const AnalysisMethod & method) override;

  /* ------------------------------------------------------------------------ */
  /* Solver interface                                                         */
  /* ------------------------------------------------------------------------ */
public:
  /// assembles the stiffness matrix,
  virtual void assembleStiffnessMatrix();
  /// assembles the internal forces in the array internal_forces
  virtual void assembleInternalForces();

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

  /// callback for the solver, this is called at beginning of solve
  void predictor() override;
  /// callback for the solver, this is called at end of solve
  void corrector() override;

  /// callback for the solver, this is called at beginning of solve
  void beforeSolveStep() override;
  /// callback for the solver, this is called at end of solve
  void afterSolveStep(bool converted = true) override;

  /// Callback for the model to instantiate the matricees when needed
  void initSolver(TimeStepSolverType time_step_solver_type,
                  NonLinearSolverType non_linear_solver_type) override;

protected:
  /* ------------------------------------------------------------------------ */
  TimeStepSolverType getDefaultSolverType() const override;
  /* ------------------------------------------------------------------------ */
  ModelSolverOptions
  getDefaultSolverOptions(const TimeStepSolverType & type) const override;

public:
  bool isDefaultSolverExplicit() {
    return method == _explicit_lumped_mass ||
           method == _explicit_consistent_mass;
  }

protected:
  /// update the current position vector
  void updateCurrentPosition();

  /* ------------------------------------------------------------------------ */
  /* Materials (solid_mechanics_model_material.cc)                            */
  /* ------------------------------------------------------------------------ */
public:
  /// register an empty material of a given type
  Material & registerNewMaterial(const ID & mat_name, const ID & mat_type,
                                 const ID & opt_param);

  /// reassigns materials depending on the material selector
  virtual void reassignMaterial();

  /// apply a constant eigen_grad_u on all quadrature points of a given material
  virtual void applyEigenGradU(const Matrix<Real> & prescribed_eigen_grad_u,
                               const ID & material_name,
                               const GhostType ghost_type = _not_ghost);

protected:
  /// register a material in the dynamic database
  Material & registerNewMaterial(const ParserSection & mat_section);

  /// read the material files to instantiate all the materials
  void instantiateMaterials();

  /// set the element_id_by_material and add the elements to the good materials
  virtual void
  assignMaterialToElements(const ElementTypeMapArray<UInt> * filter = nullptr);

  /* ------------------------------------------------------------------------ */
  /* Mass (solid_mechanics_model_mass.cc)                                     */
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

  /// fill a vector of rho
  void computeRho(Array<Real> & rho, ElementType type, GhostType ghost_type);

  /// compute the kinetic energy
  Real getKineticEnergy();
  Real getKineticEnergy(const ElementType & type, UInt index);

  /// compute the external work (for impose displacement, the velocity should be
  /// given too)
  Real getExternalWork();

  /* ------------------------------------------------------------------------ */
  /* NonLocalManager inherited members                                        */
  /* ------------------------------------------------------------------------ */
protected:
  void initializeNonLocal() override;

  void updateDataForNonLocalCriterion(ElementTypeMapReal & criterion) override;

  void computeNonLocalStresses(const GhostType & ghost_type) override;

  void
  insertIntegrationPointsInNeighborhoods(const GhostType & ghost_type) override;

  /// update the values of the non local internal
  void updateLocalInternal(ElementTypeMapReal & internal_flat,
                           const GhostType & ghost_type,
                           const ElementKind & kind) override;

  /// copy the results of the averaging in the materials
  void updateNonLocalInternal(ElementTypeMapReal & internal_flat,
                              const GhostType & ghost_type,
                              const ElementKind & kind) override;

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

  UInt getNbData(const Array<UInt> & dofs,
                 const SynchronizationTag & tag) const override;

  void packData(CommunicationBuffer & buffer, const Array<UInt> & dofs,
                const SynchronizationTag & tag) const override;

  void unpackData(CommunicationBuffer & buffer, const Array<UInt> & dofs,
                  const SynchronizationTag & tag) override;

protected:
  void
  splitElementByMaterial(const Array<Element> & elements,
                         std::vector<Array<Element>> & elements_per_mat) const;

  template <typename Operation>
  void splitByMaterial(const Array<Element> & elements, Operation && op) const;

  /* ------------------------------------------------------------------------ */
  /* Mesh Event Handler inherited members                                     */
  /* ------------------------------------------------------------------------ */
protected:
  void onNodesAdded(const Array<UInt> & nodes_list,
                    const NewNodesEvent & event) override;
  void onNodesRemoved(const Array<UInt> & element_list,
                      const Array<UInt> & new_numbering,
                      const RemovedNodesEvent & event) override;
  void onElementsAdded(const Array<Element> & nodes_list,
                       const NewElementsEvent & event) override;
  void onElementsRemoved(const Array<Element> & element_list,
                         const ElementTypeMapArray<UInt> & new_numbering,
                         const RemovedElementsEvent & event) override;
  void onElementsChanged(const Array<Element> &, const Array<Element> &,
                         const ElementTypeMapArray<UInt> &,
                         const ChangedElementsEvent &) override{};

  /* ------------------------------------------------------------------------ */
  /* Dumpable interface (kept for convenience) and dumper relative functions  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void onDump();

  //! decide wether a field is a material internal or not
  bool isInternal(const std::string & field_name,
                  const ElementKind & element_kind);
  //! give the amount of data per element
  virtual ElementTypeMap<UInt>
  getInternalDataPerElem(const std::string & field_name,
                         const ElementKind & kind);

  //! flatten a given material internal field
  ElementTypeMapArray<Real> &
  flattenInternal(const std::string & field_name, const ElementKind & kind,
                  const GhostType ghost_type = _not_ghost);
  //! flatten all the registered material internals
  void flattenAllRegisteredInternals(const ElementKind & kind);

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
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// return the dimension of the system space
  AKANTU_GET_MACRO(SpatialDimension, Model::spatial_dimension, UInt);

  /// set the value of the time step
  void setTimeStep(Real time_step, const ID & solver_id = "") override;

  /// get the value of the conversion from forces/ mass to acceleration
  AKANTU_GET_MACRO(F_M2A, f_m2a, Real);

  /// set the value of the conversion from forces/ mass to acceleration
  AKANTU_SET_MACRO(F_M2A, f_m2a, Real);

  /// get the SolidMechanicsModel::displacement vector
  AKANTU_GET_MACRO_DEREF_PTR(Displacement, displacement);

  /// get the SolidMechanicsModel::previous_displacement vector
  AKANTU_GET_MACRO_DEREF_PTR(PreviousDisplacement, previous_displacement);

  /// get the SolidMechanicsModel::current_position vector \warn only consistent
  /// after a call to SolidMechanicsModel::updateCurrentPosition
  const Array<Real> & getCurrentPosition();

  /// get  the SolidMechanicsModel::increment  vector \warn  only  consistent if
  AKANTU_GET_MACRO_DEREF_PTR(Increment, displacement_increment);

  /// get the lumped SolidMechanicsModel::mass vector
  AKANTU_GET_MACRO_DEREF_PTR(Mass, mass);

  /// get the SolidMechanicsModel::velocity vector
  AKANTU_GET_MACRO_DEREF_PTR(Velocity, velocity);

  /// get    the    SolidMechanicsModel::acceleration    vector,   updated    by
  /// SolidMechanicsModel::updateAcceleration
  AKANTU_GET_MACRO_DEREF_PTR(Acceleration, acceleration);

  /// get the SolidMechanicsModel::external_force vector (external forces)
  AKANTU_GET_MACRO_DEREF_PTR(ExternalForce, external_force);

  /// get the SolidMechanicsModel::force vector (external forces)
  Array<Real> & getForce() {
    AKANTU_DEBUG_WARNING("getForce was maintained for backward compatibility, "
                         "use getExternalForce instead");
    return getExternalForce();
  }

  /// get the SolidMechanicsModel::internal_force vector (internal forces)
  AKANTU_GET_MACRO_DEREF_PTR(InternalForce, internal_force);

  /// get the SolidMechanicsModel::blocked_dofs vector
  AKANTU_GET_MACRO_DEREF_PTR(BlockedDOFs, blocked_dofs);

  /// get an iterable on the materials
  inline decltype(auto) getMaterials();

  /// get an iterable on the materials
  inline decltype(auto) getMaterials() const;

  /// get a particular material (by material index)
  inline Material & getMaterial(UInt mat_index);

  /// get a particular material (by material index)
  inline const Material & getMaterial(UInt mat_index) const;

  /// get a particular material (by material name)
  inline Material & getMaterial(const std::string & name);

  /// get a particular material (by material name)
  inline const Material & getMaterial(const std::string & name) const;

  /// get a particular material id from is name
  inline UInt getMaterialIndex(const std::string & name) const;

  /// give the number of materials
  inline UInt getNbMaterials() const { return materials.size(); }

  /// give the material internal index from its id
  Int getInternalIndexFromID(const ID & id) const;

  /// compute the stable time step
  Real getStableTimeStep();

  /// get the energies
  Real getEnergy(const std::string & energy_id);

  /// compute the energy for energy
  Real getEnergy(const std::string & energy_id, const ElementType & type,
                 UInt index);

  AKANTU_GET_MACRO(MaterialByElement, material_index,
                   const ElementTypeMapArray<UInt> &);
  AKANTU_GET_MACRO(MaterialLocalNumbering, material_local_numbering,
                   const ElementTypeMapArray<UInt> &);

  /// vectors containing local material element index for each global element
  /// index
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(MaterialByElement, material_index,
                                         UInt);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(MaterialByElement, material_index, UInt);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(MaterialLocalNumbering,
                                         material_local_numbering, UInt);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(MaterialLocalNumbering,
                                   material_local_numbering, UInt);

  AKANTU_GET_MACRO_NOT_CONST(MaterialSelector, *material_selector,
                             MaterialSelector &);
  AKANTU_SET_MACRO(MaterialSelector, material_selector,
                   std::shared_ptr<MaterialSelector>);

  /// Access the non_local_manager interface
  AKANTU_GET_MACRO(NonLocalManager, *non_local_manager, NonLocalManager &);

  /// get the FEEngine object to integrate or interpolate on the boundary
  FEEngine & getFEEngineBoundary(const ID & name = "") override;

protected:
  /// compute the stable time step
  Real getStableTimeStep(const GhostType & ghost_type);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// release version of the displacement array
  UInt displacement_release{0};

  /// release version of the current_position array
  UInt current_position_release{0};

  /// Check if materials need to recompute the mass array
  bool need_to_reassemble_lumped_mass{true};
  /// Check if materials need to recompute the mass matrix
  bool need_to_reassemble_mass{true};

  /// mapping between material name and material internal id
  std::map<std::string, UInt> materials_names_to_id;

protected:
  /// conversion coefficient form force/mass to acceleration
  Real f_m2a{1.0};

  /// displacements array
  std::unique_ptr<Array<Real>> displacement;

  /// displacements array at the previous time step (used in finite deformation)
  std::unique_ptr<Array<Real>> previous_displacement;

  /// increment of displacement
  std::unique_ptr<Array<Real>> displacement_increment;

  /// lumped mass array
  std::unique_ptr<Array<Real>> mass;

  /// velocities array
  std::unique_ptr<Array<Real>> velocity;

  /// accelerations array
  std::unique_ptr<Array<Real>> acceleration;

  /// external forces array
  std::unique_ptr<Array<Real>> external_force;

  /// internal forces array
  std::unique_ptr<Array<Real>> internal_force;

  /// array specifing if a degree of freedom is blocked or not
  std::unique_ptr<Array<bool>> blocked_dofs;

  /// array of current position used during update residual
  std::unique_ptr<Array<Real>> current_position;

  /// Arrays containing the material index for each element
  ElementTypeMapArray<UInt> material_index;

  /// Arrays containing the position in the element filter of the material
  /// (material's local numbering)
  ElementTypeMapArray<UInt> material_local_numbering;

  /// list of used materials
  std::vector<std::unique_ptr<Material>> materials;

  /// class defining of to choose a material
  std::shared_ptr<MaterialSelector> material_selector;

  using flatten_internal_map =
      std::map<std::pair<std::string, ElementKind>,
               std::unique_ptr<ElementTypeMapArray<Real>>>;

  /// map a registered internals to be flattened for dump purposes
  flatten_internal_map registered_internals;

  /// non local manager
  std::unique_ptr<NonLocalManager> non_local_manager;

  /// tells if the material are instantiated
  bool are_materials_instantiated{false};
};

/* -------------------------------------------------------------------------- */
namespace BC {
  namespace Neumann {
    using FromStress = FromHigherDim;
    using FromTraction = FromSameDim;
  } // namespace Neumann
} // namespace BC

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "material.hh"
#include "parser.hh"

#include "solid_mechanics_model_inline_impl.cc"
#include "solid_mechanics_model_tmpl.hh"
/* -------------------------------------------------------------------------- */

#endif /* __AKANTU_SOLID_MECHANICS_MODEL_HH__ */
