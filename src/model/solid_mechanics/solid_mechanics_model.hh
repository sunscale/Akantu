/**
 * @file   solid_mechanics_model.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Jul 27 2010
 * @date last modification: Tue Jan 19 2016
 *
 * @brief  Model of Solid Mechanics
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "boundary_condition.hh"
#include "data_accessor.hh"
#include "model.hh"
#include "non_local_manager.hh"
#include "solid_mechanics_model_event_handler.hh"
/* -------------------------------------------------------------------------- */
#include "integrator_gauss.hh"
#include "shape_lagrange.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SOLID_MECHANICS_MODEL_HH__
#define __AKANTU_SOLID_MECHANICS_MODEL_HH__

namespace akantu {
class Material;
class MaterialSelector;
class DumperIOHelper;
class NonLocalManager;
} // namespace akantu

/* -------------------------------------------------------------------------- */
namespace akantu {

struct SolidMechanicsModelOptions : public ModelOptions {
  explicit SolidMechanicsModelOptions(
      AnalysisMethod analysis_method = _explicit_lumped_mass);

  template <typename... pack>
  SolidMechanicsModelOptions(use_named_args_t, pack &&... _pack);

  AnalysisMethod analysis_method;
  //bool no_init_materials;
};

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
class SolidMechanicsModel
    : public Model,
      public DataAccessor<Element>,
      public DataAccessor<UInt>,
      public BoundaryCondition<SolidMechanicsModel>,
      public NonLocalManagerCallback,
      public MeshEventHandler,
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
  SolidMechanicsModel(Mesh & mesh, UInt spatial_dimension = _all_dimensions,
                      const ID & id = "solid_mechanics_model",
                      const MemoryID & memory_id = 0);

  virtual ~SolidMechanicsModel();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  template <typename P, typename T, typename... pack>
  void initFull(named_argument::param_t<P, T &&> && first, pack &&... _pack) {
    this->initFull(SolidMechanicsModelOptions{use_named_args, first, _pack...});
  }

  /// initialize completely the model
  void initFull(
      const ModelOptions & options = SolidMechanicsModelOptions()) override;

  /// initialize the fem object needed for boundary conditions
  void initFEEngineBoundary();

  /// register the tags associated with the parallel synchronizer
  // virtual void initParallel(MeshPartition * partition,
  //                           DataAccessor<Element> * data_accessor = NULL);

  /// allocate all vectors
  virtual void initArrays();

  /// allocate all vectors
  // void initArraysPreviousDisplacment();

  /// initialize all internal arrays for materials
  virtual void initMaterials();

  /// initialize the model
  void initModel() override;

  /// init PBC synchronizer
  // void initPBC();

  /// initialize a new solver and sets it as the default one to use
  void initNewSolver(const AnalysisMethod & method);

  /// function to print the containt of the class
  void printself(std::ostream & stream, int indent = 0) const override;

protected:
  /// allocate an array if needed
  template <typename T>
  void allocNodalField(Array<T> *& array, const ID & name);

  /* ------------------------------------------------------------------------ */
  /* PBC                                                                      */
  /* ------------------------------------------------------------------------ */
public:
  /// change the equation number for proper assembly when using PBC
  //  void changeEquationNumberforPBC(std::map <UInt, UInt> & pbc_pair);

  /// synchronize Residual for output
  // void synchronizeResidual();

protected:
  /// register PBC synchronizer
  //  void registerPBCSynchronizer();

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

  /// callback for the solver, this assembles the stiffness matrix
  void assembleJacobian() override;

  /// callback for the solver, this is called at beginning of solve
  void predictor() override;
  /// callback for the solver, this is called at end of solve
  void corrector() override;

  /// Callback for the model to instantiate the matricees when needed
  void initSolver(TimeStepSolverType time_step_solver_type,
                  NonLinearSolverType non_linear_solver_type) override;

protected:
  /* ------------------------------------------------------------------------ */
  TimeStepSolverType getDefaultSolverType() const override;
  /* ------------------------------------------------------------------------ */
  ModelSolverOptions
  getDefaultSolverOptions(const TimeStepSolverType & type) const override;

  /* ------------------------------------------------------------------------ */
  /* Explicit                                                                 */
  /* ------------------------------------------------------------------------ */
  // public:
  //   /// initialize the stuff for the explicit scheme
  //   void initExplicit(AnalysisMethod analysis_method =
  //   _explicit_lumped_mass);

public:
  bool isDefaultSolverExplicit() {
    return method == _explicit_lumped_mass ||
           method == _explicit_consistent_mass;
  }

  //   /// initialize the array needed by updateResidual (residual,
  //   current_position)
  //   void initializeUpdateResidualData();

protected:
  /// update the current position vector
  void updateCurrentPosition();

  //   /// assemble the residual for the explicit scheme
  //   virtual void updateResidual(bool need_initialize = true);

  //   /**
  //    * \brief compute the acceleration from the residual
  //    * this function is the explicit equivalent to solveDynamic in implicit
  //    * In the case of lumped mass just divide the residual by the mass
  //    * In the case of not lumped mass call
  //    solveDynamic<_acceleration_corrector>
  //    */
  //   void updateAcceleration();

  /// Update the increment of displacement
  // void updateIncrement();

  // /// Copy the actuel displacement into previous displacement
  // void updatePreviousDisplacement();

  //   /// Save stress and strain through EventManager
  //   void saveStressAndStrainBeforeDamage();
  //   /// Update energies through EventManager
  //   void updateEnergiesAfterDamage();

  //   /// Solve the system @f[ A x = \alpha b @f] with A a lumped matrix
  //   void solveLumped(Array<Real> & x, const Array<Real> & A,
  //                    const Array<Real> & b, const Array<bool> & blocked_dofs,
  //                    Real alpha);

  //   /// explicit integration predictor
  //   void explicitPred();

  //   /// explicit integration corrector
  //   void explicitCorr();

  // public:
  //   void solveStep();

  //   /*
  //   ------------------------------------------------------------------------
  //   */
  //   /* Implicit */
  //   /*
  //   ------------------------------------------------------------------------
  //   */
  // public:
  //   /// initialize the solver and the jacobian_matrix (called by
  //   initImplicit)
  //   void initSolver();

  //   /// initialize the stuff for the implicit solver
  //   void initImplicit(bool            dynamic = false);

  //   /// solve Ma = f to get the initial acceleration
  //   void initialAcceleration();

  //   /// assemble the stiffness matrix
  //   void assembleStiffnessMatrix();

  // public:
  //   /**
  //    * solve a step (predictor + convergence loop + corrector) using the
  //    * the given convergence method (see akantu::SolveConvergenceMethod)
  //    * and the given convergence criteria (see
  //    * akantu::SolveConvergenceCriteria)
  //    **/
  //   template <SolveConvergenceMethod method, SolveConvergenceCriteria
  //   criteria>
  //   bool solveStep(Real tolerance, UInt max_iteration = 100);

  //   template <SolveConvergenceMethod method, SolveConvergenceCriteria
  //   criteria>
  //   bool solveStep(Real tolerance, Real & error, UInt max_iteration = 100,
  //                  bool do_not_factorize = false);

  // public:
  //   /**
  //    * solve Ku = f using the the given convergence method (see
  //    * akantu::SolveConvergenceMethod) and the given convergence
  //    * criteria (see akantu::SolveConvergenceCriteria)
  //    **/
  //   template <SolveConvergenceMethod cmethod, SolveConvergenceCriteria
  //   criteria>
  //   bool solveStatic(Real tolerance, UInt max_iteration,
  //                    bool do_not_factorize = false);

  //   /// create and return the velocity damping matrix
  //   SparseMatrix & initVelocityDampingMatrix();

  //   /// implicit time integration predictor
  //   void implicitPred();

  //   /// implicit time integration corrector
  //   void implicitCorr();

  //   /// compute the Cauchy stress on user demand.
  //   void computeCauchyStresses();

  //   // /// compute A and solve @f[ A\delta u = f_ext - f_int @f]
  //   // template <NewmarkBeta::IntegrationSchemeCorrectorType type>
  //   // void solve(Array<Real> &increment, Real block_val = 1.,
  //   //            bool need_factorize = true, bool has_profile_changed =
  //   false);

  // protected:
  //   /// finish the computation of residual to solve in increment
  //   void updateResidualInternal();

  //   /// compute the support reaction and store it in force
  //   void updateSupportReaction();

  // private:
  //   /// re-initialize the J matrix (to use if the profile of K changed)
  //   void initJacobianMatrix();

  /* ------------------------------------------------------------------------ */
  // public:
  /// Update the stresses for the computation of the residual of the Stiffness
  /// matrix in the case of finite deformation
  // void computeStresses();

  /// synchronize the ghost element boundaries values
  // void synchronizeBoundaries();

  /* ------------------------------------------------------------------------ */
  /* Materials (solid_mechanics_model_material.cc)                            */
  /* ------------------------------------------------------------------------ */
public:
  /// registers all the custom materials of a given type present in the input
  /// file
  // template <typename M> void registerNewCustomMaterials(const ID & mat_type);

  /// register an empty material of a given type
  Material & registerNewMaterial(const ID & mat_name, const ID & mat_type,
                                 const ID & opt_param);
  // /// Use a UIntData in the mesh to specify the material to use per
  // element void setMaterialIDsFromIntData(const std::string & data_name);

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
  void
  assignMaterialToElements(const ElementTypeMapArray<UInt> * filter = NULL);

  /// reinitialize dof_synchronizer and solver (either in implicit or
  /// explicit) when cohesive elements are inserted
  // void reinitializeSolver();

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

  /// assemble the lumped mass only for a given list of nodes
  void assembleMassLumped(const Array<UInt> & node_list);

  /// assemble the lumped mass only for a given list of nodes
  void assembleMass(const Array<UInt> & node_list);

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
  inline UInt getNbData(const Array<Element> & elements,
                        const SynchronizationTag & tag) const override;

  inline void packData(CommunicationBuffer & buffer,
                       const Array<Element> & elements,
                       const SynchronizationTag & tag) const override;

  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<Element> & elements,
                         const SynchronizationTag & tag) override;

  inline UInt getNbData(const Array<UInt> & dofs,
                        const SynchronizationTag & tag) const override;

  inline void packData(CommunicationBuffer & buffer, const Array<UInt> & dofs,
                       const SynchronizationTag & tag) const override;

  inline void unpackData(CommunicationBuffer & buffer, const Array<UInt> & dofs,
                         const SynchronizationTag & tag) override;

protected:
  inline void splitElementByMaterial(const Array<Element> & elements,
                                     Array<Element> * elements_per_mat) const;

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
#ifndef SWIG
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
#endif

  dumper::Field * createNodalFieldReal(const std::string & field_name,
                                       const std::string & group_name,
                                       bool padding_flag) override;

  dumper::Field * createNodalFieldBool(const std::string & field_name,
                                       const std::string & group_name,
                                       bool padding_flag) override;

  dumper::Field * createElementalField(const std::string & field_name,
                                       const std::string & group_name,
                                       bool padding_flag,
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

  /// get the current value of the time step
  // AKANTU_GET_MACRO(TimeStep, time_step, Real);

  /// set the value of the time step
  void setTimeStep(Real time_step, const ID & solver_id = "") override;
  /// void setTimeStep(Real time_step);

  /// return the of iterations done in the last solveStep
  // AKANTU_GET_MACRO(NumberIter, n_iter, UInt);

  /// get the value of the conversion from forces/ mass to acceleration
  AKANTU_GET_MACRO(F_M2A, f_m2a, Real);

  /// set the value of the conversion from forces/ mass to acceleration
  AKANTU_SET_MACRO(F_M2A, f_m2a, Real);

  /// get the SolidMechanicsModel::displacement vector
  AKANTU_GET_MACRO(Displacement, *displacement, Array<Real> &);

  /// get the SolidMechanicsModel::previous_displacement vector
  AKANTU_GET_MACRO(PreviousDisplacement, *previous_displacement, Array<Real> &);

  /// get the SolidMechanicsModel::current_position vector \warn only consistent
  /// after a call to SolidMechanicsModel::updateCurrentPosition
  const Array<Real> & getCurrentPosition();

  /// get  the SolidMechanicsModel::increment  vector \warn  only  consistent if
  /// SolidMechanicsModel::setIncrementFlagOn has been called before
  AKANTU_GET_MACRO(Increment, *displacement_increment, Array<Real> &);

  /// get the lumped SolidMechanicsModel::mass vector
  AKANTU_GET_MACRO(Mass, *mass, Array<Real> &);

  /// get the SolidMechanicsModel::velocity vector
  AKANTU_GET_MACRO(Velocity, *velocity, Array<Real> &);

  /// get    the    SolidMechanicsModel::acceleration    vector,   updated    by
  /// SolidMechanicsModel::updateAcceleration
  AKANTU_GET_MACRO(Acceleration, *acceleration, Array<Real> &);

  /// get the SolidMechanicsModel::force vector (external forces)
  AKANTU_GET_MACRO(Force, *external_force, Array<Real> &);

  /// get the SolidMechanicsModel::internal_force vector (internal forces)
  AKANTU_GET_MACRO(InternalForce, *internal_force, Array<Real> &);

  /// get the SolidMechanicsModel::blocked_dofs vector
  AKANTU_GET_MACRO(BlockedDOFs, *blocked_dofs, Array<bool> &);

  /// get the SolidMechnicsModel::incrementAcceleration vector
  // AKANTU_GET_MACRO(IncrementAcceleration, *increment_acceleration,
  //                  Array<Real> &);

  /// get the value of the SolidMechanicsModel::increment_flag
  AKANTU_GET_MACRO(IncrementFlag, increment_flag, bool);

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

  inline void setMaterialSelector(MaterialSelector & selector);

  /// give the material internal index from its id
  Int getInternalIndexFromID(const ID & id) const;

  /// compute the stable time step
  Real getStableTimeStep();

  /// get the energies
  Real getEnergy(const std::string & energy_id);

  /// compute the energy for energy
  Real getEnergy(const std::string & energy_id, const ElementType & type,
                 UInt index);

  /**
   * @brief set the SolidMechanicsModel::increment_flag  to on, the activate the
   * update of the SolidMechanicsModel::increment vector
   */
  // void setIncrementFlagOn();

  // /// get the stiffness matrix
  // AKANTU_GET_MACRO(StiffnessMatrix, *stiffness_matrix, SparseMatrix &);

  // /// get the global jacobian matrix of the system
  // AKANTU_GET_MACRO(GlobalJacobianMatrix, *jacobian_matrix, const SparseMatrix
  // &);

  // /// get the mass matrix
  // AKANTU_GET_MACRO(MassMatrix, *mass_matrix, SparseMatrix &);

  // /// get the velocity damping matrix
  // AKANTU_GET_MACRO(VelocityDampingMatrix, *velocity_damping_matrix,
  // SparseMatrix &);

  /// get the FEEngine object to integrate or interpolate on the boundary
  inline FEEngine & getFEEngineBoundary(const ID & name = "");

  // /// get integrator
  // AKANTU_GET_MACRO(Integrator, *integrator, const IntegrationScheme2ndOrder
  // &);

  // /// get synchronizer
  // AKANTU_GET_MACRO(Synchronizer, *synch_parallel, const ElementSynchronizer
  // &);

  AKANTU_GET_MACRO(MaterialByElement, material_index,
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

  /// Get the type of analysis method used
  AKANTU_GET_MACRO(AnalysisMethod, method, AnalysisMethod);

  /// Access the non_local_manager interface
  AKANTU_GET_MACRO(NonLocalManager, *non_local_manager, NonLocalManager &);
protected:
  friend class Material;

protected:
  /// compute the stable time step
  Real getStableTimeStep(const GhostType & ghost_type);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// number of iterations
  // UInt n_iter;

  /// time step
  //  Real time_step;

  /// conversion coefficient form force/mass to acceleration
  Real f_m2a;

  /// displacements array
  Array<Real> * displacement;
  UInt displacement_release{0};

  /// displacements array at the previous time step (used in finite deformation)
  Array<Real> * previous_displacement;

  /// increment of displacement
  Array<Real> * displacement_increment;

  /// lumped mass array
  Array<Real> * mass;

  /// velocities array
  Array<Real> * velocity;

  /// accelerations array
  Array<Real> * acceleration;

  /// accelerations array
  // Array<Real> * increment_acceleration;

  /// external forces array
  Array<Real> * external_force;

  /// internal forces array
  Array<Real> * internal_force;

  /// array specifing if a degree of freedom is blocked or not
  Array<bool> * blocked_dofs;

  /// array of current position used during update residual
  Array<Real> * current_position;
  UInt current_position_release{0};

  /// mass matrix
  SparseMatrix * mass_matrix;

  /// velocity damping matrix
  SparseMatrix * velocity_damping_matrix;

  /// stiffness matrix
  SparseMatrix * stiffness_matrix;

  /// jacobian matrix @f[A = cM + dD + K@f] with @f[c = \frac{1}{\beta \Delta
  /// t^2}, d = \frac{\gamma}{\beta \Delta t} @f]
  SparseMatrix * jacobian_matrix;

  /// Arrays containing the material index for each element
  ElementTypeMapArray<UInt> material_index;

  /// Arrays containing the position in the element filter of the material
  /// (material's local numbering)
  ElementTypeMapArray<UInt> material_local_numbering;

#ifdef SWIGPYTHON
public:
#endif
  /// list of used materials
  std::vector<std::unique_ptr<Material>> materials;

  /// mapping between material name and material internal id
  std::map<std::string, UInt> materials_names_to_id;
#ifdef SWIGPYTHON
protected:
#endif

  /// class defining of to choose a material
  MaterialSelector * material_selector;

  /// define if it is the default selector or not
  bool is_default_material_selector;

  // /// integration scheme of second order used
  // IntegrationScheme2ndOrder *integrator;

  /// flag defining if the increment must be computed or not
  bool increment_flag;

  /// analysis method check the list in akantu::AnalysisMethod
  AnalysisMethod method;

  /// tells if the material are instantiated
  bool are_materials_instantiated;

  typedef std::map<std::pair<std::string, ElementKind>,
                   ElementTypeMapArray<Real> *>
      flatten_internal_map;

  /// map a registered internals to be flattened for dump purposes
  flatten_internal_map registered_internals;

  /// non local manager
  std::unique_ptr<NonLocalManager> non_local_manager;

  /// pointer to the pbc synchronizer
  // PBCSynchronizer * pbc_synch;
};

/* -------------------------------------------------------------------------- */
namespace BC {
  namespace Neumann {
    typedef FromHigherDim FromStress;
    typedef FromSameDim FromTraction;
  } // namespace Neumann
} // namespace BC

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream,
                                 const SolidMechanicsModel & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "material.hh"
#include "parser.hh"

#include "solid_mechanics_model_inline_impl.cc"
#include "solid_mechanics_model_tmpl.hh"

#include "material_selector_tmpl.hh"

#endif /* __AKANTU_SOLID_MECHANICS_MODEL_HH__ */
