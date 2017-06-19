/**
 * @file   heat_transfer_model.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Srinivasa Babu Ramisetti <srinivasa.ramisetti@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Rui Wang <rui.wang@epfl.ch>
 *
 * @date creation: Sun May 01 2011
 * @date last modification: Tue Dec 08 2015
 *
 * @brief  Model of Heat Transfer
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

#ifndef __AKANTU_HEAT_TRANSFER_MODEL_HH__
#define __AKANTU_HEAT_TRANSFER_MODEL_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_types.hh"
#include "aka_voigthelper.hh"
#include "aka_memory.hh"
#include "model.hh"
#include "integrator_gauss.hh"
#include "shape_lagrange.hh"
#include "dumpable.hh"
#include "parsable.hh"
#include "solver.hh"
#include "generalized_trapezoidal.hh"

namespace akantu {
  class IntegrationScheme1stOrder;
}

namespace akantu {

struct HeatTransferModelOptions : public ModelOptions {
  HeatTransferModelOptions(AnalysisMethod analysis_method = _explicit_lumped_capacity ) : analysis_method(analysis_method) {}
  AnalysisMethod analysis_method;
};

extern const HeatTransferModelOptions default_heat_transfer_model_options;

class HeatTransferModel : public Model, public DataAccessor, public Parsable {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  typedef FEEngineTemplate<IntegratorGauss,ShapeLagrange> MyFEEngineType;

  HeatTransferModel(Mesh & mesh,
		    UInt spatial_dimension = _all_dimensions,
		    const ID & id = "heat_transfer_model",
		    const MemoryID & memory_id = 0);

  virtual ~HeatTransferModel() ;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

public:

  /// generic function to initialize everything ready for explicit dynamics
  void initFull(const ModelOptions & options = default_heat_transfer_model_options);

  /// initialize the fem object of the boundary
  void initFEEngineBoundary(bool create_surface = true);

  /// read one material file to instantiate all the materials
  void readMaterials();

  /// allocate all vectors
  void initArrays();

  /// register the tags associated with the parallel synchronizer
  void initParallel(MeshPartition * partition, DataAccessor * data_accessor=NULL);

  /// initialize the model
  void initModel();

  /// init PBC synchronizer
  void initPBC();

  /// initialize the solver and the jacobian_matrix (called by initImplicit)
  void initSolver(SolverOptions & solver_options);

  /// initialize the stuff for the implicit solver
  void initImplicit(bool dynamic, SolverOptions & solver_options = _solver_no_options);

  /// function to print the contain of the class
  virtual void printself(__attribute__ ((unused)) std::ostream & stream,
			 __attribute__ ((unused)) int indent = 0) const {};

  /* ------------------------------------------------------------------------ */
  /* Methods for explicit                                                     */
  /* ------------------------------------------------------------------------ */
public:

  /// compute and get the stable time step
  Real getStableTimeStep();

  /// compute the heat flux
  void updateResidual(bool compute_conductivity = false);

  /// calculate the lumped capacity vector for heat transfer problem
  void assembleCapacityLumped();

  /// update the temperature from the temperature rate
  void explicitPred();

  /// update the temperature rate from the increment
  void explicitCorr();

  /// implicit time integration predictor
  void implicitPred();

  /// implicit time integration corrector
  void implicitCorr();

  
  /// solve the system in temperature rate  @f$C\delta \dot T = q_{n+1} - C \dot T_{n}@f$
  /// this function needs to be run for dynamics
  void solveExplicitLumped();

  // /// initialize the heat flux
  // void initializeResidual(Array<Real> &temp);
  // /// initialize temperature
  // void initializeTemperature(Array<Real> &temp);

  /* ------------------------------------------------------------------------ */
  /* Methods for implicit                                                     */
  /* ------------------------------------------------------------------------ */
public:

  /**
   * solve Kt = q
   **/
  void solveStatic();

  /// test if the system is converged
  template <SolveConvergenceCriteria criteria>
  bool testConvergence(Real tolerance, Real & error);

  
    /**
   * solve a step (predictor + convergence loop + corrector) using the
   * the given convergence method (see akantu::SolveConvergenceMethod)
   * and the given convergence criteria (see
   * akantu::SolveConvergenceCriteria)
   **/
  template <SolveConvergenceMethod method, SolveConvergenceCriteria criteria>
  bool solveStep(Real tolerance, UInt max_iteration = 100);

  template <SolveConvergenceMethod method, SolveConvergenceCriteria criteria>
  bool solveStep(Real   tolerance,
		 Real & error,
		 UInt   max_iteration = 100,
		 bool   do_not_factorize = false);

  /// assemble the conductivity matrix
  void assembleConductivityMatrix(bool compute_conductivity = true);

  /// assemble the conductivity matrix
  void assembleCapacity();

  /// assemble the conductivity matrix
  void assembleCapacity(GhostType ghost_type);

  /// compute the capacity on quadrature points
  void computeRho(Array<Real> & rho, ElementType type,
		  GhostType ghost_type);


protected:
  /// compute A and solve @f[ A\delta u = f_ext - f_int @f]
  template <GeneralizedTrapezoidal::IntegrationSchemeCorrectorType type>
  void solve(Array<Real> &increment, Real block_val = 1.,
             bool need_factorize = true, bool has_profile_changed = false);

  /// computation of the residual for the convergence loop
  void updateResidualInternal();
  
private:

  /// compute the heat flux on ghost types
  void updateResidual(const GhostType & ghost_type, bool compute_conductivity = false);

  /// calculate the lumped capacity vector for heat transfer problem (w ghosttype)
  void assembleCapacityLumped(const GhostType & ghost_type);

  /// assemble the conductivity matrix (w/ ghost type)
  template <UInt dim>
  void assembleConductivityMatrix(const GhostType & ghost_type,
				  bool compute_conductivity = true);

  /// assemble the conductivity matrix
  template <UInt dim>
  void assembleConductivityMatrix(const ElementType & type,
				  const GhostType & ghost_type,
				  bool compute_conductivity = true);

  /// compute the conductivity tensor for each quadrature point in an array
  void computeConductivityOnQuadPoints(const GhostType & ghost_type);

  /// compute vector k \grad T for each quadrature point
  void computeKgradT(const GhostType & ghost_type,bool compute_conductivity);

  /// compute the thermal energy
  Real computeThermalEnergyByNode();

  /* ------------------------------------------------------------------------ */
  /* Data Accessor inherited members                                          */
  /* ------------------------------------------------------------------------ */
public:
  inline UInt getNbDataForElements(const Array<Element> & elements,
				   SynchronizationTag tag) const;
  inline void packElementData(CommunicationBuffer & buffer,
			      const Array<Element> & elements,
			      SynchronizationTag tag) const;
  inline void unpackElementData(CommunicationBuffer & buffer,
				const Array<Element> & elements,
				SynchronizationTag tag);

  inline UInt getNbDataToPack(SynchronizationTag tag) const;
  inline UInt getNbDataToUnpack(SynchronizationTag tag) const;
  inline void packData(CommunicationBuffer & buffer,
		       const UInt index,
		       SynchronizationTag tag) const;
  inline void unpackData(CommunicationBuffer & buffer,
			 const UInt index,
			 SynchronizationTag tag);

  /* ------------------------------------------------------------------------ */
  /* Dumpable interface                                                       */
  /* ------------------------------------------------------------------------ */
public:

  virtual dumper::Field * createNodalFieldReal(const std::string & field_name,
					       const std::string & group_name,
					       bool padding_flag);

  virtual dumper::Field * createNodalFieldBool(const std::string & field_name,
					       const std::string & group_name,
					       bool padding_flag);


  virtual dumper::Field * createElementalField(const std::string & field_name, 
					       const std::string & group_name,
					       bool padding_flag,
					       const UInt & spatial_dimension,
					       const ElementKind & kind);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  inline FEEngine & getFEEngineBoundary(const std::string & name = "");

  AKANTU_GET_MACRO(Density, density, Real);
  AKANTU_GET_MACRO(Capacity, capacity, Real);
  /// get the dimension of the system space
  AKANTU_GET_MACRO(SpatialDimension, spatial_dimension, UInt);
  /// get the current value of the time step
  AKANTU_GET_MACRO(TimeStep, time_step, Real);
  /// set the value of the time step
  AKANTU_SET_MACRO(TimeStep, time_step, Real);
  /// get the assembled heat flux
  AKANTU_GET_MACRO(Residual, *residual, Array<Real>&);
  /// get the lumped capacity
  AKANTU_GET_MACRO(CapacityLumped, *capacity_lumped, Array<Real>&);
  /// get the boundary vector
  AKANTU_GET_MACRO(BlockedDOFs, *blocked_dofs, Array<bool>&);
  /// get conductivity matrix
  AKANTU_GET_MACRO(ConductivityMatrix, *conductivity_matrix, const SparseMatrix&);
  /// get the external heat rate vector
  AKANTU_GET_MACRO(ExternalHeatRate, *external_heat_rate, Array<Real>&);
  /// get the temperature gradient
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(TemperatureGradient, temperature_gradient, Real);
  /// get the conductivity on q points
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(ConductivityOnQpoints, conductivity_on_qpoints, Real);
  /// get the conductivity on q points
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(TemperatureOnQpoints, temperature_on_qpoints, Real);
  /// internal variables
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(KGradtOnQpoints,  k_gradt_on_qpoints, Real);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(IntBtKgT, int_bt_k_gT, Real);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(BtKgT, bt_k_gT, Real);
  /// get the temperature
  AKANTU_GET_MACRO(Temperature, *temperature, Array<Real> &);
  /// get the temperature derivative
  AKANTU_GET_MACRO(TemperatureRate, *temperature_rate, Array<Real> &);
  /// get the equation number Array<Int>
  AKANTU_GET_MACRO(EquationNumber, *equation_number, const Array<Int> &);

  /// get the energy denominated by thermal
  Real getEnergy(const std::string & energy_id, const ElementType & type, UInt index);
  /// get the energy denominated by thermal
  Real getEnergy(const std::string & energy_id);

  /// get the thermal energy for a given element
  Real getThermalEnergy(const ElementType & type, UInt index);
  /// get the thermal energy for a given element
  Real getThermalEnergy();

protected:
  /* ----------------------------------------------------------------------- */
  template<class iterator>
  void getThermalEnergy(iterator Eth,
			Array<Real>::const_iterator<Real> T_it,
			Array<Real>::const_iterator<Real> T_end) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

    /// number of iterations
  UInt n_iter;

  IntegrationScheme1stOrder * integrator;

  /// time step
  Real time_step;

  /// temperatures array
  Array<Real> * temperature;

  /// temperatures derivatives array
  Array<Real> * temperature_rate;

  /// increment array (@f$\delta \dot T@f$ or @f$\delta T@f$)
  Array<Real> * increment;

  /// conductivity matrix
  SparseMatrix * conductivity_matrix;

  /// capacity matrix
  SparseMatrix *capacity_matrix;

  /// jacobian matrix
  SparseMatrix * jacobian_matrix;

  /// the density
  Real density;

  /// the speed of the changing temperature
  ElementTypeMapArray<Real> temperature_gradient;

  /// temperature field on quadrature points
  ElementTypeMapArray<Real> temperature_on_qpoints;

  /// conductivity tensor on quadrature points
  ElementTypeMapArray<Real> conductivity_on_qpoints;

  /// vector k \grad T on quad points
  ElementTypeMapArray<Real> k_gradt_on_qpoints;

  /// vector \int \grad N k \grad T
  ElementTypeMapArray<Real> int_bt_k_gT;

  /// vector \grad N k \grad T
  ElementTypeMapArray<Real> bt_k_gT;

  /// external flux vector
  Array<Real> * external_heat_rate;

  /// residuals array
  Array<Real> * residual;

  /// position of a dof in the K matrix
  Array<Int> * equation_number;

  //lumped vector
  Array<Real> * capacity_lumped;

  /// boundary vector
  Array<bool> * blocked_dofs;

  //realtime
  Real time;

  ///capacity
  Real capacity;

  //conductivity matrix
  Matrix<Real> conductivity;

  //linear variation of the conductivity (for temperature dependent conductivity)
  Real conductivity_variation;

  // reference temperature for the interpretation of temperature variation
  Real T_ref;

  //the biggest parameter of conductivity matrix
  Real conductivitymax;

  /// thermal energy by element
  ElementTypeMapArray<Real> thermal_energy;

  /// Solver
  Solver * solver;

  /// analysis method
  AnalysisMethod method;

  /// pointer to the pbc synchronizer
  PBCSynchronizer * pbc_synch;

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "heat_transfer_model_inline_impl.cc"
#endif

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const HeatTransferModel & _this)
{
  _this.printself(stream);
  return stream;
}


} // akantu



#endif /* __AKANTU_HEAT_TRANSFER_MODEL_HH__ */
