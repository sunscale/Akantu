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
 * @date last modification: Mon Feb 05 2018
 *
 * @brief  Model of Heat Transfer
 *
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
#include "data_accessor.hh"
#include "fe_engine.hh"
#include "model.hh"
/* -------------------------------------------------------------------------- */
#include <array>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_HEAT_TRANSFER_MODEL_HH__
#define __AKANTU_HEAT_TRANSFER_MODEL_HH__

namespace akantu {
template <ElementKind kind, class IntegrationOrderFunctor>
class IntegratorGauss;
template <ElementKind kind> class ShapeLagrange;
} // namespace akantu

namespace akantu {

class HeatTransferModel : public Model,
                          public DataAccessor<Element>,
                          public DataAccessor<UInt> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using FEEngineType = FEEngineTemplate<IntegratorGauss, ShapeLagrange>;

  HeatTransferModel(Mesh & mesh, UInt spatial_dimension = _all_dimensions,
                    const ID & id = "heat_transfer_model",
                    const MemoryID & memory_id = 0);

  virtual ~HeatTransferModel();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  /// generic function to initialize everything ready for explicit dynamics
  void initFullImpl(const ModelOptions & options) override;

  /// read one material file to instantiate all the materials
  void readMaterials();

  /// allocate all vectors
  void initSolver(TimeStepSolverType, NonLinearSolverType) override;

  /// initialize the model
  void initModel() override;

  void predictor() override;

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
  /* ------------------------------------------------------------------------ */
  /* Methods for explicit                                                     */
  /* ------------------------------------------------------------------------ */
public:
  /// compute and get the stable time step
  Real getStableTimeStep();

  /// set the stable timestep
  void setTimeStep(Real time_step, const ID & solver_id = "") override;

  // temporary protection to prevent bad usage: should check for bug
protected:
  /// compute the internal heat flux \todo Need code review: currently not
  /// public method
  void assembleInternalHeatRate();

public:
  /// calculate the lumped capacity vector for heat transfer problem
  void assembleCapacityLumped();

public:
  /// assemble the conductivity matrix
  void assembleConductivityMatrix();

  /// assemble the conductivity matrix
  void assembleCapacity();

  /// compute the capacity on quadrature points
  void computeRho(Array<Real> & rho, ElementType type, GhostType ghost_type);

private:
  /// calculate the lumped capacity vector for heat transfer problem (w
  /// ghost type)
  void assembleCapacityLumped(const GhostType & ghost_type);

  /// compute the conductivity tensor for each quadrature point in an array
  void computeConductivityOnQuadPoints(const GhostType & ghost_type);

  /// compute vector \f[k \grad T\f] for each quadrature point
  void computeKgradT(const GhostType & ghost_type);

  /// compute the thermal energy
  Real computeThermalEnergyByNode();

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

  inline UInt getNbData(const Array<UInt> & indexes,
                        const SynchronizationTag & tag) const override;
  inline void packData(CommunicationBuffer & buffer,
                       const Array<UInt> & indexes,
                       const SynchronizationTag & tag) const override;
  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<UInt> & indexes,
                         const SynchronizationTag & tag) override;

  /* ------------------------------------------------------------------------ */
  /* Dumpable interface                                                       */
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
  AKANTU_GET_MACRO(Density, density, Real);
  AKANTU_GET_MACRO(Capacity, capacity, Real);
  /// get the dimension of the system space
  AKANTU_GET_MACRO(SpatialDimension, spatial_dimension, UInt);
  /// get the current value of the time step
  AKANTU_GET_MACRO(TimeStep, time_step, Real);
  /// get the assembled heat flux
  AKANTU_GET_MACRO(InternalHeatRate, *internal_heat_rate, Array<Real> &);
  /// get the boundary vector
  AKANTU_GET_MACRO(BlockedDOFs, *blocked_dofs, Array<bool> &);
  /// get the external heat rate vector
  AKANTU_GET_MACRO(ExternalHeatRate, *external_heat_rate, Array<Real> &);
  /// get the temperature gradient
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(TemperatureGradient,
                                         temperature_gradient, Real);
  /// get the conductivity on q points
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(ConductivityOnQpoints,
                                         conductivity_on_qpoints, Real);
  /// get the conductivity on q points
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(TemperatureOnQpoints,
                                         temperature_on_qpoints, Real);
  /// internal variables
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(KgradT, k_gradt_on_qpoints, Real);
  /// get the temperature
  AKANTU_GET_MACRO(Temperature, *temperature, Array<Real> &);
  /// get the temperature derivative
  AKANTU_GET_MACRO(TemperatureRate, *temperature_rate, Array<Real> &);

  /// get the energy denominated by thermal
  Real getEnergy(const std::string & energy_id, const ElementType & type,
                 UInt index);
  /// get the energy denominated by thermal
  Real getEnergy(const std::string & energy_id);

  /// get the thermal energy for a given element
  Real getThermalEnergy(const ElementType & type, UInt index);
  /// get the thermal energy for a given element
  Real getThermalEnergy();

protected:
  /* ------------------------------------------------------------------------ */
  FEEngine & getFEEngineBoundary(const ID & name = "") override;

  /* ----------------------------------------------------------------------- */
  template <class iterator>
  void getThermalEnergy(iterator Eth, Array<Real>::const_iterator<Real> T_it,
                        Array<Real>::const_iterator<Real> T_end) const;

  template <typename T>
  void allocNodalField(Array<T> *& array, const ID & name);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// number of iterations
  //UInt n_iter;

  /// time step
  Real time_step;

  /// temperatures array
  Array<Real> * temperature{nullptr};

  /// temperatures derivatives array
  Array<Real> * temperature_rate{nullptr};

  /// increment array (@f$\delta \dot T@f$ or @f$\delta T@f$)
  Array<Real> * increment{nullptr};

  /// the density
  Real density;

  /// the speed of the changing temperature
  ElementTypeMapArray<Real> temperature_gradient;

  /// temperature field on quadrature points
  ElementTypeMapArray<Real> temperature_on_qpoints;

  /// conductivity tensor on quadrature points
  ElementTypeMapArray<Real> conductivity_on_qpoints;

  /// vector \f[k \grad T\f] on quad points
  ElementTypeMapArray<Real> k_gradt_on_qpoints;

  /// external flux vector
  Array<Real> * external_heat_rate{nullptr};

  /// residuals array
  Array<Real> * internal_heat_rate{nullptr};

  /// boundary vector
  Array<bool> * blocked_dofs{nullptr};

  // realtime
  //Real time;

  /// capacity
  Real capacity;

  // conductivity matrix
  Matrix<Real> conductivity;

  // linear variation of the conductivity (for temperature dependent
  // conductivity)
  Real conductivity_variation;

  // reference temperature for the interpretation of temperature variation
  Real T_ref;

  // the biggest parameter of conductivity matrix
  //Real conductivitymax;

  bool need_to_reassemble_capacity{true};
  bool need_to_reassemble_capacity_lumped{true};
  UInt temperature_release{0};
  UInt conductivity_matrix_release{UInt(-1)};
  std::unordered_map<GhostType, bool> initial_conductivity{{_not_ghost, true},
                                                           {_ghost, true}};
  std::unordered_map<GhostType, UInt> conductivity_release{{_not_ghost, 0},
                                                           {_ghost, 0}};
};

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "heat_transfer_model_inline_impl.hh"

#endif /* __AKANTU_HEAT_TRANSFER_MODEL_HH__ */
