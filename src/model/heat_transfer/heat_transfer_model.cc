/**
 * @file   heat_transfer_model.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Srinivasa Babu Ramisetti <srinivasa.ramisetti@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Rui Wang <rui.wang@epfl.ch>
 *
 * @date creation: Sun May 01 2011
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Implementation of HeatTransferModel class
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
#include "heat_transfer_model.hh"
#include "dumpable_inline_impl.hh"
#include "element_synchronizer.hh"
#include "fe_engine_template.hh"
#include "generalized_trapezoidal.hh"
#include "group_manager_inline_impl.hh"
#include "integrator_gauss.hh"
#include "mesh.hh"
#include "parser.hh"
#include "shape_lagrange.hh"

#ifdef AKANTU_USE_IOHELPER
#include "dumper_element_partition.hh"
#include "dumper_elemental_field.hh"
#include "dumper_internal_material_field.hh"
#include "dumper_iohelper_paraview.hh"
#endif

/* -------------------------------------------------------------------------- */
namespace akantu {

namespace heat_transfer {
  namespace details {
    class ComputeRhoFunctor {
    public:
      ComputeRhoFunctor(const HeatTransferModel & model) : model(model){};

      void operator()(Matrix<Real> & rho, const Element & /*unused*/) {
        rho.set(model.getCapacity() * model.getDensity());
      }

    private:
      const HeatTransferModel & model;
    };
  } // namespace details
} // namespace heat_transfer

/* -------------------------------------------------------------------------- */
HeatTransferModel::HeatTransferModel(Mesh & mesh, UInt dim, const ID & id)
    : Model(mesh, ModelType::_heat_transfer_model, dim, id),
      temperature_gradient("temperature_gradient", id),
      temperature_on_qpoints("temperature_on_qpoints", id),
      conductivity_on_qpoints("conductivity_on_qpoints", id),
      k_gradt_on_qpoints("k_gradt_on_qpoints", id) {
  AKANTU_DEBUG_IN();

  conductivity = Matrix<Real>(this->spatial_dimension, this->spatial_dimension);

  this->initDOFManager();

  this->registerDataAccessor(*this);

  if (this->mesh.isDistributed()) {
    auto & synchronizer = this->mesh.getElementSynchronizer();
    this->registerSynchronizer(synchronizer,
                               SynchronizationTag::_htm_temperature);
    this->registerSynchronizer(synchronizer,
                               SynchronizationTag::_htm_gradient_temperature);
  }

  registerFEEngineObject<FEEngineType>(id + ":fem", mesh, spatial_dimension);

#ifdef AKANTU_USE_IOHELPER
  this->mesh.registerDumper<DumperParaview>("heat_transfer", id, true);
  this->mesh.addDumpMesh(mesh, spatial_dimension, _not_ghost, _ek_regular);
#endif

  this->registerParam("conductivity", conductivity, _pat_parsmod);
  this->registerParam("conductivity_variation", conductivity_variation, 0.,
                      _pat_parsmod);
  this->registerParam("temperature_reference", T_ref, 0., _pat_parsmod);
  this->registerParam("capacity", capacity, _pat_parsmod);
  this->registerParam("density", density, _pat_parsmod);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::initModel() {
  auto & fem = this->getFEEngine();
  fem.initShapeFunctions(_not_ghost);
  fem.initShapeFunctions(_ghost);

  temperature_on_qpoints.initialize(fem, _nb_component = 1);
  temperature_gradient.initialize(fem, _nb_component = spatial_dimension);
  conductivity_on_qpoints.initialize(fem, _nb_component = spatial_dimension *
                                                          spatial_dimension);
  k_gradt_on_qpoints.initialize(fem, _nb_component = spatial_dimension);
}

/* -------------------------------------------------------------------------- */
FEEngine & HeatTransferModel::getFEEngineBoundary(const ID & name) {
  return aka::as_type<FEEngine>(getFEEngineClassBoundary<FEEngineType>(name));
}

/* -------------------------------------------------------------------------- */
HeatTransferModel::~HeatTransferModel() = default;

/* -------------------------------------------------------------------------- */
void HeatTransferModel::assembleCapacityLumped(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto & fem = getFEEngineClass<FEEngineType>();
  heat_transfer::details::ComputeRhoFunctor compute_rho(*this);

  for (auto && type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_regular)) {
    fem.assembleFieldLumped(compute_rho, "M", "temperature",
                            this->getDOFManager(), type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
MatrixType HeatTransferModel::getMatrixType(const ID & matrix_id) {
  if (matrix_id == "K" or matrix_id == "M") {
    return _symmetric;
  }

  return _mt_not_defined;
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::assembleMatrix(const ID & matrix_id) {
  if (matrix_id == "K") {
    this->assembleConductivityMatrix();
  } else if (matrix_id == "M" and need_to_reassemble_capacity) {
    this->assembleCapacity();
  }
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::assembleLumpedMatrix(const ID & matrix_id) {
  if (matrix_id == "M" and need_to_reassemble_capacity) {
    this->assembleCapacityLumped();
  }
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::assembleResidual() {
  AKANTU_DEBUG_IN();

  this->assembleInternalHeatRate();

  this->getDOFManager().assembleToResidual("temperature",
                                           *this->external_heat_rate, 1);
  this->getDOFManager().assembleToResidual("temperature",
                                           *this->internal_heat_rate, 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::predictor() { ++temperature_release; }

/* -------------------------------------------------------------------------- */
void HeatTransferModel::assembleCapacityLumped() {
  AKANTU_DEBUG_IN();

  if (!this->getDOFManager().hasLumpedMatrix("M")) {
    this->getDOFManager().getNewLumpedMatrix("M");
  }

  this->getDOFManager().zeroLumpedMatrix("M");

  assembleCapacityLumped(_not_ghost);
  assembleCapacityLumped(_ghost);

  need_to_reassemble_capacity_lumped = false;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::initSolver(TimeStepSolverType time_step_solver_type,
                                   NonLinearSolverType /*unused*/) {
  DOFManager & dof_manager = this->getDOFManager();

  this->allocNodalField(this->temperature, 1, "temperature");
  this->allocNodalField(this->external_heat_rate, 1, "external_heat_rate");
  this->allocNodalField(this->internal_heat_rate, 1, "internal_heat_rate");
  this->allocNodalField(this->blocked_dofs, 1, "blocked_dofs");

  if (!dof_manager.hasDOFs("temperature")) {
    dof_manager.registerDOFs("temperature", *this->temperature, _dst_nodal);
    dof_manager.registerBlockedDOFs("temperature", *this->blocked_dofs);
  }

  if (time_step_solver_type == TimeStepSolverType::_dynamic ||
      time_step_solver_type == TimeStepSolverType::_dynamic_lumped) {
    this->allocNodalField(this->temperature_rate, 1, "temperature_rate");

    if (!dof_manager.hasDOFsDerivatives("temperature", 1)) {
      dof_manager.registerDOFsDerivative("temperature", 1,
                                         *this->temperature_rate);
    }
  }
}

/* -------------------------------------------------------------------------- */
std::tuple<ID, TimeStepSolverType>
HeatTransferModel::getDefaultSolverID(const AnalysisMethod & method) {
  switch (method) {
  case _explicit_lumped_mass: {
    return std::make_tuple("explicit_lumped",
                           TimeStepSolverType::_dynamic_lumped);
  }
  case _static: {
    return std::make_tuple("static", TimeStepSolverType::_static);
  }
  case _implicit_dynamic: {
    return std::make_tuple("implicit", TimeStepSolverType::_dynamic);
  }
  default:
    return std::make_tuple("unknown", TimeStepSolverType::_not_defined);
  }
}

/* -------------------------------------------------------------------------- */
ModelSolverOptions HeatTransferModel::getDefaultSolverOptions(
    const TimeStepSolverType & type) const {
  ModelSolverOptions options;

  switch (type) {
  case TimeStepSolverType::_dynamic_lumped: {
    options.non_linear_solver_type = NonLinearSolverType::_lumped;
    options.integration_scheme_type["temperature"] =
        IntegrationSchemeType::_forward_euler;
    options.solution_type["temperature"] = IntegrationScheme::_temperature_rate;
    break;
  }
  case TimeStepSolverType::_static: {
    options.non_linear_solver_type = NonLinearSolverType::_newton_raphson;
    options.integration_scheme_type["temperature"] =
        IntegrationSchemeType::_pseudo_time;
    options.solution_type["temperature"] = IntegrationScheme::_not_defined;
    break;
  }
  case TimeStepSolverType::_dynamic: {
    if (this->method == _explicit_consistent_mass) {
      options.non_linear_solver_type = NonLinearSolverType::_newton_raphson;
      options.integration_scheme_type["temperature"] =
          IntegrationSchemeType::_forward_euler;
      options.solution_type["temperature"] =
          IntegrationScheme::_temperature_rate;
    } else {
      options.non_linear_solver_type = NonLinearSolverType::_newton_raphson;
      options.integration_scheme_type["temperature"] =
          IntegrationSchemeType::_backward_euler;
      options.solution_type["temperature"] = IntegrationScheme::_temperature;
    }
    break;
  }
  default:
    AKANTU_EXCEPTION(type << " is not a valid time step solver type");
  }

  return options;
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::assembleConductivityMatrix() {
  AKANTU_DEBUG_IN();

  this->computeConductivityOnQuadPoints(_not_ghost);

  if (conductivity_release[_not_ghost] == conductivity_matrix_release) {
    return;
  }

  AKANTU_DEBUG_ASSERT(this->getDOFManager().hasMatrix("K"),
                      "The K matrix has not been initialized yet.");

  this->getDOFManager().zeroMatrix("K");

  auto & fem = this->getFEEngine();

  for (auto && type : mesh.elementTypes(spatial_dimension)) {
    auto nb_element = mesh.getNbElement(type);
    auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    auto nb_quadrature_points = fem.getNbIntegrationPoints(type);

    auto bt_d_b = std::make_unique<Array<Real>>(
        nb_element * nb_quadrature_points,
        nb_nodes_per_element * nb_nodes_per_element, "B^t*D*B");

    fem.computeBtDB(conductivity_on_qpoints(type), *bt_d_b, 2, type);

    /// compute @f$ k_e = \int_e \mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
    auto K_e = std::make_unique<Array<Real>>(
        nb_element, nb_nodes_per_element * nb_nodes_per_element, "K_e");

    fem.integrate(*bt_d_b, *K_e, nb_nodes_per_element * nb_nodes_per_element,
                  type);

    this->getDOFManager().assembleElementalMatricesToMatrix(
        "K", "temperature", *K_e, type, _not_ghost, _symmetric);
  }

  conductivity_matrix_release = conductivity_release[_not_ghost];

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::computeConductivityOnQuadPoints(GhostType ghost_type) {
  // if already computed once check if need to compute
  if (not initial_conductivity[ghost_type]) {
    // if temperature did not change, conductivity will not vary
    if (temperature_release == conductivity_release[ghost_type]) {
      return;
    }
    // if conductivity_variation is 0 no need to recompute
    if (conductivity_variation == 0.) {
      return;
    }
  }

  for (auto && type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_regular)) {
    auto & temperature_interpolated = temperature_on_qpoints(type, ghost_type);

    // compute the temperature on quadrature points
    this->getFEEngine().interpolateOnIntegrationPoints(
        *temperature, temperature_interpolated, 1, type, ghost_type);

    auto & cond = conductivity_on_qpoints(type, ghost_type);
    for (auto && tuple :
         zip(make_view(cond, spatial_dimension, spatial_dimension),
             temperature_interpolated)) {
      auto & C = std::get<0>(tuple);
      auto & T = std::get<1>(tuple);
      C = conductivity;

      Matrix<Real> variation(spatial_dimension, spatial_dimension,
                             conductivity_variation * (T - T_ref));
      // @TODO: Guillaume are you sure ? why due you compute variation then ?
      C += conductivity_variation;
    }
  }

  conductivity_release[ghost_type] = temperature_release;
  initial_conductivity[ghost_type] = false;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::computeKgradT(GhostType ghost_type) {
  computeConductivityOnQuadPoints(ghost_type);

  for (auto && type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_regular)) {
    auto & gradient = temperature_gradient(type, ghost_type);
    this->getFEEngine().gradientOnIntegrationPoints(*temperature, gradient, 1,
                                                    type, ghost_type);

    for (auto && values :
         zip(make_view(conductivity_on_qpoints(type, ghost_type),
                       spatial_dimension, spatial_dimension),
             make_view(gradient, spatial_dimension),
             make_view(k_gradt_on_qpoints(type, ghost_type),
                       spatial_dimension))) {
      const auto & C = std::get<0>(values);
      const auto & BT = std::get<1>(values);
      auto & k_BT = std::get<2>(values);

      k_BT.mul<false>(C, BT);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::assembleInternalHeatRate() {
  AKANTU_DEBUG_IN();

  this->internal_heat_rate->zero();

  this->synchronize(SynchronizationTag::_htm_temperature);
  auto & fem = this->getFEEngine();

  for (auto ghost_type : ghost_types) {
    // compute k \grad T
    computeKgradT(ghost_type);

    for (auto type :
         mesh.elementTypes(spatial_dimension, ghost_type, _ek_regular)) {
      UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

      auto & k_gradt_on_qpoints_vect = k_gradt_on_qpoints(type, ghost_type);

      UInt nb_quad_points = k_gradt_on_qpoints_vect.size();
      Array<Real> bt_k_gT(nb_quad_points, nb_nodes_per_element);
      fem.computeBtD(k_gradt_on_qpoints_vect, bt_k_gT, type, ghost_type);

      UInt nb_elements = mesh.getNbElement(type, ghost_type);
      Array<Real> int_bt_k_gT(nb_elements, nb_nodes_per_element);

      fem.integrate(bt_k_gT, int_bt_k_gT, nb_nodes_per_element, type,
                    ghost_type);

      this->getDOFManager().assembleElementalArrayLocalArray(
          int_bt_k_gT, *this->internal_heat_rate, type, ghost_type, -1);
    }
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Real HeatTransferModel::getStableTimeStep() {
  AKANTU_DEBUG_IN();

  Real el_size;
  Real min_el_size = std::numeric_limits<Real>::max();
  Real conductivitymax = conductivity(0, 0);

  // get the biggest parameter from k11 until k33//
  for (UInt i = 0; i < spatial_dimension; i++) {
    for (UInt j = 0; j < spatial_dimension; j++) {
      conductivitymax = std::max(conductivity(i, j), conductivitymax);
    }
  }
  for (auto && type :
       mesh.elementTypes(spatial_dimension, _not_ghost, _ek_regular)) {

    UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);

    Array<Real> coord(0, nb_nodes_per_element * spatial_dimension);
    FEEngine::extractNodalToElementField(mesh, mesh.getNodes(), coord, type,
                                         _not_ghost);

    auto el_coord = coord.begin(spatial_dimension, nb_nodes_per_element);
    UInt nb_element = mesh.getNbElement(type);

    for (UInt el = 0; el < nb_element; ++el, ++el_coord) {
      el_size = getFEEngine().getElementInradius(*el_coord, type);
      min_el_size = std::min(min_el_size, el_size);
    }

    AKANTU_DEBUG_INFO("The minimum element size : "
                      << min_el_size
                      << " and the max conductivity is : " << conductivitymax);
  }

  Real min_dt = 2. * min_el_size * min_el_size / 4. * density * capacity /
                conductivitymax;

  mesh.getCommunicator().allReduce(min_dt, SynchronizerOperation::_min);

  AKANTU_DEBUG_OUT();

  return min_dt;
}
/* -------------------------------------------------------------------------- */

void HeatTransferModel::setTimeStep(Real time_step, const ID & solver_id) {
  Model::setTimeStep(time_step, solver_id);

#if defined(AKANTU_USE_IOHELPER)
  this->mesh.getDumper("heat_transfer").setTimeStep(time_step);
#endif
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::readMaterials() {
  auto sect = this->getParserSection();

  if (not std::get<1>(sect)) {
    const auto & section = std::get<0>(sect);
    this->parseSection(section);
  }

  conductivity_on_qpoints.set(conductivity);
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::initFullImpl(const ModelOptions & options) {
  Model::initFullImpl(options);

  readMaterials();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::assembleCapacity() {
  AKANTU_DEBUG_IN();
  auto ghost_type = _not_ghost;

  this->getDOFManager().zeroMatrix("M");

  auto & fem = getFEEngineClass<FEEngineType>();

  heat_transfer::details::ComputeRhoFunctor rho_functor(*this);

  for (auto && type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_regular)) {
    fem.assembleFieldMatrix(rho_functor, "M", "temperature",
                            this->getDOFManager(), type, ghost_type);
  }

  need_to_reassemble_capacity = false;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::computeRho(Array<Real> & rho, ElementType type,
                                   GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  FEEngine & fem = this->getFEEngine();
  UInt nb_element = mesh.getNbElement(type, ghost_type);
  UInt nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);

  rho.resize(nb_element * nb_quadrature_points);
  rho.set(this->capacity);

  // Real * rho_1_val = rho.storage();
  // /// compute @f$ rho @f$ for each nodes of each element
  // for (UInt el = 0; el < nb_element; ++el) {
  //   for (UInt n = 0; n < nb_quadrature_points; ++n) {
  //     *rho_1_val++ = this->capacity;
  //   }
  // }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Real HeatTransferModel::computeThermalEnergyByNode() {
  AKANTU_DEBUG_IN();

  Real ethermal = 0.;

  for (auto && pair : enumerate(make_view(
           *internal_heat_rate, internal_heat_rate->getNbComponent()))) {
    auto n = std::get<0>(pair);
    auto & heat_rate = std::get<1>(pair);

    Real heat = 0.;
    bool is_local_node = mesh.isLocalOrMasterNode(n);
    bool count_node = is_local_node;

    for (UInt i = 0; i < heat_rate.size(); ++i) {
      if (count_node) {
        heat += heat_rate[i] * time_step;
      }
    }
    ethermal += heat;
  }

  mesh.getCommunicator().allReduce(ethermal, SynchronizerOperation::_sum);

  AKANTU_DEBUG_OUT();
  return ethermal;
}

/* -------------------------------------------------------------------------- */
template <class iterator>
void HeatTransferModel::getThermalEnergy(
    iterator Eth, Array<Real>::const_iterator<Real> T_it,
    const Array<Real>::const_iterator<Real> & T_end) const {
  for (; T_it != T_end; ++T_it, ++Eth) {
    *Eth = capacity * density * *T_it;
  }
}

/* -------------------------------------------------------------------------- */
Real HeatTransferModel::getThermalEnergy(ElementType type, UInt index) {
  AKANTU_DEBUG_IN();

  UInt nb_quadrature_points = getFEEngine().getNbIntegrationPoints(type);
  Vector<Real> Eth_on_quarature_points(nb_quadrature_points);

  auto T_it = this->temperature_on_qpoints(type).begin();
  T_it += index * nb_quadrature_points;

  auto T_end = T_it + nb_quadrature_points;

  getThermalEnergy(Eth_on_quarature_points.storage(), T_it, T_end);

  return getFEEngine().integrate(Eth_on_quarature_points, type, index);
}

/* -------------------------------------------------------------------------- */
Real HeatTransferModel::getThermalEnergy() {
  Real Eth = 0;

  auto & fem = getFEEngine();

  for (auto && type :
       mesh.elementTypes(spatial_dimension, _not_ghost, _ek_regular)) {
    auto nb_element = mesh.getNbElement(type, _not_ghost);
    auto nb_quadrature_points = fem.getNbIntegrationPoints(type, _not_ghost);
    Array<Real> Eth_per_quad(nb_element * nb_quadrature_points, 1);

    auto & temperature_interpolated = temperature_on_qpoints(type);

    // compute the temperature on quadrature points
    this->getFEEngine().interpolateOnIntegrationPoints(
        *temperature, temperature_interpolated, 1, type);

    auto T_it = temperature_interpolated.begin();
    auto T_end = temperature_interpolated.end();
    getThermalEnergy(Eth_per_quad.begin(), T_it, T_end);

    Eth += fem.integrate(Eth_per_quad, type);
  }

  return Eth;
}

/* -------------------------------------------------------------------------- */
Real HeatTransferModel::getEnergy(const std::string & id) {
  AKANTU_DEBUG_IN();
  Real energy = 0;

  if (id == "thermal") {
    energy = getThermalEnergy();
  }
  // reduction sum over all processors
  mesh.getCommunicator().allReduce(energy, SynchronizerOperation::_sum);

  AKANTU_DEBUG_OUT();
  return energy;
}

/* -------------------------------------------------------------------------- */
Real HeatTransferModel::getEnergy(const std::string & id, ElementType type,
                                  UInt index) {
  AKANTU_DEBUG_IN();

  Real energy = 0.;

  if (id == "thermal") {
    energy = getThermalEnergy(type, index);
  }

  AKANTU_DEBUG_OUT();
  return energy;
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER

std::shared_ptr<dumpers::Field> HeatTransferModel::createNodalFieldBool(
    const std::string & field_name, const std::string & group_name,
    __attribute__((unused)) bool padding_flag) {

  std::map<std::string, Array<bool> *> uint_nodal_fields;
  uint_nodal_fields["blocked_dofs"] = blocked_dofs.get();

  auto field = mesh.createNodalField(uint_nodal_fields[field_name], group_name);
  return field;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field> HeatTransferModel::createNodalFieldReal(
    const std::string & field_name, const std::string & group_name,
    __attribute__((unused)) bool padding_flag) {

  if (field_name == "capacity_lumped") {
    AKANTU_EXCEPTION(
        "Capacity lumped is a nodal field now stored in the DOF manager."
        "Therefore it cannot be used by a dumper anymore");
  }

  std::map<std::string, Array<Real> *> real_nodal_fields;
  real_nodal_fields["temperature"] = temperature.get();
  real_nodal_fields["temperature_rate"] = temperature_rate.get();
  real_nodal_fields["external_heat_rate"] = external_heat_rate.get();
  real_nodal_fields["internal_heat_rate"] = internal_heat_rate.get();
  real_nodal_fields["increment"] = increment.get();

  std::shared_ptr<dumpers::Field> field =
      mesh.createNodalField(real_nodal_fields[field_name], group_name);

  return field;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field> HeatTransferModel::createElementalField(
    const std::string & field_name, const std::string & group_name,
    bool /*padding_flag*/, UInt /*spatial_dimension*/,
    ElementKind element_kind) {

  std::shared_ptr<dumpers::Field> field;

  if (field_name == "partitions") {
    field = mesh.createElementalField<UInt, dumpers::ElementPartitionField>(
        mesh.getConnectivities(), group_name, this->spatial_dimension,
        element_kind);
  } else if (field_name == "temperature_gradient") {
    ElementTypeMap<UInt> nb_data_per_elem =
        this->mesh.getNbDataPerElem(temperature_gradient);

    field = mesh.createElementalField<Real, dumpers::InternalMaterialField>(
        temperature_gradient, group_name, this->spatial_dimension, element_kind,
        nb_data_per_elem);
  } else if (field_name == "conductivity") {
    ElementTypeMap<UInt> nb_data_per_elem =
        this->mesh.getNbDataPerElem(conductivity_on_qpoints);

    field = mesh.createElementalField<Real, dumpers::InternalMaterialField>(
        conductivity_on_qpoints, group_name, this->spatial_dimension,
        element_kind, nb_data_per_elem);
  }

  return field;
}

/* -------------------------------------------------------------------------- */
#else
/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field> HeatTransferModel::createElementalField(
    const std::string & /* field_name*/, const std::string & /*group_name*/,
    bool /*padding_flag*/, ElementKind /*element_kind*/) {
  return nullptr;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field>
HeatTransferModel::createNodalFieldBool(const std::string & /*field_name*/,
                                        const std::string & /*group_name*/,
                                        bool /*padding_flag*/) {
  return nullptr;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field>
HeatTransferModel::createNodalFieldReal(const std::string & /*field_name*/,
                                        const std::string & /*group_name*/,
                                        bool /*padding_flag*/) {
  return nullptr;
}
#endif

/* -------------------------------------------------------------------------- */
inline UInt HeatTransferModel::getNbData(const Array<UInt> & indexes,
                                         const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  UInt nb_nodes = indexes.size();

  switch (tag) {
  case SynchronizationTag::_htm_temperature: {
    size += nb_nodes * sizeof(Real);
    break;
  }
  default: {
    AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
inline void HeatTransferModel::packData(CommunicationBuffer & buffer,
                                        const Array<UInt> & indexes,
                                        const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  for (auto index : indexes) {
    switch (tag) {
    case SynchronizationTag::_htm_temperature: {
      buffer << (*temperature)(index);
      break;
    }
    default: {
      AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
    }
    }
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void HeatTransferModel::unpackData(CommunicationBuffer & buffer,
                                          const Array<UInt> & indexes,
                                          const SynchronizationTag & tag) {
  AKANTU_DEBUG_IN();

  for (auto index : indexes) {
    switch (tag) {
    case SynchronizationTag::_htm_temperature: {
      buffer >> (*temperature)(index);
      break;
    }
    default: {
      AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
    }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline UInt HeatTransferModel::getNbData(const Array<Element> & elements,
                                         const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  UInt nb_nodes_per_element = 0;
  Array<Element>::const_iterator<Element> it = elements.begin();
  Array<Element>::const_iterator<Element> end = elements.end();
  for (; it != end; ++it) {
    const Element & el = *it;
    nb_nodes_per_element += Mesh::getNbNodesPerElement(el.type);
  }

  switch (tag) {
  case SynchronizationTag::_htm_temperature: {
    size += nb_nodes_per_element * sizeof(Real); // temperature
    break;
  }
  case SynchronizationTag::_htm_gradient_temperature: {
    // temperature gradient
    size += getNbIntegrationPoints(elements) * spatial_dimension * sizeof(Real);
    size += nb_nodes_per_element * sizeof(Real); // nodal temperatures
    break;
  }
  default: {
    AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
inline void HeatTransferModel::packData(CommunicationBuffer & buffer,
                                        const Array<Element> & elements,
                                        const SynchronizationTag & tag) const {
  switch (tag) {
  case SynchronizationTag::_htm_temperature: {
    packNodalDataHelper(*temperature, buffer, elements, mesh);
    break;
  }
  case SynchronizationTag::_htm_gradient_temperature: {
    packElementalDataHelper(temperature_gradient, buffer, elements, true,
                            getFEEngine());
    packNodalDataHelper(*temperature, buffer, elements, mesh);
    break;
  }
  default: {
    AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }
}

/* -------------------------------------------------------------------------- */
inline void HeatTransferModel::unpackData(CommunicationBuffer & buffer,
                                          const Array<Element> & elements,
                                          const SynchronizationTag & tag) {
  switch (tag) {
  case SynchronizationTag::_htm_temperature: {
    unpackNodalDataHelper(*temperature, buffer, elements, mesh);
    break;
  }
  case SynchronizationTag::_htm_gradient_temperature: {
    unpackElementalDataHelper(temperature_gradient, buffer, elements, true,
                              getFEEngine());
    unpackNodalDataHelper(*temperature, buffer, elements, mesh);

    break;
  }
  default: {
    AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }
}

/* -------------------------------------------------------------------------- */
} // namespace akantu
