/**
 * @file   heat_transfer_model.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Srinivasa Babu Ramisetti <srinivasa.ramisetti@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Rui Wang <rui.wang@epfl.ch>
 *
 * @date creation: Sun May 01 2011
 * @date last modification: Mon Nov 30 2015
 *
 * @brief  Implementation of HeatTransferModel class
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
#include "heat_transfer_model.hh"
#include "dumpable_inline_impl.hh"
#include "element_synchronizer.hh"
#include "fe_engine_template.hh"
#include "generalized_trapezoidal.hh"
#include "group_manager_inline_impl.cc"
#include "mesh.hh"
#include "parser.hh"

#ifdef AKANTU_USE_IOHELPER
#include "dumper_element_partition.hh"
#include "dumper_elemental_field.hh"
#include "dumper_internal_material_field.hh"
#include "dumper_paraview.hh"
#endif

/* -------------------------------------------------------------------------- */
namespace akantu {

class ComputeRhoFunctor {
public:
  ComputeRhoFunctor(const HeatTransferModel & model) : model(model){};

  void operator()(Matrix<Real> & rho, const Element &) const {
    rho.set(model.getCapacity());
  }

private:
  const HeatTransferModel & model;
};

/* -------------------------------------------------------------------------- */
HeatTransferModel::HeatTransferModel(Mesh & mesh, UInt dim, const ID & id,
                                     const MemoryID & memory_id)
    : Model(mesh, ModelType::_heat_transfer_model, dim, id, memory_id),
      integrator(nullptr), temperature_gradient("temperature_gradient", id),
      temperature_on_qpoints("temperature_on_qpoints", id),
      conductivity_on_qpoints("conductivity_on_qpoints", id),
      k_gradt_on_qpoints("k_gradt_on_qpoints", id), conductivity(dim, dim) {
  AKANTU_DEBUG_IN();

  this->initDOFManager();

  this->registerDataAccessor(*this);

  if (this->mesh.isDistributed()) {
    auto & synchronizer = this->mesh.getElementSynchronizer();
    this->registerSynchronizer(synchronizer, _gst_htm_capacity);
    this->registerSynchronizer(synchronizer, _gst_htm_temperature);
    this->registerSynchronizer(synchronizer, _gst_htm_gradient_temperature);
  }

  std::stringstream sstr;
  sstr << id << ":fem";
  registerFEEngineObject<MyFEEngineType>(sstr.str(), mesh, spatial_dimension);

  this->temperature = nullptr;
  this->internal_heat_rate = nullptr;
  this->blocked_dofs = nullptr;

#ifdef AKANTU_USE_IOHELPER
  this->mesh.registerDumper<DumperParaview>("paraview_all", id, true);
  this->mesh.addDumpMesh(mesh, spatial_dimension, _not_ghost, _ek_regular);
#endif

  this->registerParam("conductivity", conductivity, _pat_parsmod);
  this->registerParam("conductivity_variation", conductivity_variation, 0.,
                      _pat_parsmod);
  this->registerParam("temperature_reference", T_ref, 0., _pat_parsmod);
  this->registerParam("capacity", capacity, _pat_parsmod);
  // this->registerParam("density", density, _pat_parsmod);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::initModel() {
  getFEEngine().initShapeFunctions(_not_ghost);
  getFEEngine().initShapeFunctions(_ghost);
}

/* -------------------------------------------------------------------------- */
FEEngine & HeatTransferModel::getFEEngineBoundary(const ID & name) {
  return dynamic_cast<FEEngine &>(
      getFEEngineClassBoundary<MyFEEngineType>(name));
}



/* -------------------------------------------------------------------------- */
template <typename T>
void HeatTransferModel::allocNodalField(Array<T> *& array, const ID & name) {
  if (array == nullptr) {
    UInt nb_nodes = mesh.getNbNodes();
    std::stringstream sstr_disp;
    sstr_disp << id << ":" << name;

    array = &(alloc<T>(sstr_disp.str(), nb_nodes, 1, T()));
  }
}

/* -------------------------------------------------------------------------- */
HeatTransferModel::~HeatTransferModel() = default;

/* -------------------------------------------------------------------------- */
void HeatTransferModel::assembleCapacityLumped(const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  auto & fem = getFEEngineClass<MyFEEngineType>();
  ComputeRhoFunctor compute_rho(*this);

  for (auto & type : mesh.elementTypes(spatial_dimension, ghost_type)) {

    fem.assembleFieldLumped(compute_rho, "M", "temperature",
                            this->getDOFManager(), type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::assembleCapacityLumped() {
  AKANTU_DEBUG_IN();

  if (!this->getDOFManager().hasLumpedMatrix("M")) {
    this->getDOFManager().getNewLumpedMatrix("M");
  }

  this->getDOFManager().clearLumpedMatrix("M");

  assembleCapacityLumped(_not_ghost);
  assembleCapacityLumped(_ghost);

  AKANTU_DEBUG_OUT();
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
void HeatTransferModel::initSolver(TimeStepSolverType time_step_solver_type,
                                   NonLinearSolverType) {

  DOFManager & dof_manager = this->getDOFManager();

  this->allocNodalField(this->temperature, "temperature");
  this->allocNodalField(this->external_heat_rate, "external_heat_rate");
  this->allocNodalField(this->internal_heat_rate, "internal_heat_rate");
  this->allocNodalField(this->blocked_dofs, "blocked_dofs");

  if (!dof_manager.hasDOFs("temperature")) {
    dof_manager.registerDOFs("temperature", *this->temperature, _dst_nodal);
    dof_manager.registerBlockedDOFs("temperature", *this->blocked_dofs);
  }

  if (time_step_solver_type == _tsst_dynamic ||
      time_step_solver_type == _tsst_dynamic_lumped) {

    this->allocNodalField(this->temperature_rate, "temperature_rate");

    if (!dof_manager.hasDOFsDerivatives("temperature", 1)) {
      dof_manager.registerDOFsDerivative("temperature", 1,
                                         *this->temperature_rate);
    }
  }

  /* --------------------------------------------------------------------------
   */
  // byelementtype vectors
  temperature_on_qpoints.initialize(this->getFEEngine(), _nb_component = 1);

  temperature_gradient.initialize(this->getFEEngine(),
                                  _nb_component = spatial_dimension);

  conductivity_on_qpoints.initialize(this->getFEEngine(),
                                     _nb_component =
                                         spatial_dimension * spatial_dimension);

  k_gradt_on_qpoints.initialize(this->getFEEngine(),
                                _nb_component = spatial_dimension);
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::assembleJacobian() {
  AKANTU_DEBUG_IN();

  if (!this->getDOFManager().hasLumpedMatrix("K")) {
    this->getDOFManager().getNewLumpedMatrix("K");
  }
  this->getDOFManager().clearLumpedMatrix("K");

  switch (mesh.getSpatialDimension()) {
  case 1:
    this->assembleConductivityMatrix<1>(_not_ghost);
    break;
  case 2:
    this->assembleConductivityMatrix<2>(_not_ghost);
    break;
  case 3:
    this->assembleConductivityMatrix<3>(_not_ghost);
    break;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void HeatTransferModel::assembleConductivityMatrix(const GhostType & ghost_type,
                                                   bool compute_conductivity) {
  AKANTU_DEBUG_IN();

  Mesh & mesh = this->getFEEngine().getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type);
  for (; it != last_type; ++it) {
    this->assembleConductivityMatrix<dim>(*it, ghost_type,
                                          compute_conductivity);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void HeatTransferModel::assembleConductivityMatrix(const ElementType & type,
                                                   const GhostType & ghost_type,
                                                   bool compute_conductivity) {
  AKANTU_DEBUG_IN();

  const Array<Real> & shapes_derivatives =
      this->getFEEngine().getShapesDerivatives(type, ghost_type);

  UInt nb_element = mesh.getNbElement(type, ghost_type);
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points =
      getFEEngine().getNbIntegrationPoints(type, ghost_type);

  /// compute @f$\mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  UInt bt_d_b_size = nb_nodes_per_element;

  Array<Real> * bt_d_b = new Array<Real>(nb_element * nb_quadrature_points,
                                         bt_d_b_size * bt_d_b_size, "B^t*D*B");

  Matrix<Real> Bt_D(nb_nodes_per_element, dim);

  Array<Real>::const_iterator<Matrix<Real>> shapes_derivatives_it =
      shapes_derivatives.begin(dim, nb_nodes_per_element);

  auto Bt_D_B_it = bt_d_b->begin(bt_d_b_size, bt_d_b_size);

  if (compute_conductivity)
    this->computeConductivityOnQuadPoints(ghost_type);
  auto D_it = conductivity_on_qpoints(type, ghost_type).begin(dim, dim);
  auto D_end = conductivity_on_qpoints(type, ghost_type).end(dim, dim);

  for (; D_it != D_end; ++D_it, ++Bt_D_B_it, ++shapes_derivatives_it) {
    Matrix<Real> & D = *D_it;
    const Matrix<Real> & B = *shapes_derivatives_it;
    Matrix<Real> & Bt_D_B = *Bt_D_B_it;

    Bt_D.mul<true, false>(B, D);
    Bt_D_B.mul<false, false>(Bt_D, B);
  }

  /// compute @f$ k_e = \int_e \mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  Array<Real> * K_e =
      new Array<Real>(nb_element, bt_d_b_size * bt_d_b_size, "K_e");

  this->getFEEngine().integrate(*bt_d_b, *K_e, bt_d_b_size * bt_d_b_size, type,
                                ghost_type);

  delete bt_d_b;

  this->getDOFManager().assembleElementalMatricesToMatrix(
      "K", "temperature", *K_e, type, ghost_type, _symmetric);
  delete K_e;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::computeConductivityOnQuadPoints(
    const GhostType & ghost_type) {

  for (auto & type : mesh.elementTypes(spatial_dimension, ghost_type)) {

    Array<Real> & temperature_interpolated =
        temperature_on_qpoints(type, ghost_type);

    // compute the temperature on quadrature points
    this->getFEEngine().interpolateOnIntegrationPoints(
        *temperature, temperature_interpolated, 1, type, ghost_type);

    auto C_it = conductivity_on_qpoints(type, ghost_type)
                    .begin(spatial_dimension, spatial_dimension);
    auto C_end = conductivity_on_qpoints(type, ghost_type)
                     .end(spatial_dimension, spatial_dimension);

    auto T_it = temperature_interpolated.begin();

    for (; C_it != C_end; ++C_it, ++T_it) {
      Matrix<Real> & C = *C_it;
      Real & T = *T_it;
      C = conductivity;

      Matrix<Real> variation(spatial_dimension, spatial_dimension,
                             conductivity_variation * (T - T_ref));
      C += conductivity_variation;
    }
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::computeKgradT(const GhostType & ghost_type) {

  if (this->compute_conductivity)
    computeConductivityOnQuadPoints(ghost_type);

  for (auto & type : mesh.elementTypes(spatial_dimension, ghost_type)) {

    Array<Real> & gradient = temperature_gradient(type, ghost_type);
    this->getFEEngine().gradientOnIntegrationPoints(*temperature, gradient, 1,
                                                    type, ghost_type);

    Array<Real>::matrix_iterator C_it =
        conductivity_on_qpoints(type, ghost_type)
            .begin(spatial_dimension, spatial_dimension);

    Array<Real>::vector_iterator BT_it = gradient.begin(spatial_dimension);

    Array<Real>::vector_iterator k_BT_it =
        k_gradt_on_qpoints(type, ghost_type).begin(spatial_dimension);
    Array<Real>::vector_iterator k_BT_end =
        k_gradt_on_qpoints(type, ghost_type).end(spatial_dimension);

    for (; k_BT_it != k_BT_end; ++k_BT_it, ++BT_it, ++C_it) {
      Vector<Real> & k_BT = *k_BT_it;
      Vector<Real> & BT = *BT_it;
      Matrix<Real> & C = *C_it;

      k_BT.mul<false>(C, BT);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::assembleInternalHeatRate() {
  AKANTU_DEBUG_IN();

  this->synchronize(_gst_htm_temperature);

  for (auto ghost_type : ghost_types) {
    for (auto type : mesh.elementTypes(spatial_dimension, ghost_type)) {

      Array<Real> & shapes_derivatives = const_cast<Array<Real> &>(
          getFEEngine().getShapesDerivatives(type, ghost_type));

      UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

      // compute k \grad T
      computeKgradT(ghost_type);
      auto & k_gradt_on_qpoints_vect = k_gradt_on_qpoints(type, ghost_type);
      auto k_BT_it = k_gradt_on_qpoints_vect.begin(spatial_dimension);
      auto B_it =
          shapes_derivatives.begin(spatial_dimension, nb_nodes_per_element);

      UInt nb_quad_points = k_gradt_on_qpoints_vect.size();
      Array<Real> bt_k_gT(nb_quad_points, nb_nodes_per_element);

      auto Bt_k_BT_it = bt_k_gT.begin(nb_nodes_per_element);
      auto Bt_k_BT_end = bt_k_gT.end(nb_nodes_per_element);

      for (; Bt_k_BT_it != Bt_k_BT_end; ++Bt_k_BT_it, ++B_it, ++k_BT_it) {
        Vector<Real> & k_BT = *k_BT_it;
        Vector<Real> & Bt_k_BT = *Bt_k_BT_it;
        Matrix<Real> & B = *B_it;

        Bt_k_BT.mul<true>(B, k_BT);
      }

      UInt nb_elements = mesh.getNbElement(type, ghost_type);
      Array<Real> int_bt_k_gT(nb_elements, nb_nodes_per_element);

      this->getFEEngine().integrate(bt_k_gT, int_bt_k_gT, nb_nodes_per_element,
                                    type, ghost_type);

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
  for (UInt i = 0; i < spatial_dimension; i++)
    for (UInt j = 0; j < spatial_dimension; j++)
      conductivitymax = std::max(conductivity(i, j), conductivitymax);

  for (auto & type : mesh.elementTypes(spatial_dimension, _not_ghost)) {

    UInt nb_nodes_per_element =
        getFEEngine().getMesh().getNbNodesPerElement(type);

    Array<Real> coord(0, nb_nodes_per_element * spatial_dimension);
    FEEngine::extractNodalToElementField(getFEEngine().getMesh(),
                                         getFEEngine().getMesh().getNodes(),
                                         coord, type, _not_ghost);

    auto el_coord = coord.begin(spatial_dimension, nb_nodes_per_element);
    UInt nb_element = getFEEngine().getMesh().getNbElement(type);

    for (UInt el = 0; el < nb_element; ++el, ++el_coord) {
      el_size = getFEEngine().getElementInradius(*el_coord, type);
      min_el_size = std::min(min_el_size, el_size);
    }

    AKANTU_DEBUG_INFO("The minimum element size : "
                      << min_el_size
                      << " and the max conductivity is : " << conductivitymax);
  }

  Real min_dt =
      2 * min_el_size * min_el_size * density * capacity / conductivitymax;

  mesh.getCommunicator().allReduce(min_dt, SynchronizerOperation::_min);

  AKANTU_DEBUG_OUT();

  return min_dt;
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::readMaterials() {
  auto sect = this->getParserSection();

  if (std::get<1>(sect)) {
    const auto & section = std::get<0>(sect);
    this->parseSection(section);
  }
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

  auto & fem = getFEEngineClass<MyFEEngineType>();

  ComputeRhoFunctor rho_functor(*this);

  for (auto && type : mesh.elementTypes(spatial_dimension, ghost_type)) {
    fem.assembleFieldMatrix(rho_functor, "M", "temperature",
                            this->getDOFManager(), type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::computeRho(Array<Real> & rho, ElementType type,
                                   GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  FEEngine & fem = this->getFEEngine();
  UInt nb_element = fem.getMesh().getNbElement(type, ghost_type);
  UInt nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);

  rho.resize(nb_element * nb_quadrature_points);
  Real * rho_1_val = rho.storage();

  /// compute @f$ rho @f$ for each nodes of each element
  for (UInt el = 0; el < nb_element; ++el) {

    for (UInt n = 0; n < nb_quadrature_points; ++n) {
      *rho_1_val++ = this->capacity;
    }
  }

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
    bool is_not_pbc_slave_node = !isPBCSlaveNode(n);
    bool count_node = is_local_node && is_not_pbc_slave_node;

    for (UInt i = 0; i < heat_rate.size(); ++i) {
      if (count_node)
        heat += heat_rate[i] * time_step;
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
    Array<Real>::const_iterator<Real> T_end) const {
  for (; T_it != T_end; ++T_it, ++Eth) {
    *Eth = capacity * density * *T_it;
  }
}

/* -------------------------------------------------------------------------- */
Real HeatTransferModel::getThermalEnergy(const ElementType & type, UInt index) {
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

  Mesh & mesh = getFEEngine().getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension);
  for (; it != last_type; ++it) {
    UInt nb_element = getFEEngine().getMesh().getNbElement(*it, _not_ghost);
    UInt nb_quadrature_points =
        getFEEngine().getNbIntegrationPoints(*it, _not_ghost);
    Array<Real> Eth_per_quad(nb_element * nb_quadrature_points, 1);

    auto T_it = this->temperature_on_qpoints(*it).begin();
    auto T_end = this->temperature_on_qpoints(*it).end();
    getThermalEnergy(Eth_per_quad.begin(), T_it, T_end);

    Eth += getFEEngine().integrate(Eth_per_quad, *it);
  }

  return Eth;
}

/* -------------------------------------------------------------------------- */
Real HeatTransferModel::getEnergy(const std::string & id) {
  AKANTU_DEBUG_IN();
  Real energy = 0;

  if (id == "thermal")
    energy = getThermalEnergy();

  // reduction sum over all processors
  mesh.getCommunicator().allReduce(energy, SynchronizerOperation::_sum);

  AKANTU_DEBUG_OUT();
  return energy;
}

/* -------------------------------------------------------------------------- */
Real HeatTransferModel::getEnergy(const std::string & id,
                                  const ElementType & type, UInt index) {
  AKANTU_DEBUG_IN();

  Real energy = 0.;

  if (id == "thermal")
    energy = getThermalEnergy(type, index);

  AKANTU_DEBUG_OUT();
  return energy;
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER

dumper::Field * HeatTransferModel::createNodalFieldBool(
    const std::string & field_name, const std::string & group_name,
    __attribute__((unused)) bool padding_flag) {

  std::map<std::string, Array<bool> *> uint_nodal_fields;
  uint_nodal_fields["blocked_dofs"] = blocked_dofs;

  dumper::Field * field = nullptr;
  field = mesh.createNodalField(uint_nodal_fields[field_name], group_name);
  return field;
}

/* -------------------------------------------------------------------------- */
dumper::Field * HeatTransferModel::createNodalFieldReal(
    const std::string & field_name, const std::string & group_name,
    __attribute__((unused)) bool padding_flag) {

  std::map<std::string, Array<Real> *> real_nodal_fields;
  real_nodal_fields["temperature"] = temperature;
  real_nodal_fields["temperature_rate"] = temperature_rate;
  real_nodal_fields["external_heat_rate"] = external_heat_rate;
  real_nodal_fields["internal_heat_rate"] = internal_heat_rate;
  real_nodal_fields["capacity_lumped"] = capacity_lumped;
  real_nodal_fields["increment"] = increment;

  dumper::Field * field =
      mesh.createNodalField(real_nodal_fields[field_name], group_name);

  return field;
}

/* -------------------------------------------------------------------------- */
dumper::Field * HeatTransferModel::createElementalField(
    const std::string & field_name, const std::string & group_name,
    __attribute__((unused)) bool padding_flag,
    __attribute__((unused)) const UInt & spatial_dimension,
    const ElementKind & element_kind) {

  dumper::Field * field = nullptr;

  if (field_name == "partitions")
    field = mesh.createElementalField<UInt, dumper::ElementPartitionField>(
        mesh.getConnectivities(), group_name, this->spatial_dimension,
        element_kind);
  else if (field_name == "temperature_gradient") {
    ElementTypeMap<UInt> nb_data_per_elem =
        this->mesh.getNbDataPerElem(temperature_gradient, element_kind);

    field = mesh.createElementalField<Real, dumper::InternalMaterialField>(
        temperature_gradient, group_name, this->spatial_dimension, element_kind,
        nb_data_per_elem);
  } else if (field_name == "conductivity") {
    ElementTypeMap<UInt> nb_data_per_elem =
        this->mesh.getNbDataPerElem(conductivity_on_qpoints, element_kind);

    field = mesh.createElementalField<Real, dumper::InternalMaterialField>(
        conductivity_on_qpoints, group_name, this->spatial_dimension,
        element_kind, nb_data_per_elem);
  }

  return field;
}

/* -------------------------------------------------------------------------- */
#else
/* -------------------------------------------------------------------------- */
dumper::Field * HeatTransferModel::createElementalField(
    __attribute__((unused)) const std::string & field_name,
    __attribute__((unused)) const std::string & group_name,
    __attribute__((unused)) bool padding_flag,
    __attribute__((unused)) const ElementKind & element_kind) {
  return nullptr;
}

/* -------------------------------------------------------------------------- */
dumper::Field * HeatTransferModel::createNodalFieldBool(
    __attribute__((unused)) const std::string & field_name,
    __attribute__((unused)) const std::string & group_name,
    __attribute__((unused)) bool padding_flag) {
  return nullptr;
}

/* -------------------------------------------------------------------------- */
dumper::Field * HeatTransferModel::createNodalFieldReal(
    __attribute__((unused)) const std::string & field_name,
    __attribute__((unused)) const std::string & group_name,
    __attribute__((unused)) bool padding_flag) {
  return nullptr;
}

#endif

} // akantu
