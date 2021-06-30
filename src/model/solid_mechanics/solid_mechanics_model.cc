/**
 * @file   solid_mechanics_model.cc
 *
 * @author Ramin Aghababaei <ramin.aghababaei@epfl.ch>
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Clement Roux <clement.roux@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Tue Jul 27 2010
 * @date last modification: Wed Feb 21 2018
 *
 * @brief  Implementation of the SolidMechanicsModel class
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
#include "solid_mechanics_model.hh"
#include "integrator_gauss.hh"
#include "shape_lagrange.hh"
#include "solid_mechanics_model_tmpl.hh"

#include "communicator.hh"
#include "element_synchronizer.hh"
#include "sparse_matrix.hh"
#include "synchronizer_registry.hh"

#include "dumpable_inline_impl.hh"
#ifdef AKANTU_USE_IOHELPER
#include "dumper_iohelper_paraview.hh"
#endif

#include "material_non_local.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
/**
 * A solid mechanics model need a mesh  and a dimension to be created. the model
 * by it  self can not  do a lot,  the good init  functions should be  called in
 * order to configure the model depending on what we want to do.
 *
 * @param  mesh mesh  representing  the model  we  want to  simulate
 * @param dim spatial  dimension of the problem, if dim =  0 (default value) the
 * dimension of the problem is assumed to be the on of the mesh
 * @param id an id to identify the model
 * @param model_type this is an internal parameter for inheritance purposes
 */
SolidMechanicsModel::SolidMechanicsModel(Mesh & mesh, UInt dim, const ID & id,
                                         const ModelType model_type)
    : Model(mesh, model_type, dim, id),
      material_index("material index", id),
      material_local_numbering("material local numbering", id) {
  AKANTU_DEBUG_IN();

  this->registerFEEngineObject<MyFEEngineType>("SolidMechanicsFEEngine", mesh,
                                               Model::spatial_dimension);

#if defined(AKANTU_USE_IOHELPER)
  this->mesh.registerDumper<DumperParaview>("solid_mechanics_model", id, true);
  this->mesh.addDumpMesh(mesh, Model::spatial_dimension, _not_ghost,
                         _ek_regular);
#endif

  material_selector = std::make_shared<DefaultMaterialSelector>(material_index);

  this->initDOFManager();

  this->registerDataAccessor(*this);

  if (this->mesh.isDistributed()) {
    auto & synchronizer = this->mesh.getElementSynchronizer();
    this->registerSynchronizer(synchronizer, SynchronizationTag::_material_id);
    this->registerSynchronizer(synchronizer, SynchronizationTag::_smm_mass);
    this->registerSynchronizer(synchronizer, SynchronizationTag::_smm_stress);
    this->registerSynchronizer(synchronizer, SynchronizationTag::_for_dump);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SolidMechanicsModel::~SolidMechanicsModel() = default;

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::setTimeStep(Real time_step, const ID & solver_id) {
  Model::setTimeStep(time_step, solver_id);

#if defined(AKANTU_USE_IOHELPER)
  this->mesh.getDumper().setTimeStep(time_step);
#endif
}

/* -------------------------------------------------------------------------- */
/* Initialization                                                             */
/* -------------------------------------------------------------------------- */
/**
 * This function groups  many of the initialization in on  function. For most of
 * basics  case the  function should  be  enough. The  functions initialize  the
 * model, the internal  vectors, set them to 0, and  depending on the parameters
 * it also initialize the explicit or implicit solver.
 *
 * @param options
 * \parblock
 * contains the different options to initialize the model
 * \li \c analysis_method specify the type of solver to use
 * \endparblock
 */
void SolidMechanicsModel::initFullImpl(const ModelOptions & options) {
  material_index.initialize(mesh, _element_kind = _ek_not_defined,
                            _default_value = UInt(-1), _with_nb_element = true);
  material_local_numbering.initialize(mesh, _element_kind = _ek_not_defined,
                                      _with_nb_element = true);

  Model::initFullImpl(options);

  // initialize the materials
  if (not this->parser.getLastParsedFile().empty()) {
    this->instantiateMaterials();
    this->initMaterials();
  }

  this->initBC(*this, *displacement, *displacement_increment, *external_force);
}

/* -------------------------------------------------------------------------- */
TimeStepSolverType SolidMechanicsModel::getDefaultSolverType() const {
  return TimeStepSolverType::_dynamic_lumped;
}

/* -------------------------------------------------------------------------- */
ModelSolverOptions SolidMechanicsModel::getDefaultSolverOptions(
    const TimeStepSolverType & type) const {
  ModelSolverOptions options;

  switch (type) {
  case TimeStepSolverType::_dynamic_lumped: {
    options.non_linear_solver_type = NonLinearSolverType::_lumped;
    options.integration_scheme_type["displacement"] =
        IntegrationSchemeType::_central_difference;
    options.solution_type["displacement"] = IntegrationScheme::_acceleration;
    break;
  }
  case TimeStepSolverType::_static: {
    options.non_linear_solver_type = NonLinearSolverType::_newton_raphson;
    options.integration_scheme_type["displacement"] =
        IntegrationSchemeType::_pseudo_time;
    options.solution_type["displacement"] = IntegrationScheme::_not_defined;
    break;
  }
  case TimeStepSolverType::_dynamic: {
    if (this->method == _explicit_consistent_mass) {
      options.non_linear_solver_type = NonLinearSolverType::_newton_raphson;
      options.integration_scheme_type["displacement"] =
          IntegrationSchemeType::_central_difference;
      options.solution_type["displacement"] = IntegrationScheme::_acceleration;
    } else {
      options.non_linear_solver_type = NonLinearSolverType::_newton_raphson;
      options.integration_scheme_type["displacement"] =
          IntegrationSchemeType::_trapezoidal_rule_2;
      options.solution_type["displacement"] = IntegrationScheme::_displacement;
    }
    break;
  }
  default:
    AKANTU_EXCEPTION(type << " is not a valid time step solver type");
  }

  return options;
}

/* -------------------------------------------------------------------------- */
std::tuple<ID, TimeStepSolverType>
SolidMechanicsModel::getDefaultSolverID(const AnalysisMethod & method) {
  switch (method) {
  case _explicit_lumped_mass: {
    return std::make_tuple("explicit_lumped",
                           TimeStepSolverType::_dynamic_lumped);
  }
  case _explicit_consistent_mass: {
    return std::make_tuple("explicit", TimeStepSolverType::_dynamic);
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
void SolidMechanicsModel::initSolver(TimeStepSolverType time_step_solver_type,
                                     NonLinearSolverType /*unused*/) {
  auto & dof_manager = this->getDOFManager();

  /* ------------------------------------------------------------------------ */
  // for alloc type of solvers
  this->allocNodalField(this->displacement, spatial_dimension, "displacement");
  this->allocNodalField(this->previous_displacement, spatial_dimension,
                        "previous_displacement");
  this->allocNodalField(this->displacement_increment, spatial_dimension,
                        "displacement_increment");
  this->allocNodalField(this->internal_force, spatial_dimension,
                        "internal_force");
  this->allocNodalField(this->external_force, spatial_dimension,
                        "external_force");
  this->allocNodalField(this->blocked_dofs, spatial_dimension, "blocked_dofs");
  this->allocNodalField(this->current_position, spatial_dimension,
                        "current_position");

  // initialize the current positions
  this->current_position->copy(this->mesh.getNodes());

  /* ------------------------------------------------------------------------ */
  if (!dof_manager.hasDOFs("displacement")) {
    dof_manager.registerDOFs("displacement", *this->displacement, _dst_nodal);
    dof_manager.registerBlockedDOFs("displacement", *this->blocked_dofs);
    dof_manager.registerDOFsIncrement("displacement",
                                      *this->displacement_increment);
    dof_manager.registerDOFsPrevious("displacement",
                                     *this->previous_displacement);
  }

  /* ------------------------------------------------------------------------ */
  // for dynamic
  if (time_step_solver_type == TimeStepSolverType::_dynamic ||
      time_step_solver_type == TimeStepSolverType::_dynamic_lumped) {
    this->allocNodalField(this->velocity, spatial_dimension, "velocity");
    this->allocNodalField(this->acceleration, spatial_dimension,
                          "acceleration");

    if (!dof_manager.hasDOFsDerivatives("displacement", 1)) {
      dof_manager.registerDOFsDerivative("displacement", 1, *this->velocity);
      dof_manager.registerDOFsDerivative("displacement", 2,
                                         *this->acceleration);
    }
  }
}

/* -------------------------------------------------------------------------- */
/**
 * Initialize the model,basically it  pre-compute the shapes, shapes derivatives
 * and jacobian
 */
void SolidMechanicsModel::initModel() {
  /// \todo add  the current position  as a parameter to  initShapeFunctions for
  /// large deformation
  getFEEngine().initShapeFunctions(_not_ghost);
  getFEEngine().initShapeFunctions(_ghost);
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleResidual() {
  AKANTU_DEBUG_IN();

  /* ------------------------------------------------------------------------ */
  // computes the internal forces
  this->assembleInternalForces();

  /* ------------------------------------------------------------------------ */
  this->getDOFManager().assembleToResidual("displacement",
                                           *this->external_force, 1);
  this->getDOFManager().assembleToResidual("displacement",
                                           *this->internal_force, 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleResidual(const ID & residual_part) {
  AKANTU_DEBUG_IN();

  if ("external" == residual_part) {
    this->getDOFManager().assembleToResidual("displacement",
                                             *this->external_force, 1);
    AKANTU_DEBUG_OUT();
    return;
  }

  if ("internal" == residual_part) {
    this->assembleInternalForces();
    this->getDOFManager().assembleToResidual("displacement",
                                             *this->internal_force, 1);
    AKANTU_DEBUG_OUT();
    return;
  }

  AKANTU_CUSTOM_EXCEPTION(
      debug::SolverCallbackResidualPartUnknown(residual_part));

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
MatrixType SolidMechanicsModel::getMatrixType(const ID & matrix_id) {
  // \TODO check the materials to know what is the correct answer
  if (matrix_id == "C") {
    return _mt_not_defined;
  }

  if (matrix_id == "K") {
    auto matrix_type = _unsymmetric;

    for (auto & material : materials) {
      matrix_type = std::max(matrix_type, material->getMatrixType(matrix_id));
    }
  }
  return _symmetric;
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleMatrix(const ID & matrix_id) {
  if (matrix_id == "K") {
    this->assembleStiffnessMatrix();
  } else if (matrix_id == "M") {
    this->assembleMass();
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleLumpedMatrix(const ID & matrix_id) {
  if (matrix_id == "M") {
    this->assembleMassLumped();
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::beforeSolveStep() {
  for (auto & material : materials) {
    material->beforeSolveStep();
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::afterSolveStep(bool converged) {
  for (auto & material : materials) {
    material->afterSolveStep(converged);
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::predictor() { ++displacement_release; }

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::corrector() { ++displacement_release; }

/* -------------------------------------------------------------------------- */
/**
 * This function computes the internal forces as \f$F_{int} = \int_{\Omega} N
 * \sigma d\Omega@\f$
 */
void SolidMechanicsModel::assembleInternalForces() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Assemble the internal forces");

  this->internal_force->zero();

  // compute the stresses of local elements
  AKANTU_DEBUG_INFO("Compute local stresses");
  for (auto & material : materials) {
    material->computeAllStresses(_not_ghost);
  }

  /* ------------------------------------------------------------------------ */
  /* Computation of the non local part */
  if (this->non_local_manager) {
    this->non_local_manager->computeAllNonLocalStresses();
  }

  // communicate the stresses
  AKANTU_DEBUG_INFO("Send data for residual assembly");
  this->asynchronousSynchronize(SynchronizationTag::_smm_stress);

  // assemble the forces due to local stresses
  AKANTU_DEBUG_INFO("Assemble residual for local elements");
  for (auto & material : materials) {
    material->assembleInternalForces(_not_ghost);
  }

  // finalize communications
  AKANTU_DEBUG_INFO("Wait distant stresses");
  this->waitEndSynchronize(SynchronizationTag::_smm_stress);

  // assemble the stresses due to ghost elements
  AKANTU_DEBUG_INFO("Assemble residual for ghost elements");
  for (auto & material : materials) {
    material->assembleInternalForces(_ghost);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleStiffnessMatrix() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Assemble the new stiffness matrix.");

  // Check if materials need to recompute the matrix
  bool need_to_reassemble = false;

  for (auto & material : materials) {
    need_to_reassemble |= material->hasMatrixChanged("K");
  }

  if (need_to_reassemble) {
    this->getDOFManager().getMatrix("K").zero();

    // call compute stiffness matrix on each local elements
    for (auto & material : materials) {
      material->assembleStiffnessMatrix(_not_ghost);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::updateCurrentPosition() {
  if (this->current_position_release == this->displacement_release) {
    return;
  }

  this->current_position->copy(this->mesh.getNodes());

  auto cpos_it = this->current_position->begin(Model::spatial_dimension);
  auto cpos_end = this->current_position->end(Model::spatial_dimension);
  auto disp_it = this->displacement->begin(Model::spatial_dimension);

  for (; cpos_it != cpos_end; ++cpos_it, ++disp_it) {
    *cpos_it += *disp_it;
  }

  this->current_position_release = this->displacement_release;
}

/* -------------------------------------------------------------------------- */
const Array<Real> & SolidMechanicsModel::getCurrentPosition() {
  this->updateCurrentPosition();
  return *this->current_position;
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::updateDataForNonLocalCriterion(
    ElementTypeMapReal & criterion) {
  const ID field_name = criterion.getName();
  for (auto & material : materials) {
    if (!material->isInternal<Real>(field_name, _ek_regular)) {
      continue;
    }

    for (auto ghost_type : ghost_types) {
      material->flattenInternal(field_name, criterion, ghost_type, _ek_regular);
    }
  }
}

/* -------------------------------------------------------------------------- */
/* Information                                                                */
/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getStableTimeStep() {
  AKANTU_DEBUG_IN();

  Real min_dt = getStableTimeStep(_not_ghost);

  /// reduction min over all processors
  mesh.getCommunicator().allReduce(min_dt, SynchronizerOperation::_min);

  AKANTU_DEBUG_OUT();
  return min_dt;
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getStableTimeStep(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real min_dt = std::numeric_limits<Real>::max();

  this->updateCurrentPosition();

  Element elem;
  elem.ghost_type = ghost_type;

  for (auto type :
       mesh.elementTypes(Model::spatial_dimension, ghost_type, _ek_regular)) {
    elem.type = type;
    UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);
    UInt nb_element = mesh.getNbElement(type);

    auto mat_indexes = material_index(type, ghost_type).begin();
    auto mat_loc_num = material_local_numbering(type, ghost_type).begin();

    Array<Real> X(0, nb_nodes_per_element * Model::spatial_dimension);
    FEEngine::extractNodalToElementField(mesh, *current_position, X, type,
                                         _not_ghost);

    auto X_el = X.begin(Model::spatial_dimension, nb_nodes_per_element);

    for (UInt el = 0; el < nb_element;
         ++el, ++X_el, ++mat_indexes, ++mat_loc_num) {
      elem.element = *mat_loc_num;
      Real el_h = getFEEngine().getElementInradius(*X_el, type);
      Real el_c = this->materials[*mat_indexes]->getCelerity(elem);
      Real el_dt = el_h / el_c;

      min_dt = std::min(min_dt, el_dt);
    }
  }

  AKANTU_DEBUG_OUT();
  return min_dt;
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getKineticEnergy() {
  AKANTU_DEBUG_IN();

  Real ekin = 0.;
  UInt nb_nodes = mesh.getNbNodes();

  if (this->getDOFManager().hasLumpedMatrix("M")) {
    auto m_it = this->mass->begin(Model::spatial_dimension);
    auto m_end = this->mass->end(Model::spatial_dimension);
    auto v_it = this->velocity->begin(Model::spatial_dimension);

    for (UInt n = 0; m_it != m_end; ++n, ++m_it, ++v_it) {
      const auto & v = *v_it;
      const auto & m = *m_it;

      Real mv2 = 0.;
      auto is_local_node = mesh.isLocalOrMasterNode(n);
      // bool is_not_pbc_slave_node = !isPBCSlaveNode(n);
      auto count_node = is_local_node; // && is_not_pbc_slave_node;
      if (count_node) {
        for (UInt i = 0; i < Model::spatial_dimension; ++i) {
          if (m(i) > std::numeric_limits<Real>::epsilon()) {
            mv2 += v(i) * v(i) * m(i);
          }
        }
      }

      ekin += mv2;
    }
  } else if (this->getDOFManager().hasMatrix("M")) {
    Array<Real> Mv(nb_nodes, Model::spatial_dimension);
    this->getDOFManager().assembleMatMulVectToArray("displacement", "M",
                                                    *this->velocity, Mv);

    for (auto && data : zip(arange(nb_nodes), make_view(Mv, spatial_dimension),
                            make_view(*this->velocity, spatial_dimension))) {
      ekin += std::get<2>(data).dot(std::get<1>(data)) *
              static_cast<Real>(mesh.isLocalOrMasterNode(std::get<0>(data)));
    }
  } else {
    AKANTU_ERROR("No function called to assemble the mass matrix.");
  }

  mesh.getCommunicator().allReduce(ekin, SynchronizerOperation::_sum);

  AKANTU_DEBUG_OUT();
  return ekin * .5;
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getKineticEnergy(ElementType type, UInt index) {
  AKANTU_DEBUG_IN();

  UInt nb_quadrature_points = getFEEngine().getNbIntegrationPoints(type);

  Array<Real> vel_on_quad(nb_quadrature_points, Model::spatial_dimension);
  Array<UInt> filter_element(1, 1, index);

  getFEEngine().interpolateOnIntegrationPoints(*velocity, vel_on_quad,
                                               Model::spatial_dimension, type,
                                               _not_ghost, filter_element);

  auto vit = vel_on_quad.begin(Model::spatial_dimension);
  auto vend = vel_on_quad.end(Model::spatial_dimension);

  Vector<Real> rho_v2(nb_quadrature_points);

  Real rho = materials[material_index(type)(index)]->getRho();

  for (UInt q = 0; vit != vend; ++vit, ++q) {
    rho_v2(q) = rho * vit->dot(*vit);
  }

  AKANTU_DEBUG_OUT();

  return .5 * getFEEngine().integrate(rho_v2, type, index);
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getExternalWork() {
  AKANTU_DEBUG_IN();

  auto ext_force_it = external_force->begin(Model::spatial_dimension);
  auto int_force_it = internal_force->begin(Model::spatial_dimension);
  auto boun_it = blocked_dofs->begin(Model::spatial_dimension);

  decltype(ext_force_it) incr_or_velo_it;
  if (this->method == _static) {
    incr_or_velo_it =
        this->displacement_increment->begin(Model::spatial_dimension);
  } else {
    incr_or_velo_it = this->velocity->begin(Model::spatial_dimension);
  }

  Real work = 0.;

  UInt nb_nodes = this->mesh.getNbNodes();

  for (UInt n = 0; n < nb_nodes;
       ++n, ++ext_force_it, ++int_force_it, ++boun_it, ++incr_or_velo_it) {
    const auto & int_force = *int_force_it;
    const auto & ext_force = *ext_force_it;
    const auto & boun = *boun_it;
    const auto & incr_or_velo = *incr_or_velo_it;

    bool is_local_node = this->mesh.isLocalOrMasterNode(n);
    // bool is_not_pbc_slave_node = !this->isPBCSlaveNode(n);
    bool count_node = is_local_node; // && is_not_pbc_slave_node;

    if (count_node) {
      for (UInt i = 0; i < Model::spatial_dimension; ++i) {
        if (boun(i)) {
          work -= int_force(i) * incr_or_velo(i);
        } else {
          work += ext_force(i) * incr_or_velo(i);
        }
      }
    }
  }

  mesh.getCommunicator().allReduce(work, SynchronizerOperation::_sum);

  if (this->method != _static) {
    work *= this->getTimeStep();
  }

  AKANTU_DEBUG_OUT();
  return work;
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getEnergy(const std::string & energy_id) {
  AKANTU_DEBUG_IN();

  if (energy_id == "kinetic") {
    return getKineticEnergy();
  }

  if (energy_id == "external work") {
    return getExternalWork();
  }

  Real energy = 0.;
  for (auto & material : materials) {
    energy += material->getEnergy(energy_id);
  }

  /// reduction sum over all processors
  mesh.getCommunicator().allReduce(energy, SynchronizerOperation::_sum);

  AKANTU_DEBUG_OUT();
  return energy;
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getEnergy(const std::string & energy_id,
                                    ElementType type, UInt index) {
  AKANTU_DEBUG_IN();

  if (energy_id == "kinetic") {
    return getKineticEnergy(type, index);
  }

  UInt mat_index = this->material_index(type, _not_ghost)(index);
  UInt mat_loc_num = this->material_local_numbering(type, _not_ghost)(index);
  Real energy =
      this->materials[mat_index]->getEnergy(energy_id, type, mat_loc_num);

  AKANTU_DEBUG_OUT();
  return energy;
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getEnergy(const ID & energy_id,
                                    const ID & group_id) {
  auto && group = mesh.getElementGroup(group_id);
  auto energy = 0.;
  for(auto && type : group.elementTypes()) {
    for(auto el : group.getElementsIterable(type)) {
      energy += getEnergy(energy_id, el);
    }
  }

  /// reduction sum over all processors
  mesh.getCommunicator().allReduce(energy, SynchronizerOperation::_sum);

  return energy;
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::onElementsAdded(const Array<Element> & element_list,
                                          const NewElementsEvent & event) {
  AKANTU_DEBUG_IN();

  this->material_index.initialize(mesh, _element_kind = _ek_not_defined,
                                  _with_nb_element = true,
                                  _default_value = UInt(-1));
  this->material_local_numbering.initialize(
      mesh, _element_kind = _ek_not_defined, _with_nb_element = true,
      _default_value = UInt(-1));

  ElementTypeMapArray<UInt> filter("new_element_filter", this->getID());

  for (const auto & elem : element_list) {
    if (mesh.getSpatialDimension(elem.type) != spatial_dimension) {
      continue;
    }

    if (!filter.exists(elem.type, elem.ghost_type)) {
      filter.alloc(0, 1, elem.type, elem.ghost_type);
    }
    filter(elem.type, elem.ghost_type).push_back(elem.element);
  }

  // this fails in parallel if the event is sent on facet between constructor
  // and initFull \todo: to debug...
  this->assignMaterialToElements(&filter);

  for (auto & material : materials) {
    material->onElementsAdded(element_list, event);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::onElementsRemoved(
    const Array<Element> & element_list,
    const ElementTypeMapArray<UInt> & new_numbering,
    const RemovedElementsEvent & event) {
  for (auto & material : materials) {
    material->onElementsRemoved(element_list, new_numbering, event);
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::onNodesAdded(const Array<UInt> & nodes_list,
                                       const NewNodesEvent & event) {
  AKANTU_DEBUG_IN();
  UInt nb_nodes = mesh.getNbNodes();

  if (displacement) {
    displacement->resize(nb_nodes, 0.);
    ++displacement_release;
  }
  if (mass) {
    mass->resize(nb_nodes, 0.);
  }
  if (velocity) {
    velocity->resize(nb_nodes, 0.);
  }
  if (acceleration) {
    acceleration->resize(nb_nodes, 0.);
  }
  if (external_force) {
    external_force->resize(nb_nodes, 0.);
  }
  if (internal_force) {
    internal_force->resize(nb_nodes, 0.);
  }
  if (blocked_dofs) {
    blocked_dofs->resize(nb_nodes, false);
  }
  if (current_position) {
    current_position->resize(nb_nodes, 0.);
  }

  if (previous_displacement) {
    previous_displacement->resize(nb_nodes, 0.);
  }
  if (displacement_increment) {
    displacement_increment->resize(nb_nodes, 0.);
  }

  for (auto & material : materials) {
    material->onNodesAdded(nodes_list, event);
  }

  need_to_reassemble_lumped_mass = true;
  need_to_reassemble_mass = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::onNodesRemoved(const Array<UInt> & /*element_list*/,
                                         const Array<UInt> & new_numbering,
                                         const RemovedNodesEvent & /*event*/) {
  if (displacement) {
    mesh.removeNodesFromArray(*displacement, new_numbering);
    ++displacement_release;
  }
  if (mass) {
    mesh.removeNodesFromArray(*mass, new_numbering);
  }
  if (velocity) {
    mesh.removeNodesFromArray(*velocity, new_numbering);
  }
  if (acceleration) {
    mesh.removeNodesFromArray(*acceleration, new_numbering);
  }
  if (internal_force) {
    mesh.removeNodesFromArray(*internal_force, new_numbering);
  }
  if (external_force) {
    mesh.removeNodesFromArray(*external_force, new_numbering);
  }
  if (blocked_dofs) {
    mesh.removeNodesFromArray(*blocked_dofs, new_numbering);
  }

  // if (increment_acceleration)
  //   mesh.removeNodesFromArray(*increment_acceleration, new_numbering);
  if (displacement_increment) {
    mesh.removeNodesFromArray(*displacement_increment, new_numbering);
  }

  if (previous_displacement) {
    mesh.removeNodesFromArray(*previous_displacement, new_numbering);
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::printself(std::ostream & stream, int indent) const {
  std::string space(indent, AKANTU_INDENT);

  stream << space << "Solid Mechanics Model [" << std::endl;
  stream << space << " + id                : " << id << std::endl;
  stream << space << " + spatial dimension : " << Model::spatial_dimension
         << std::endl;

  stream << space << " + fem [" << std::endl;
  getFEEngine().printself(stream, indent + 2);
  stream << space << " ]" << std::endl;

  stream << space << " + nodals information [" << std::endl;
  displacement->printself(stream, indent + 2);
  if (velocity) {
    velocity->printself(stream, indent + 2);
  }
  if (acceleration) {
    acceleration->printself(stream, indent + 2);
  }
  if (mass) {
    mass->printself(stream, indent + 2);
  }
  external_force->printself(stream, indent + 2);
  internal_force->printself(stream, indent + 2);
  blocked_dofs->printself(stream, indent + 2);
  stream << space << " ]" << std::endl;

  stream << space << " + material information [" << std::endl;
  material_index.printself(stream, indent + 2);
  stream << space << " ]" << std::endl;

  stream << space << " + materials [" << std::endl;
  for (const auto & material : materials) {
    material->printself(stream, indent + 2);
  }
  stream << space << " ]" << std::endl;

  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initializeNonLocal() {
  this->non_local_manager->synchronize(*this, SynchronizationTag::_material_id);
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::insertIntegrationPointsInNeighborhoods(
    GhostType ghost_type) {
  for (auto & mat : materials) {
    MaterialNonLocalInterface * mat_non_local;
    if ((mat_non_local =
             dynamic_cast<MaterialNonLocalInterface *>(mat.get())) == nullptr) {
      continue;
    }

    ElementTypeMapArray<Real> quadrature_points_coordinates(
        "quadrature_points_coordinates_tmp_nl", this->id);
    quadrature_points_coordinates.initialize(this->getFEEngine(),
                                             _nb_component = spatial_dimension,
                                             _ghost_type = ghost_type);

    for (const auto & type : quadrature_points_coordinates.elementTypes(
             Model::spatial_dimension, ghost_type)) {
      this->getFEEngine().computeIntegrationPointsCoordinates(
          quadrature_points_coordinates(type, ghost_type), type, ghost_type);
    }

    mat_non_local->initMaterialNonLocal();

    mat_non_local->insertIntegrationPointsInNeighborhoods(
        ghost_type, quadrature_points_coordinates);
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::computeNonLocalStresses(GhostType ghost_type) {
  for (auto & mat : materials) {
    if (not aka::is_of_type<MaterialNonLocalInterface>(*mat)) {
      continue;
    }

    auto & mat_non_local = dynamic_cast<MaterialNonLocalInterface &>(*mat);
    mat_non_local.computeNonLocalStresses(ghost_type);
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::updateLocalInternal(
    ElementTypeMapReal & internal_flat, GhostType ghost_type,
    ElementKind kind) {
  const ID field_name = internal_flat.getName();
  for (auto & material : materials) {
    if (material->isInternal<Real>(field_name, kind)) {
      material->flattenInternal(field_name, internal_flat, ghost_type, kind);
    }
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::updateNonLocalInternal(
    ElementTypeMapReal & internal_flat, GhostType ghost_type,
    ElementKind kind) {

  const ID field_name = internal_flat.getName();

  for (auto & mat : materials) {
    if (not aka::is_of_type<MaterialNonLocalInterface>(*mat)) {
      continue;
    }

    auto & mat_non_local = dynamic_cast<MaterialNonLocalInterface &>(*mat);
    mat_non_local.updateNonLocalInternals(internal_flat, field_name, ghost_type,
                                          kind);
  }
}

/* -------------------------------------------------------------------------- */
FEEngine & SolidMechanicsModel::getFEEngineBoundary(const ID & name) {
  return getFEEngineClassBoundary<MyFEEngineType>(name);
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::splitElementByMaterial(
    const Array<Element> & elements,
    std::vector<Array<Element>> & elements_per_mat) const {
  for (const auto & el : elements) {
    Element mat_el = el;
    mat_el.element = this->material_local_numbering(el);
    elements_per_mat[this->material_index(el)].push_back(mat_el);
  }
}

/* -------------------------------------------------------------------------- */
UInt SolidMechanicsModel::getNbData(const Array<Element> & elements,
                                    const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  UInt nb_nodes_per_element = 0;

  for (const Element & el : elements) {
    nb_nodes_per_element += Mesh::getNbNodesPerElement(el.type);
  }

  switch (tag) {
  case SynchronizationTag::_material_id: {
    size += elements.size() * sizeof(UInt);
    break;
  }
  case SynchronizationTag::_smm_mass: {
    size += nb_nodes_per_element * sizeof(Real) *
            Model::spatial_dimension; // mass vector
    break;
  }
  case SynchronizationTag::_smm_for_gradu: {
    size += nb_nodes_per_element * Model::spatial_dimension *
            sizeof(Real); // displacement
    break;
  }
  case SynchronizationTag::_smm_boundary: {
    // force, displacement, boundary
    size += nb_nodes_per_element * Model::spatial_dimension *
            (2 * sizeof(Real) + sizeof(bool));
    break;
  }
  case SynchronizationTag::_for_dump: {
    // displacement, velocity, acceleration, residual, force
    size += nb_nodes_per_element * Model::spatial_dimension * sizeof(Real) * 5;
    break;
  }
  default: {
  }
  }

  if (tag != SynchronizationTag::_material_id) {
    splitByMaterial(elements, [&](auto && mat, auto && elements) {
      size += mat.getNbData(elements, tag);
    });
  }

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::packData(CommunicationBuffer & buffer,
                                   const Array<Element> & elements,
                                   const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  switch (tag) {
  case SynchronizationTag::_material_id: {
    packElementalDataHelper(
        material_index, buffer, elements, false, getFEEngine());
    break;
  }
  case SynchronizationTag::_smm_mass: {
    packNodalDataHelper(*mass, buffer, elements, mesh);
    break;
  }
  case SynchronizationTag::_smm_for_gradu: {
    packNodalDataHelper(*displacement, buffer, elements, mesh);
    break;
  }
  case SynchronizationTag::_for_dump: {
    packNodalDataHelper(*displacement, buffer, elements, mesh);
    packNodalDataHelper(*velocity, buffer, elements, mesh);
    packNodalDataHelper(*acceleration, buffer, elements, mesh);
    packNodalDataHelper(*internal_force, buffer, elements, mesh);
    packNodalDataHelper(*external_force, buffer, elements, mesh);
    break;
  }
  case SynchronizationTag::_smm_boundary: {
    packNodalDataHelper(*external_force, buffer, elements, mesh);
    packNodalDataHelper(*velocity, buffer, elements, mesh);
    packNodalDataHelper(*blocked_dofs, buffer, elements, mesh);
    break;
  }
  default: {
  }
  }

  if (tag != SynchronizationTag::_material_id) {
    splitByMaterial(elements, [&](auto && mat, auto && elements) {
      mat.packData(buffer, elements, tag);
    });
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::unpackData(CommunicationBuffer & buffer,
                                     const Array<Element> & elements,
                                     const SynchronizationTag & tag) {
  AKANTU_DEBUG_IN();

  switch (tag) {
  case SynchronizationTag::_material_id: {
    for (auto && element : elements) {
      UInt recv_mat_index;
      buffer >> recv_mat_index;
      UInt & mat_index = material_index(element);
      if (mat_index != UInt(-1)) {
        continue;
      }

      // add ghosts element to the correct material
      mat_index = recv_mat_index;
      UInt index = materials[mat_index]->addElement(element);
      material_local_numbering(element) = index;
    }
    break;
  }
  case SynchronizationTag::_smm_mass: {
    unpackNodalDataHelper(*mass, buffer, elements, mesh);
    break;
  }
  case SynchronizationTag::_smm_for_gradu: {
    unpackNodalDataHelper(*displacement, buffer, elements, mesh);
    break;
  }
  case SynchronizationTag::_for_dump: {
    unpackNodalDataHelper(*displacement, buffer, elements, mesh);
    unpackNodalDataHelper(*velocity, buffer, elements, mesh);
    unpackNodalDataHelper(*acceleration, buffer, elements, mesh);
    unpackNodalDataHelper(*internal_force, buffer, elements, mesh);
    unpackNodalDataHelper(*external_force, buffer, elements, mesh);
    break;
  }
  case SynchronizationTag::_smm_boundary: {
    unpackNodalDataHelper(*external_force, buffer, elements, mesh);
    unpackNodalDataHelper(*velocity, buffer, elements, mesh);
    unpackNodalDataHelper(*blocked_dofs, buffer, elements, mesh);
    break;
  }
  default: {
  }
  }

  if (tag != SynchronizationTag::_material_id) {
    splitByMaterial(elements, [&](auto && mat, auto && elements) {
      mat.unpackData(buffer, elements, tag);
    });
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
UInt SolidMechanicsModel::getNbData(const Array<UInt> & dofs,
                                    const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  //  UInt nb_nodes = mesh.getNbNodes();

  switch (tag) {
  case SynchronizationTag::_smm_uv: {
    size += sizeof(Real) * Model::spatial_dimension * 2;
    break;
  }
  case SynchronizationTag::_smm_res: /* FALLTHRU */
  case SynchronizationTag::_smm_mass: {
    size += sizeof(Real) * Model::spatial_dimension;
    break;
  }
  case SynchronizationTag::_for_dump: {
    size += sizeof(Real) * Model::spatial_dimension * 5;
    break;
  }
  default: {
    AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
  return size * dofs.size();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::packData(CommunicationBuffer & buffer,
                                   const Array<UInt> & dofs,
                                   const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  switch (tag) {
  case SynchronizationTag::_smm_uv: {
    packDOFDataHelper(*displacement, buffer, dofs);
    packDOFDataHelper(*velocity, buffer, dofs);
    break;
  }
  case SynchronizationTag::_smm_res: {
    packDOFDataHelper(*internal_force, buffer, dofs);
    break;
  }
  case SynchronizationTag::_smm_mass: {
    packDOFDataHelper(*mass, buffer, dofs);
    break;
  }
  case SynchronizationTag::_for_dump: {
    packDOFDataHelper(*displacement, buffer, dofs);
    packDOFDataHelper(*velocity, buffer, dofs);
    packDOFDataHelper(*acceleration, buffer, dofs);
    packDOFDataHelper(*internal_force, buffer, dofs);
    packDOFDataHelper(*external_force, buffer, dofs);
    break;
  }
  default: {
    AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::unpackData(CommunicationBuffer & buffer,
                                     const Array<UInt> & dofs,
                                     const SynchronizationTag & tag) {
  AKANTU_DEBUG_IN();

  switch (tag) {
  case SynchronizationTag::_smm_uv: {
    unpackDOFDataHelper(*displacement, buffer, dofs);
    unpackDOFDataHelper(*velocity, buffer, dofs);
    break;
  }
  case SynchronizationTag::_smm_res: {
    unpackDOFDataHelper(*internal_force, buffer, dofs);
    break;
  }
  case SynchronizationTag::_smm_mass: {
    unpackDOFDataHelper(*mass, buffer, dofs);
    break;
  }
  case SynchronizationTag::_for_dump: {
    unpackDOFDataHelper(*displacement, buffer, dofs);
    unpackDOFDataHelper(*velocity, buffer, dofs);
    unpackDOFDataHelper(*acceleration, buffer, dofs);
    unpackDOFDataHelper(*internal_force, buffer, dofs);
    unpackDOFDataHelper(*external_force, buffer, dofs);
    break;
  }
  default: {
    AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
}

} // namespace akantu
