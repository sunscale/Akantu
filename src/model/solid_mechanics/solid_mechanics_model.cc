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
 * @date last modification: Tue Jan 19 2016
 *
 * @brief  Implementation of the SolidMechanicsModel class
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
#include "solid_mechanics_model.hh"
#include "integrator_gauss.hh"
#include "shape_lagrange.hh"

#include "element_synchronizer.hh"
#include "sparse_matrix.hh"
#include "static_communicator.hh"
#include "synchronizer_registry.hh"

#include "dumpable_inline_impl.hh"
#ifdef AKANTU_USE_IOHELPER
#include "dumper_paraview.hh"
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
 */
SolidMechanicsModel::SolidMechanicsModel(Mesh & mesh, UInt dim, const ID & id,
                                         const MemoryID & memory_id)
    : Model(mesh, dim, id, memory_id), BoundaryCondition<SolidMechanicsModel>(),
      f_m2a(1.0), displacement(nullptr), previous_displacement(nullptr),
      displacement_increment(nullptr), mass(nullptr), velocity(nullptr),
      acceleration(nullptr), external_force(nullptr), internal_force(nullptr),
      blocked_dofs(nullptr), current_position(nullptr),
      material_index("material index", id, memory_id),
      material_local_numbering("material local numbering", id, memory_id),
      material_selector(new DefaultMaterialSelector(material_index)),
      is_default_material_selector(true), increment_flag(false),
      are_materials_instantiated(false) { //, pbc_synch(nullptr) {
  AKANTU_DEBUG_IN();

  this->registerFEEngineObject<MyFEEngineType>("SolidMechanicsFEEngine", mesh,
                                               Model::spatial_dimension);

#if defined(AKANTU_USE_IOHELPER)
  this->mesh.registerDumper<DumperParaview>("paraview_all", id, true);
  this->mesh.addDumpMesh(mesh, Model::spatial_dimension, _not_ghost,
                         _ek_regular);
#endif

  this->initDOFManager();

  this->registerDataAccessor(*this);

  if (this->mesh.isDistributed()) {
    auto & synchronizer = this->mesh.getElementSynchronizer();
    this->registerSynchronizer(synchronizer, _gst_material_id);
    this->registerSynchronizer(synchronizer, _gst_smm_mass);
    this->registerSynchronizer(synchronizer, _gst_smm_stress);
    this->registerSynchronizer(synchronizer, _gst_smm_boundary);
    this->registerSynchronizer(synchronizer, _gst_for_dump);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SolidMechanicsModel::~SolidMechanicsModel() {
  AKANTU_DEBUG_IN();

  if (is_default_material_selector) {
    delete material_selector;
    material_selector = nullptr;
  }

  for (auto & internal : this->registered_internals) {
    delete internal.second;
  }

  //  delete pbc_synch;

  AKANTU_DEBUG_OUT();
}

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
 * @param material_file the  file containing the materials to  use
 * @param method the analysis method wanted.  See the akantu::AnalysisMethod for
 * the different possibilities
 */
void SolidMechanicsModel::initFull(const ModelOptions & options) {
  Model::initFull(options);

  const SolidMechanicsModelOptions & smm_options =
      dynamic_cast<const SolidMechanicsModelOptions &>(options);

  this->method = smm_options.analysis_method;

  // initialize the vectors
  for (auto ghost_type : ghost_types) {
    for (auto type : mesh.elementTypes(Model::spatial_dimension, ghost_type,
                                       _ek_not_defined)) {
      UInt nb_element = mesh.getNbElement(type, ghost_type);
      this->material_index.alloc(nb_element, 1, type, ghost_type);
      this->material_local_numbering.alloc(nb_element, 1, type, ghost_type);
    }
  }

  if (!this->hasDefaultSolver())
    this->initNewSolver(this->method);

  // initialize pbc
  if (this->pbc_pair.size() != 0)
    this->initPBC();

  // initialize the materials
  if (this->parser->getLastParsedFile() != "") {
    this->instantiateMaterials();
  }

  this->initMaterials();

  this->initBC(*this, *displacement, *displacement_increment, *external_force);
}

/* -------------------------------------------------------------------------- */
TimeStepSolverType SolidMechanicsModel::getDefaultSolverType() const {
  return _tsst_dynamic_lumped;
}

/* -------------------------------------------------------------------------- */
ModelSolverOptions SolidMechanicsModel::getDefaultSolverOptions(
    const TimeStepSolverType & type) const {
  ModelSolverOptions options;

  switch (type) {
  case _tsst_dynamic_lumped: {
    options.non_linear_solver_type = _nls_lumped;
    options.integration_scheme_type["displacement"] = _ist_central_difference;
    options.solution_type["displacement"] = IntegrationScheme::_acceleration;
    break;
  }
  case _tsst_static: {
    options.non_linear_solver_type = _nls_newton_raphson;
    options.integration_scheme_type["displacement"] = _ist_pseudo_time;
    options.solution_type["displacement"] = IntegrationScheme::_not_defined;
    break;
  }
  case _tsst_dynamic: {
    if (this->method == _explicit_consistent_mass) {
      options.non_linear_solver_type = _nls_newton_raphson;
      options.integration_scheme_type["displacement"] = _ist_central_difference;
      options.solution_type["displacement"] = IntegrationScheme::_acceleration;
    } else {
      options.non_linear_solver_type = _nls_newton_raphson;
      options.integration_scheme_type["displacement"] = _ist_trapezoidal_rule_2;
      options.solution_type["displacement"] = IntegrationScheme::_displacement;
    }
    break;
  }
  }

  return options;
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initNewSolver(const AnalysisMethod & method) {
  ID solver_name;
  TimeStepSolverType tss_type;

  this->method = method;

  switch (this->method) {
  case _explicit_lumped_mass: {
    solver_name = "explicit_lumped";
    tss_type = _tsst_dynamic_lumped;
    break;
  }
  case _explicit_consistent_mass: {
    solver_name = "explicit";
    tss_type = _tsst_dynamic;
    break;
  }
  case _static: {
    solver_name = "static";
    tss_type = _tsst_static;
    break;
  }
  case _implicit_dynamic: {
    solver_name = "implicit";
    tss_type = _tsst_dynamic;
    break;
  }
  }

  if (!this->hasSolver(solver_name)) {
    ModelSolverOptions options = this->getDefaultSolverOptions(tss_type);
    this->getNewSolver(solver_name, tss_type, options.non_linear_solver_type);
    this->setIntegrationScheme(solver_name, "displacement",
                               options.integration_scheme_type["displacement"],
                               options.solution_type["displacement"]);
    this->setDefaultSolver(solver_name);
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initSolver(TimeStepSolverType time_step_solver_type,
                                     NonLinearSolverType) {
  DOFManager & dof_manager = this->getDOFManager();

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
  if (time_step_solver_type == _tsst_dynamic ||
      time_step_solver_type == _tsst_dynamic_lumped) {
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
void SolidMechanicsModel::initFEEngineBoundary() {
  FEEngine & fem_boundary = getFEEngineBoundary();
  fem_boundary.initShapeFunctions(_not_ghost);
  fem_boundary.initShapeFunctions(_ghost);

  fem_boundary.computeNormalsOnIntegrationPoints(_not_ghost);
  fem_boundary.computeNormalsOnIntegrationPoints(_ghost);
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
MatrixType SolidMechanicsModel::getMatrixType(const ID & matrix_id) {
  // \TODO check the materials to know what is the correct answer
  if (matrix_id == "C")
    return _mt_not_defined;

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
  for (auto & material : materials)
    material->beforeSolveStep();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::afterSolveStep() {
  for (auto & material : materials)
    material->afterSolveStep();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::predictor() { ++displacement_release; }

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::corrector() { ++displacement_release; }

/* -------------------------------------------------------------------------- */
/**
 * This function computes the internal forces as F_{int} = \int_{\Omega} N
 * \sigma d\Omega@f$
 */
void SolidMechanicsModel::assembleInternalForces() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Assemble the internal forces");

  this->internal_force->clear();

  // compute the stresses of local elements
  AKANTU_DEBUG_INFO("Compute local stresses");
  for (auto & material : materials) {
    material->computeAllStresses(_not_ghost);
  }

  /* ------------------------------------------------------------------------ */
  /* Computation of the non local part */
  if (this->non_local_manager)
    this->non_local_manager->computeAllNonLocalStresses();

  // communicate the stresses
  AKANTU_DEBUG_INFO("Send data for residual assembly");
  this->asynchronousSynchronize(_gst_smm_stress);

  // assemble the forces due to local stresses
  AKANTU_DEBUG_INFO("Assemble residual for local elements");
  for (auto & material : materials) {
    material->assembleInternalForces(_not_ghost);
  }

  // finalize communications
  AKANTU_DEBUG_INFO("Wait distant stresses");
  this->waitEndSynchronize(_gst_smm_stress);

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
    need_to_reassemble |= material->hasStiffnessMatrixChanged();
  }

  if (need_to_reassemble) {
    this->getDOFManager().getMatrix("K").clear();

    // call compute stiffness matrix on each local elements
    for (auto & material : materials) {
      material->assembleStiffnessMatrix(_not_ghost);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::updateCurrentPosition() {
  if (this->current_position_release == this->displacement_release)
    return;

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
    if (!material->isInternal<Real>(field_name, _ek_regular))
      continue;

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
  StaticCommunicator::getStaticCommunicator().allReduce(min_dt, _so_min);

  AKANTU_DEBUG_OUT();
  return min_dt;
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getStableTimeStep(const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  Real min_dt = std::numeric_limits<Real>::max();

  this->updateCurrentPosition();

  Element elem;
  elem.ghost_type = ghost_type;
  elem.kind = _ek_regular;

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
      const Vector<Real> & v = *v_it;
      const Vector<Real> & m = *m_it;

      Real mv2 = 0;
      bool is_local_node = mesh.isLocalOrMasterNode(n);
      // bool is_not_pbc_slave_node = !isPBCSlaveNode(n);
      bool count_node = is_local_node; // && is_not_pbc_slave_node;
      if (count_node) {
        for (UInt i = 0; i < Model::spatial_dimension; ++i) {
          if (m(i) > std::numeric_limits<Real>::epsilon())
            mv2 += v(i) * v(i) * m(i);
        }
      }

      ekin += mv2;
    }
  } else if (this->getDOFManager().hasMatrix("M")) {
    Array<Real> Mv(nb_nodes, Model::spatial_dimension);
    this->getDOFManager().getMatrix("M").matVecMul(*this->velocity, Mv);

    auto mv_it = Mv.begin(Model::spatial_dimension);
    auto mv_end = Mv.end(Model::spatial_dimension);
    auto v_it = this->velocity->begin(Model::spatial_dimension);

    for (; mv_it != mv_end; ++mv_it, ++v_it) {
      ekin += v_it->dot(*mv_it);
    }
  } else {
    AKANTU_DEBUG_ERROR("No function called to assemble the mass matrix.");
  }

  StaticCommunicator::getStaticCommunicator().allReduce(ekin, _so_sum);

  AKANTU_DEBUG_OUT();
  return ekin * .5;
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getKineticEnergy(const ElementType & type,
                                           UInt index) {
  AKANTU_DEBUG_IN();

  UInt nb_quadrature_points = getFEEngine().getNbIntegrationPoints(type);

  Array<Real> vel_on_quad(nb_quadrature_points, Model::spatial_dimension);
  Array<UInt> filter_element(1, 1, index);

  getFEEngine().interpolateOnIntegrationPoints(*velocity, vel_on_quad,
                                               Model::spatial_dimension, type,
                                               _not_ghost, filter_element);

  Array<Real>::vector_iterator vit =
      vel_on_quad.begin(Model::spatial_dimension);
  Array<Real>::vector_iterator vend = vel_on_quad.end(Model::spatial_dimension);

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

  auto incr_or_velo_it = this->velocity->begin(Model::spatial_dimension);
  if (this->method == _static) {
    incr_or_velo_it =
        this->displacement_increment->begin(Model::spatial_dimension);
  }

  auto ext_force_it = external_force->begin(Model::spatial_dimension);
  auto int_force_it = internal_force->begin(Model::spatial_dimension);
  auto boun_it = blocked_dofs->begin(Model::spatial_dimension);

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
        if (boun(i))
          work -= int_force(i) * incr_or_velo(i);
        else
          work += ext_force(i) * incr_or_velo(i);
      }
    }
  }

  StaticCommunicator::getStaticCommunicator().allReduce(work, _so_sum);

  if (this->method != _static)
    work *= this->getTimeStep();
  AKANTU_DEBUG_OUT();
  return work;
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getEnergy(const std::string & energy_id) {
  AKANTU_DEBUG_IN();

  if (energy_id == "kinetic") {
    return getKineticEnergy();
  } else if (energy_id == "external work") {
    return getExternalWork();
  }

  Real energy = 0.;
  for (auto & material : materials)
    energy += material->getEnergy(energy_id);

  /// reduction sum over all processors
  StaticCommunicator::getStaticCommunicator().allReduce(energy, _so_sum);

  AKANTU_DEBUG_OUT();
  return energy;
}
/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getEnergy(const std::string & energy_id,
                                    const ElementType & type, UInt index) {
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
void SolidMechanicsModel::onElementsAdded(const Array<Element> & element_list,
                                          const NewElementsEvent & event) {
  AKANTU_DEBUG_IN();

  for (auto ghost_type : ghost_types) {
    for (auto type : mesh.elementTypes(Model::spatial_dimension, ghost_type,
                                       _ek_not_defined)) {
      UInt nb_element = this->mesh.getNbElement(type, ghost_type);

      if (!material_index.exists(type, ghost_type)) {
        this->material_index.alloc(nb_element, 1, type, ghost_type);
        this->material_local_numbering.alloc(nb_element, 1, type, ghost_type);
      } else {
        this->material_index(type, ghost_type).resize(nb_element);
        this->material_local_numbering(type, ghost_type).resize(nb_element);
      }
    }
  }

  ElementTypeMapArray<UInt> filter("new_element_filter", this->getID(),
                                   this->getMemoryID());

  for (auto & elem : element_list) {
    if (!filter.exists(elem.type, elem.ghost_type))
      filter.alloc(0, 1, elem.type, elem.ghost_type);

    filter(elem.type, elem.ghost_type).push_back(elem.element);
  }

  this->assignMaterialToElements(&filter);

  for (auto & material : materials)
    material->onElementsAdded(element_list, event);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::onElementsRemoved(
    __attribute__((unused)) const Array<Element> & element_list,
    const ElementTypeMapArray<UInt> & new_numbering,
    const RemovedElementsEvent & event) {
  this->getFEEngine().initShapeFunctions(_not_ghost);
  this->getFEEngine().initShapeFunctions(_ghost);

  for (auto & material : materials) {
    material->onElementsRemoved(element_list, new_numbering, event);
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::onNodesAdded(const Array<UInt> & nodes_list,
                                       __attribute__((unused))
                                       const NewNodesEvent & event) {
  AKANTU_DEBUG_IN();
  UInt nb_nodes = mesh.getNbNodes();

  if (displacement) {
    displacement->resize(nb_nodes, 0.);
    ++displacement_release;
  }
  if (mass)
    mass->resize(nb_nodes, 0.);
  if (velocity)
    velocity->resize(nb_nodes, 0.);
  if (acceleration)
    acceleration->resize(nb_nodes, 0.);
  if (external_force)
    external_force->resize(nb_nodes, 0.);
  if (internal_force)
    internal_force->resize(nb_nodes, 0.);
  if (blocked_dofs)
    blocked_dofs->resize(nb_nodes, 0.);
  if (current_position)
    current_position->resize(nb_nodes, 0.);

  if (previous_displacement)
    previous_displacement->resize(nb_nodes, 0.);
  if (displacement_increment)
    displacement_increment->resize(nb_nodes, 0.);

  for (auto & material : materials) {
    material->onNodesAdded(nodes_list, event);
  }

  need_to_reassemble_lumped_mass = true;
  need_to_reassemble_mass = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::onNodesRemoved(__attribute__((unused))
                                         const Array<UInt> & element_list,
                                         const Array<UInt> & new_numbering,
                                         __attribute__((unused))
                                         const RemovedNodesEvent & event) {
  if (displacement) {
    mesh.removeNodesFromArray(*displacement, new_numbering);
    ++displacement_release;
  }
  if (mass)
    mesh.removeNodesFromArray(*mass, new_numbering);
  if (velocity)
    mesh.removeNodesFromArray(*velocity, new_numbering);
  if (acceleration)
    mesh.removeNodesFromArray(*acceleration, new_numbering);
  if (internal_force)
    mesh.removeNodesFromArray(*internal_force, new_numbering);
  if (external_force)
    mesh.removeNodesFromArray(*external_force, new_numbering);
  if (blocked_dofs)
    mesh.removeNodesFromArray(*blocked_dofs, new_numbering);

  // if (increment_acceleration)
  //   mesh.removeNodesFromArray(*increment_acceleration, new_numbering);
  if (displacement_increment)
    mesh.removeNodesFromArray(*displacement_increment, new_numbering);

  if (previous_displacement)
    mesh.removeNodesFromArray(*previous_displacement, new_numbering);

  // if (method != _explicit_lumped_mass) {
  //   this->initSolver();
  // }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::printself(std::ostream & stream, int indent) const {
  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
    ;

  stream << space << "Solid Mechanics Model [" << std::endl;
  stream << space << " + id                : " << id << std::endl;
  stream << space << " + spatial dimension : " << Model::spatial_dimension
         << std::endl;
  stream << space << " + fem [" << std::endl;
  getFEEngine().printself(stream, indent + 2);
  stream << space << AKANTU_INDENT << "]" << std::endl;
  stream << space << " + nodals information [" << std::endl;
  displacement->printself(stream, indent + 2);
  mass->printself(stream, indent + 2);
  velocity->printself(stream, indent + 2);
  acceleration->printself(stream, indent + 2);
  external_force->printself(stream, indent + 2);
  internal_force->printself(stream, indent + 2);
  blocked_dofs->printself(stream, indent + 2);
  stream << space << AKANTU_INDENT << "]" << std::endl;

  stream << space << " + material information [" << std::endl;
  material_index.printself(stream, indent + 2);
  stream << space << AKANTU_INDENT << "]" << std::endl;

  stream << space << " + materials [" << std::endl;
  for (auto & material : materials) {
    material->printself(stream, indent + 1);
  }
  stream << space << AKANTU_INDENT << "]" << std::endl;

  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::insertIntegrationPointsInNeighborhoods(
    const GhostType & ghost_type) {
  for (auto & mat : materials) {
    MaterialNonLocalInterface * mat_non_local;
    if ((mat_non_local =
             dynamic_cast<MaterialNonLocalInterface *>(mat.get())) == nullptr)
      continue;

    ElementTypeMapArray<Real> quadrature_points_coordinates(
        "quadrature_points_coordinates_tmp_nl", this->id, this->memory_id);
    quadrature_points_coordinates.initialize(this->getFEEngine(),
                                             _nb_component = spatial_dimension,
                                             _ghost_type = ghost_type);

    for (auto & type : quadrature_points_coordinates.elementTypes(
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
void SolidMechanicsModel::computeNonLocalStresses(
    const GhostType & ghost_type) {
  for (auto & mat : materials) {
    try {
      auto & mat_non_local = dynamic_cast<MaterialNonLocalInterface &>(*mat);
      mat_non_local.computeNonLocalStresses(ghost_type);
    } catch (std::bad_cast &) {
    }
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::updateLocalInternal(
    ElementTypeMapReal & internal_flat, const GhostType & ghost_type,
    const ElementKind & kind) {
  const ID field_name = internal_flat.getName();
  for (auto & material : materials) {
    if (material->isInternal<Real>(field_name, kind))
      material->flattenInternal(field_name, internal_flat, ghost_type, kind);
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::updateNonLocalInternal(
    ElementTypeMapReal & internal_flat, const GhostType & ghost_type,
    const ElementKind & kind) {

  const ID field_name = internal_flat.getName();

  for (auto & mat : materials) {
    try {
      auto & mat_non_local = dynamic_cast<MaterialNonLocalInterface &>(*mat);
      mat_non_local.updateNonLocalInternals(internal_flat, field_name,
                                            ghost_type, kind);
    } catch (std::bad_cast &) {
    }
  }
}

/* -------------------------------------------------------------------------- */
FEEngine & SolidMechanicsModel::getFEEngineBoundary(const ID & name) {
  return dynamic_cast<FEEngine &>(
      getFEEngineClassBoundary<MyFEEngineType>(name));
}

} // namespace akantu
