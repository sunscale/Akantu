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
#include "aka_math.hh"
#include "aka_common.hh"
#include "solid_mechanics_model.hh"

#include "group_manager_inline_impl.cc"
#include "element_group.hh"
#include "static_communicator.hh"

#include "dumpable_inline_impl.hh"
#ifdef AKANTU_USE_IOHELPER
#include "dumper_field.hh"
#include "dumper_paraview.hh"
#include "dumper_homogenizing_field.hh"
#include "dumper_internal_material_field.hh"
#include "dumper_elemental_field.hh"
#include "dumper_material_padders.hh"
#include "dumper_element_partition.hh"
#include "dumper_iohelper.hh"
#endif

#ifdef AKANTU_DAMAGE_NON_LOCAL
#include "non_local_manager.hh"
#endif
/* -------------------------------------------------------------------------- */
#include <cmath>
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__

const SolidMechanicsModelOptions
    default_solid_mechanics_model_options(_explicit_lumped_mass, false);

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
      time_step(NAN), f_m2a(1.0), displacement(NULL),
      previous_displacement(NULL), increment(NULL), mass(NULL), velocity(NULL),
      acceleration(NULL), increment_acceleration(NULL), external_force(NULL),
      internal_force(NULL), blocked_dofs(NULL), current_position(NULL),
      material_index("material index", id),
      material_local_numbering("material local numbering", id), materials(0),
      material_selector(new DefaultMaterialSelector(material_index)),
      is_default_material_selector(true), increment_flag(false),
      synch_parallel(NULL), are_materials_instantiated(false) {
  AKANTU_DEBUG_IN();

  this->createSynchronizerRegistry(this);

  this->registerFEEngineObject<MyFEEngineType>("SolidMechanicsFEEngine", mesh,
                                               spatial_dimension);

  this->mesh.registerEventHandler(*this);

#if defined(AKANTU_USE_IOHELPER)
  this->mesh.registerDumper<DumperParaview>("paraview_all", id, true);
  this->mesh.addDumpMesh(mesh, spatial_dimension, _not_ghost, _ek_regular);
#endif

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SolidMechanicsModel::~SolidMechanicsModel() {
  AKANTU_DEBUG_IN();

  std::vector<Material *>::iterator mat_it;
  for (mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    delete *mat_it;
  }

  materials.clear();

  delete synch_parallel;

  if (is_default_material_selector) {
    delete material_selector;
    material_selector = NULL;
  }

  flatten_internal_map::iterator fl_it = this->registered_internals.begin();
  flatten_internal_map::iterator fl_end = this->registered_internals.end();
  for (; fl_it != fl_end; ++fl_it) {
    delete fl_it->second;
  }

#ifdef AKANTU_DAMAGE_NON_LOCAL
  delete non_local_manager;
  non_local_manager = NULL;
#endif

  delete pbc_synch;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::setTimeStep(Real time_step) {
  this->time_step = time_step;

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
  this->initArrays();

  // set the initial condition to 0
  external_force->clear();
  velocity->clear();
  acceleration->clear();
  displacement->clear();

  // initialize pbc
  if (this->pbc_pair.size() != 0)
    this->initPBC();

// initialize the time integration schemes
// switch (this->method) {
// case _explicit_lumped_mass:
//   this->initExplicit();
//   break;
// case _explicit_consistent_mass:
//   this->initSolver();
//   this->initExplicit();
//   break;
// case _implicit_dynamic:
//   this->initImplicit(true);
//   break;
// case _static:
//   this->initImplicit(false);
//   this->initArraysPreviousDisplacment();
//   break;
// default:
//   AKANTU_EXCEPTION("analysis method not recognised by SolidMechanicsModel");
//   break;
// }

#ifdef AKANTU_DAMAGE_NON_LOCAL
  /// create the non-local manager object for non-local damage computations
  this->non_local_manager = new NonLocalManager(*this);
#endif

  // initialize the materials
  if (this->parser->getLastParsedFile() != "") {
    this->instantiateMaterials();
  }

  if (!smm_options.no_init_materials) {
    this->initMaterials();
  }

  if (increment_flag)
    this->initBC(*this, *displacement, *increment, *external_force);
  else
    this->initBC(*this, *displacement, *external_force);
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initParallel(MeshPartition * partition,
                                       DataAccessor * data_accessor) {
  AKANTU_DEBUG_IN();

  if (data_accessor == NULL)
    data_accessor = this;
  synch_parallel = &createParallelSynch(partition, data_accessor);

  synch_registry->registerSynchronizer(*synch_parallel, _gst_material_id);
  synch_registry->registerSynchronizer(*synch_parallel, _gst_smm_mass);
  synch_registry->registerSynchronizer(*synch_parallel, _gst_smm_stress);
  synch_registry->registerSynchronizer(*synch_parallel, _gst_smm_boundary);
  synch_registry->registerSynchronizer(*synch_parallel, _gst_for_dump);

  AKANTU_DEBUG_OUT();
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
// void SolidMechanicsModel::initExplicit(AnalysisMethod analysis_method) {
//   AKANTU_DEBUG_IN();

//   // in case of switch from implicit to explicit
//   if (!this->isExplicit())
//     method = analysis_method;

//   // if (integrator) delete integrator;
//   // integrator = new CentralDifference();

//   UInt nb_nodes = acceleration->getSize();
//   UInt nb_degree_of_freedom = acceleration->getNbComponent();

//   std::stringstream sstr;
//   sstr << id << ":increment_acceleration";
//   increment_acceleration =
//       &(alloc<Real>(sstr.str(), nb_nodes, nb_degree_of_freedom, Real()));

//   AKANTU_DEBUG_OUT();
// }

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initArraysPreviousDisplacment() {
  AKANTU_DEBUG_IN();

  this->setIncrementFlagOn();
  if (not this->previous_displacement) {
    UInt nb_nodes = this->mesh.getNbNodes();
    std::stringstream sstr_disp_t;
    sstr_disp_t << this->id << ":previous_displacement";

    this->previous_displacement = &(this->alloc<Real>(
        sstr_disp_t.str(), nb_nodes, this->spatial_dimension, 0.));
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Allocate all the needed vectors. By  default their are not necessarily set to
 * 0
 */
void SolidMechanicsModel::initArrays() {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = mesh.getNbNodes();
  std::stringstream sstr_disp;
  sstr_disp << id << ":displacement";
  std::stringstream sstr_forc;
  sstr_forc << id << ":force";
  std::stringstream sstr_ifor;
  sstr_forc << id << ":internal_force";
  std::stringstream sstr_boun;
  sstr_boun << id << ":blocked_dofs";

  /* ------------------------------------------------------------------------ */
  // for static
  this->displacement = &(alloc<Real>(sstr_disp.str(), nb_nodes,
                                     spatial_dimension, REAL_INIT_VALUE));
  this->internal_force = &(alloc<Real>(sstr_ifor.str(), nb_nodes,
                                       spatial_dimension, REAL_INIT_VALUE));
  this->external_force = &(alloc<Real>(sstr_forc.str(), nb_nodes,
                                       spatial_dimension, REAL_INIT_VALUE));
  this->blocked_dofs =
      &(alloc<bool>(sstr_boun.str(), nb_nodes, spatial_dimension, false));

  this->getDOFManager().registerDOFs("displacements", *this->displacement,
                                     _dst_nodal);
  this->getDOFManager().registerBlockedDOFs("displacements",
                                            *this->blocked_dofs);

  std::stringstream sstr_curp;
  sstr_curp << id << ":current_position";
  this->current_position =
      &(alloc<Real>(sstr_curp.str(), 0, spatial_dimension, REAL_INIT_VALUE));

  /* ------------------------------------------------------------------------ */
  // for dynamic
  std::stringstream sstr_velo;
  sstr_velo << id << ":velocity";
  std::stringstream sstr_acce;
  sstr_acce << id << ":acceleration";
  this->velocity = &(alloc<Real>(sstr_velo.str(), nb_nodes, spatial_dimension,
                                 REAL_INIT_VALUE));
  this->acceleration = &(alloc<Real>(sstr_acce.str(), nb_nodes,
                                     spatial_dimension, REAL_INIT_VALUE));

  this->getDOFManager().registerDOFsDerivative("displacements", 1,
                                               *this->velocity);
  this->getDOFManager().registerDOFsDerivative("displacements", 2,
                                               *this->acceleration);

  for (UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType)g;
    Mesh::type_iterator it =
        mesh.firstType(spatial_dimension, gt, _ek_not_defined);
    Mesh::type_iterator end =
        mesh.lastType(spatial_dimension, gt, _ek_not_defined);
    for (; it != end; ++it) {
      UInt nb_element = mesh.getNbElement(*it, gt);
      this->material_index.alloc(nb_element, 1, *it, gt);
      this->material_local_numbering.alloc(nb_element, 1, *it, gt);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Initialize the model,basically it  pre-compute the shapes, shapes derivatives
 * and jacobian
 *
 */
void SolidMechanicsModel::initModel() {
  /// \todo add  the current position  as a parameter to  initShapeFunctions for
  /// large deformation
  getFEEngine().initShapeFunctions(_not_ghost);
  getFEEngine().initShapeFunctions(_ghost);
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initPBC() {
  Model::initPBC();
  registerPBCSynchronizer();

  // as long as there are ones on the diagonal of the matrix, we can put
  // boudandary true for slaves
  std::map<UInt, UInt>::iterator it = pbc_pair.begin();
  std::map<UInt, UInt>::iterator end = pbc_pair.end();
  UInt dim = mesh.getSpatialDimension();
  while (it != end) {
    for (UInt i = 0; i < dim; ++i)
      (*blocked_dofs)((*it).first, i) = true;
    ++it;
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::registerPBCSynchronizer() {
  pbc_synch = new PBCSynchronizer(pbc_pair);
  synch_registry->registerSynchronizer(*pbc_synch, _gst_smm_uv);
  synch_registry->registerSynchronizer(*pbc_synch, _gst_smm_mass);
  synch_registry->registerSynchronizer(*pbc_synch, _gst_smm_res);
  synch_registry->registerSynchronizer(*pbc_synch, _gst_for_dump);
  //  changeLocalEquationNumberForPBC(pbc_pair, mesh.getSpatialDimension());
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleResidual() {
  AKANTU_DEBUG_IN();

  this->assembleInternalForces();

  this->getDOFManager().assembleToResidual("displacements",
                                           *this->external_force,  1);
  this->getDOFManager().assembleToResidual("displacements",
                                           *this->internal_force, -1);


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * This function computes the internal forces as F_{int} = \int_{\Omega} N
 * \sigma d\Omega@f$
 */
void SolidMechanicsModel::assembleInternalForces() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Assemble the internal forces");

  this->internal_force->clear();

  AKANTU_DEBUG_INFO("Compute local stresses");

  std::vector<Material *>::iterator mat_it;
  for (mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    Material & mat = **mat_it;
    mat.computeAllStresses(_not_ghost);
  }

#ifdef AKANTU_DAMAGE_NON_LOCAL
  /* ------------------------------------------------------------------------ */
  /* Computation of the non local part */
  this->non_local_manager->computeAllNonLocalStresses();
#endif

  /* ------------------------------------------------------------------------ */
  /* assembling the forces internal */
  // communicate the stress
  AKANTU_DEBUG_INFO("Send data for residual assembly");
  synch_registry->asynchronousSynchronize(_gst_smm_stress);

  AKANTU_DEBUG_INFO("Assemble residual for local elements");
  for (mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    Material & mat = **mat_it;
    mat.assembleInternalForces(_not_ghost);
  }

  AKANTU_DEBUG_INFO("Wait distant stresses");
  // finalize communications
  synch_registry->waitEndSynchronize(_gst_smm_stress);

  AKANTU_DEBUG_INFO("Assemble residual for ghost elements");
  for (mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    Material & mat = **mat_it;
    mat.assembleInternalForces(_ghost);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleJacobian() {
  this->assembleStiffnessMatrix();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleStiffnessMatrix() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Assemble the new stiffness matrix.");

  // call compute stiffness matrix on each local elements
  std::vector<Material *>::iterator mat_it;
  for (mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    (*mat_it)->assembleStiffnessMatrix(_not_ghost);
  }

  AKANTU_DEBUG_OUT();
}



/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::updateCurrentPosition() {
  AKANTU_DEBUG_IN();

  this->current_position->copy(this->mesh.getNodes());

  Array<Real>::vector_iterator cpos_it =
      this->current_position->begin(spatial_dimension);
  Array<Real>::vector_iterator cpos_end =
      this->current_position->end(spatial_dimension);
  Array<Real>::const_vector_iterator disp_it =
      this->displacement->begin(spatial_dimension);

  for (; cpos_it != cpos_end; ++cpos_it, ++disp_it) {
    *cpos_it += *disp_it;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
// void SolidMechanicsModel::initializeUpdateResidualData() {
//   AKANTU_DEBUG_IN();
//   UInt nb_nodes = mesh.getNbNodes();
//   internal_force->resize(nb_nodes);

//   /// copy the forces in residual for boundary conditions
//   this->getDOFManager().assembleToResidual("displacements",
//                                            *this->external_force);

//   // start synchronization
//   synch_registry->asynchronousSynchronize(_gst_smm_uv);
//   synch_registry->waitEndSynchronize(_gst_smm_uv);

//   this->updateCurrentPosition();

//   AKANTU_DEBUG_OUT();
// }

/* -------------------------------------------------------------------------- */
/* Explicit scheme                                                            */
/* -------------------------------------------------------------------------- */


/* -------------------------------------------------------------------------- */
// void SolidMechanicsModel::computeStresses() {
//   if (isExplicit()) {
//     // start synchronization
//     synch_registry->asynchronousSynchronize(_gst_smm_uv);
//     synch_registry->waitEndSynchronize(_gst_smm_uv);

//     // compute stresses on all local elements for each materials
//     std::vector<Material *>::iterator mat_it;
//     for (mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
//       Material & mat = **mat_it;
//       mat.computeAllStresses(_not_ghost);
//     }

// #ifdef AKANTU_DAMAGE_NON_LOCAL
//     /* Computation of the non local part */
//     this->non_local_manager->computeAllNonLocalStresses();
// #endif
//   } else {
//     std::vector<Material *>::iterator mat_it;
//     for (mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
//       Material & mat = **mat_it;
//       mat.computeAllStressesFromTangentModuli(_not_ghost);
//     }
//   }
// }

/* -------------------------------------------------------------------------- */
/* Implicit scheme                                                            */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
// /**
//  * Initialize the solver and create the sparse matrices needed.
//  *
//  */
// void SolidMechanicsModel::initSolver(__attribute__((unused))
//                                      SolverOptions & options) {
//   UInt nb_global_nodes = mesh.getNbGlobalNodes();

//   jacobian_matrix = &(this->getDOFManager().getNewMatrix("jacobian",
//   _symmetric));

//   //  jacobian_matrix->buildProfile(mesh, *dof_synchronizer,
//   spatial_dimension);

//   if (!isExplicit()) {
//     delete stiffness_matrix;
//     std::stringstream sstr_sti;
//     sstr_sti << id << ":stiffness_matrix";

//     stiffness_matrix = &(this->getDOFManager().getNewMatrix("stiffness",
//     _symmetric));
//   }

//   if (solver) solver->initialize(options);
// }

// /* --------------------------------------------------------------------------
// */
// void SolidMechanicsModel::initJacobianMatrix() {
//   // @todo make it more flexible: this is an ugly patch to treat the case of
//   non
//   // fix profile of the K matrix
//   delete jacobian_matrix;

//   jacobian_matrix = &(this->getDOFManager().getNewMatrix("jacobian",
//   "stiffness"));

//   std::stringstream sstr_solv; sstr_solv << id << ":solver";
//   delete solver;
//   solver = new SolverMumps(*jacobian_matrix, sstr_solv.str());
//   if(solver)
//     solver->initialize(_solver_no_options);
// }

/* -------------------------------------------------------------------------- */
/**
 * Initialize the implicit solver, either for dynamic or static cases,
 *
 * @param dynamic
 */
// void SolidMechanicsModel::initImplicit(bool dynamic) {
//   AKANTU_DEBUG_IN();

//   method = dynamic ? _implicit_dynamic : _static;

//   if (!increment)
//     setIncrementFlagOn();

//   initSolver();

//   // if(method == _implicit_dynamic) {
//   //   if(integrator) delete integrator;
//   //   integrator = new TrapezoidalRule2();
//   // }

//   AKANTU_DEBUG_OUT();
// }

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
// SparseMatrix & SolidMechanicsModel::initVelocityDampingMatrix() {
//   return this->getDOFManager().getNewMatrix("velocity_damping", "jacobian");
// }

// /* --------------------------------------------------------------------------
// */
// void SolidMechanicsModel::implicitPred() {
//   AKANTU_DEBUG_IN();

//   if(previous_displacement) previous_displacement->copy(*displacement);

//   if(method == _implicit_dynamic)
//     integrator->integrationSchemePred(time_step,
//                                    *displacement,
//                                    *velocity,
//                                    *acceleration,
//                                    *blocked_dofs);

//   AKANTU_DEBUG_OUT();
// }

// /* --------------------------------------------------------------------------
// */
// void SolidMechanicsModel::implicitCorr() {
//   AKANTU_DEBUG_IN();

//   if(method == _implicit_dynamic) {
//     integrator->integrationSchemeCorrDispl(time_step,
//                                         *displacement,
//                                         *velocity,
//                                         *acceleration,
//                                         *blocked_dofs,
//                                         *increment);
//   } else {
//     UInt nb_nodes = displacement->getSize();
//     UInt nb_degree_of_freedom = displacement->getNbComponent() * nb_nodes;

//     Real * incr_val = increment->storage();
//     Real * disp_val = displacement->storage();
//     bool * boun_val = blocked_dofs->storage();

//     for (UInt j = 0; j < nb_degree_of_freedom; ++j, ++disp_val, ++incr_val,
//     ++boun_val){
//       *incr_val *= (1. - *boun_val);
//       *disp_val += *incr_val;
//     }
//   }

//   AKANTU_DEBUG_OUT();
// }

/* -------------------------------------------------------------------------- */
// void SolidMechanicsModel::updateIncrement() {
//   AKANTU_DEBUG_IN();

//   AKANTU_DEBUG_ASSERT(
//       previous_displacement,
//       "The previous displacement has to be initialized."
//           << " Are you working with Finite or Ineslactic deformations?");

//   UInt nb_nodes = displacement->getSize();
//   UInt nb_degree_of_freedom = displacement->getNbComponent() * nb_nodes;

//   Real * incr_val = increment->storage();
//   Real * disp_val = displacement->storage();
//   Real * prev_disp_val = previous_displacement->storage();

//   for (UInt j = 0; j < nb_degree_of_freedom;
//        ++j, ++disp_val, ++incr_val, ++prev_disp_val)
//     *incr_val = (*disp_val - *prev_disp_val);

//   AKANTU_DEBUG_OUT();
// }

/* -------------------------------------------------------------------------- */
// void SolidMechanicsModel::updatePreviousDisplacement() {
//   AKANTU_DEBUG_IN();

//   AKANTU_DEBUG_ASSERT(
//       previous_displacement,
//       "The previous displacement has to be initialized."
//           << " Are you working with Finite or Ineslactic deformations?");

//   previous_displacement->copy(*displacement);

//   AKANTU_DEBUG_OUT();
// }

/* -------------------------------------------------------------------------- */
/* Information                                                                */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::synchronizeBoundaries() {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(synch_registry,
                      "Synchronizer registry was not initialized."
                          << " Did you call initParallel?");
  synch_registry->synchronize(_gst_smm_boundary);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::synchronizeResidual() {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(synch_registry,
                      "Synchronizer registry was not initialized."
                          << " Did you call initPBC?");
  synch_registry->synchronize(_gst_smm_res);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::setIncrementFlagOn() {
  AKANTU_DEBUG_IN();

  if (!increment) {
    UInt nb_nodes = mesh.getNbNodes();
    std::stringstream sstr_inc;
    sstr_inc << id << ":increment";
    increment = &(alloc<Real>(sstr_inc.str(), nb_nodes, spatial_dimension, 0.));
  }

  increment_flag = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getStableTimeStep() {
  AKANTU_DEBUG_IN();

  Real min_dt = getStableTimeStep(_not_ghost);

  /// reduction min over all processors
  StaticCommunicator::getStaticCommunicator().allReduce(&min_dt, 1, _so_min);

  AKANTU_DEBUG_OUT();
  return min_dt;
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getStableTimeStep(const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  Material ** mat_val = &(materials.at(0));
  Real min_dt = std::numeric_limits<Real>::max();

  this->updateCurrentPosition();

  Element elem;
  elem.ghost_type = ghost_type;
  elem.kind = _ek_regular;

  Mesh::type_iterator it =
      mesh.firstType(spatial_dimension, ghost_type, _ek_regular);
  Mesh::type_iterator end =
      mesh.lastType(spatial_dimension, ghost_type, _ek_regular);
  for (; it != end; ++it) {
    elem.type = *it;
    UInt nb_nodes_per_element = mesh.getNbNodesPerElement(*it);
    UInt nb_element = mesh.getNbElement(*it);

    Array<UInt>::const_scalar_iterator mat_indexes =
        material_index(*it, ghost_type).begin();
    Array<UInt>::const_scalar_iterator mat_loc_num =
        material_local_numbering(*it, ghost_type).begin();

    Array<Real> X(0, nb_nodes_per_element * spatial_dimension);
    FEEngine::extractNodalToElementField(mesh, *current_position, X, *it,
                                         _not_ghost);

    Array<Real>::matrix_iterator X_el =
        X.begin(spatial_dimension, nb_nodes_per_element);

    for (UInt el = 0; el < nb_element;
         ++el, ++X_el, ++mat_indexes, ++mat_loc_num) {
      elem.element = *mat_loc_num;
      Real el_h = getFEEngine().getElementInradius(*X_el, *it);
      Real el_c = mat_val[*mat_indexes]->getCelerity(elem);
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

  if (!mass && !mass_matrix)
    AKANTU_DEBUG_ERROR("No function called to assemble the mass matrix.");

  Real ekin = 0.;

  UInt nb_nodes = mesh.getNbNodes();

  Real * vel_val = velocity->storage();
  Real * mass_val = mass->storage();

  for (UInt n = 0; n < nb_nodes; ++n) {
    Real mv2 = 0;
    bool is_local_node = mesh.isLocalOrMasterNode(n);
    bool is_not_pbc_slave_node = !isPBCSlaveNode(n);
    bool count_node = is_local_node && is_not_pbc_slave_node;
    for (UInt i = 0; i < spatial_dimension; ++i) {
      if (count_node)
        mv2 += *vel_val * *vel_val * *mass_val;

      vel_val++;
      mass_val++;
    }
    ekin += mv2;
  }

  StaticCommunicator::getStaticCommunicator().allReduce(&ekin, 1, _so_sum);

  AKANTU_DEBUG_OUT();
  return ekin * .5;
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getKineticEnergy(const ElementType & type,
                                           UInt index) {
  AKANTU_DEBUG_IN();

  UInt nb_quadrature_points = getFEEngine().getNbIntegrationPoints(type);

  Array<Real> vel_on_quad(nb_quadrature_points, spatial_dimension);
  Array<UInt> filter_element(1, 1, index);

  getFEEngine().interpolateOnIntegrationPoints(*velocity, vel_on_quad,
                                               spatial_dimension, type,
                                               _not_ghost, filter_element);

  Array<Real>::vector_iterator vit = vel_on_quad.begin(spatial_dimension);
  Array<Real>::vector_iterator vend = vel_on_quad.end(spatial_dimension);

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

  Real * incr_or_velo = NULL;
  if (this->method == _static) {
    incr_or_velo = this->increment->storage();
  } else {
    incr_or_velo = this->velocity->storage();
  }
  Real * forc = external_force->storage();
  Real * resi = internal_force->storage();
  bool * boun = blocked_dofs->storage();

  Real work = 0.;

  UInt nb_nodes = this->mesh.getNbNodes();

  for (UInt n = 0; n < nb_nodes; ++n) {
    bool is_local_node = this->mesh.isLocalOrMasterNode(n);
    bool is_not_pbc_slave_node = !this->isPBCSlaveNode(n);
    bool count_node = is_local_node && is_not_pbc_slave_node;

    for (UInt i = 0; i < this->spatial_dimension; ++i) {
      if (count_node) {
        if (*boun)
          work -= *resi * *incr_or_velo;
        else
          work += *forc * *incr_or_velo;
      }

      ++incr_or_velo;
      ++forc;
      ++resi;
      ++boun;
    }
  }

  StaticCommunicator::getStaticCommunicator().allReduce(&work, 1, _so_sum);

  if (this->method != _static)
    work *= this->time_step;
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
  std::vector<Material *>::iterator mat_it;
  for (mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    energy += (*mat_it)->getEnergy(energy_id);
  }

  /// reduction sum over all processors
  StaticCommunicator::getStaticCommunicator().allReduce(&energy, 1, _so_sum);

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

  std::vector<Material *>::iterator mat_it;
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

  this->getFEEngine().initShapeFunctions(_not_ghost);
  this->getFEEngine().initShapeFunctions(_ghost);

  for (UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType)g;
    Mesh::type_iterator it =
        this->mesh.firstType(spatial_dimension, gt, _ek_not_defined);
    Mesh::type_iterator end =
        this->mesh.lastType(spatial_dimension, gt, _ek_not_defined);
    for (; it != end; ++it) {
      UInt nb_element = this->mesh.getNbElement(*it, gt);
      if (!material_index.exists(*it, gt)) {
        this->material_index.alloc(nb_element, 1, *it, gt);
        this->material_local_numbering.alloc(nb_element, 1, *it, gt);
      } else {
        this->material_index(*it, gt).resize(nb_element);
        this->material_local_numbering(*it, gt).resize(nb_element);
      }
    }
  }

  Array<Element>::const_iterator<Element> it = element_list.begin();
  Array<Element>::const_iterator<Element> end = element_list.end();

  ElementTypeMapArray<UInt> filter("new_element_filter", this->getID());

  for (UInt el = 0; it != end; ++it, ++el) {
    const Element & elem = *it;

    if (!filter.exists(elem.type, elem.ghost_type))
      filter.alloc(0, 1, elem.type, elem.ghost_type);

    filter(elem.type, elem.ghost_type).push_back(elem.element);
  }

  this->assignMaterialToElements(&filter);

  std::vector<Material *>::iterator mat_it;
  for (mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    (*mat_it)->onElementsAdded(element_list, event);
  }

  if (method == _explicit_lumped_mass)
    this->assembleMassLumped();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::onElementsRemoved(
    __attribute__((unused)) const Array<Element> & element_list,
    const ElementTypeMapArray<UInt> & new_numbering,
    const RemovedElementsEvent & event) {
  this->getFEEngine().initShapeFunctions(_not_ghost);
  this->getFEEngine().initShapeFunctions(_ghost);

  std::vector<Material *>::iterator mat_it;
  for (mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    (*mat_it)->onElementsRemoved(element_list, new_numbering, event);
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::onNodesAdded(const Array<UInt> & nodes_list,
                                       __attribute__((unused))
                                       const NewNodesEvent & event) {
  AKANTU_DEBUG_IN();
  UInt nb_nodes = mesh.getNbNodes();

  if (displacement)
    displacement->resize(nb_nodes);
  if (mass)
    mass->resize(nb_nodes);
  if (velocity)
    velocity->resize(nb_nodes);
  if (acceleration)
    acceleration->resize(nb_nodes);
  if (external_force)
    external_force->resize(nb_nodes);
  if (internal_force)
    internal_force->resize(nb_nodes);
  if (blocked_dofs)
    blocked_dofs->resize(nb_nodes);

  if (previous_displacement)
    previous_displacement->resize(nb_nodes);
  if (increment_acceleration)
    increment_acceleration->resize(nb_nodes);
  if (increment)
    increment->resize(nb_nodes);

  std::vector<Material *>::iterator mat_it;
  for (mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    (*mat_it)->onNodesAdded(nodes_list, event);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::onNodesRemoved(__attribute__((unused))
                                         const Array<UInt> & element_list,
                                         const Array<UInt> & new_numbering,
                                         __attribute__((unused))
                                         const RemovedNodesEvent & event) {
  if (displacement)
    mesh.removeNodesFromArray(*displacement, new_numbering);
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
  // if (increment)
  //   mesh.removeNodesFromArray(*increment, new_numbering);

  // if (method != _explicit_lumped_mass) {
  //   this->initSolver();
  // }
}

/* -------------------------------------------------------------------------- */
bool SolidMechanicsModel::isInternal(const std::string & field_name,
                                     const ElementKind & element_kind) {

  bool is_internal = false;

  /// check if at least one material contains field_id as an internal
  for (UInt m = 0; m < materials.size() && !is_internal; ++m) {
    is_internal |= materials[m]->isInternal<Real>(field_name, element_kind);
  }

  return is_internal;
}
/* -------------------------------------------------------------------------- */
ElementTypeMap<UInt>
SolidMechanicsModel::getInternalDataPerElem(const std::string & field_name,
                                            const ElementKind & element_kind) {

  if (!(this->isInternal(field_name, element_kind)))
    AKANTU_EXCEPTION("unknown internal " << field_name);

  for (UInt m = 0; m < materials.size(); ++m) {
    if (materials[m]->isInternal<Real>(field_name, element_kind))
      return materials[m]->getInternalDataPerElem<Real>(field_name,
                                                        element_kind);
  }

  return ElementTypeMap<UInt>();
}

/* -------------------------------------------------------------------------- */
ElementTypeMapArray<Real> &
SolidMechanicsModel::flattenInternal(const std::string & field_name,
                                     const ElementKind & kind,
                                     const GhostType ghost_type) {
  std::pair<std::string, ElementKind> key(field_name, kind);
  if (this->registered_internals.count(key) == 0) {
    this->registered_internals[key] =
        new ElementTypeMapArray<Real>(field_name, this->id);
  }

  ElementTypeMapArray<Real> * internal_flat = this->registered_internals[key];

  typedef ElementTypeMapArray<Real>::type_iterator iterator;
  iterator tit = internal_flat->firstType(spatial_dimension, ghost_type, kind);
  iterator end = internal_flat->lastType(spatial_dimension, ghost_type, kind);

  for (; tit != end; ++tit) {
    ElementType type = *tit;
    (*internal_flat)(type, ghost_type).clear();
  }

  for (UInt m = 0; m < materials.size(); ++m) {
    if (materials[m]->isInternal<Real>(field_name, kind))
      materials[m]->flattenInternal(field_name, *internal_flat, ghost_type,
                                    kind);
  }

  return *internal_flat;
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::flattenAllRegisteredInternals(
    const ElementKind & kind) {

  typedef std::map<std::pair<std::string, ElementKind>,
                   ElementTypeMapArray<Real> *>::iterator iterator;
  iterator it = this->registered_internals.begin();
  iterator end = this->registered_internals.end();

  while (it != end) {
    if (kind == it->first.second)
      this->flattenInternal(it->first.first, kind);
    ++it;
  }
}

/* -------------------------------------------------------------------------- */

void SolidMechanicsModel::onDump() {
  this->flattenAllRegisteredInternals(_ek_regular);
}
/* -------------------------------------------------------------------------- */

#ifdef AKANTU_USE_IOHELPER
dumper::Field * SolidMechanicsModel::createElementalField(
    const std::string & field_name, const std::string & group_name,
    bool padding_flag, const UInt & spatial_dimension,
    const ElementKind & kind) {

  dumper::Field * field = NULL;

  if (field_name == "partitions")
    field = mesh.createElementalField<UInt, dumper::ElementPartitionField>(
        mesh.getConnectivities(), group_name, spatial_dimension, kind);
  else if (field_name == "material_index")
    field = mesh.createElementalField<UInt, Vector, dumper::ElementalField>(
        material_index, group_name, spatial_dimension, kind);
  else {
    // this copy of field_name is used to compute derivated data such as
    // strain and von mises stress that are based on grad_u and stress
    std::string field_name_copy(field_name);

    if (field_name == "strain" || field_name == "Green strain" ||
        field_name == "principal strain" ||
        field_name == "principal Green strain")
      field_name_copy = "grad_u";
    else if (field_name == "Von Mises stress")
      field_name_copy = "stress";

    bool is_internal = this->isInternal(field_name_copy, kind);

    if (is_internal) {
      ElementTypeMap<UInt> nb_data_per_elem =
          this->getInternalDataPerElem(field_name_copy, kind);
      ElementTypeMapArray<Real> & internal_flat =
          this->flattenInternal(field_name_copy, kind);
      field = mesh.createElementalField<Real, dumper::InternalMaterialField>(
          internal_flat, group_name, spatial_dimension, kind, nb_data_per_elem);
      if (field_name == "strain") {
        dumper::ComputeStrain<false> * foo =
            new dumper::ComputeStrain<false>(*this);
        field = dumper::FieldComputeProxy::createFieldCompute(field, *foo);
      } else if (field_name == "Von Mises stress") {
        dumper::ComputeVonMisesStress * foo =
            new dumper::ComputeVonMisesStress(*this);
        field = dumper::FieldComputeProxy::createFieldCompute(field, *foo);
      } else if (field_name == "Green strain") {
        dumper::ComputeStrain<true> * foo =
            new dumper::ComputeStrain<true>(*this);
        field = dumper::FieldComputeProxy::createFieldCompute(field, *foo);
      } else if (field_name == "principal strain") {
        dumper::ComputePrincipalStrain<false> * foo =
            new dumper::ComputePrincipalStrain<false>(*this);
        field = dumper::FieldComputeProxy::createFieldCompute(field, *foo);
      } else if (field_name == "principal Green strain") {
        dumper::ComputePrincipalStrain<true> * foo =
            new dumper::ComputePrincipalStrain<true>(*this);
        field = dumper::FieldComputeProxy::createFieldCompute(field, *foo);
      }

      // treat the paddings
      if (padding_flag) {
        if (field_name == "stress") {
          if (spatial_dimension == 2) {
            dumper::StressPadder<2> * foo = new dumper::StressPadder<2>(*this);
            field = dumper::FieldComputeProxy::createFieldCompute(field, *foo);
          }
        } else if (field_name == "strain" || field_name == "Green strain") {
          if (spatial_dimension == 2) {
            dumper::StrainPadder<2> * foo = new dumper::StrainPadder<2>(*this);
            field = dumper::FieldComputeProxy::createFieldCompute(field, *foo);
          }
        }
      }

      // homogenize the field
      dumper::ComputeFunctorInterface * foo =
          dumper::HomogenizerProxy::createHomogenizer(*field);

      field = dumper::FieldComputeProxy::createFieldCompute(field, *foo);
    }
  }
  return field;
}

/* -------------------------------------------------------------------------- */
dumper::Field *
SolidMechanicsModel::createNodalFieldReal(const std::string & field_name,
                                          const std::string & group_name,
                                          bool padding_flag) {

  std::map<std::string, Array<Real> *> real_nodal_fields;
  real_nodal_fields["displacement"] = this->displacement;
  real_nodal_fields["mass"] = this->mass;
  real_nodal_fields["velocity"] = this->velocity;
  real_nodal_fields["acceleration"] = this->acceleration;
  real_nodal_fields["force"] = this->external_force;
  real_nodal_fields["residual"] = this->internal_force;
  real_nodal_fields["increment"] = this->increment;

  dumper::Field * field = NULL;
  if (padding_flag)
    field = this->mesh.createNodalField(real_nodal_fields[field_name],
                                        group_name, 3);
  else
    field =
        this->mesh.createNodalField(real_nodal_fields[field_name], group_name);

  return field;
}

/* -------------------------------------------------------------------------- */
dumper::Field *
SolidMechanicsModel::createNodalFieldBool(const std::string & field_name,
                                          const std::string & group_name,
                                          bool padding_flag) {

  std::map<std::string, Array<bool> *> uint_nodal_fields;
  uint_nodal_fields["blocked_dofs"] = blocked_dofs;

  dumper::Field * field = NULL;
  field = mesh.createNodalField(uint_nodal_fields[field_name], group_name);
  return field;
}
/* -------------------------------------------------------------------------- */
#else
/* -------------------------------------------------------------------------- */
dumper::Field * SolidMechanicsModel::createElementalField(
    __attribute__((unused)) const std::string & field_name,
    __attribute__((unused)) const std::string & group_name,
    __attribute__((unused)) bool padding_flag,
    __attribute__((unused)) const UInt & spatial_dimension,
    __attribute__((unused)) const ElementKind & kind) {
  return NULL;
}
/* -------------------------------------------------------------------------- */
dumper::Field * SolidMechanicsModel::createNodalFieldReal(
    __attribute__((unused)) const std::string & field_name,
    __attribute__((unused)) const std::string & group_name,
    __attribute__((unused)) bool padding_flag) {
  return NULL;
}

/* -------------------------------------------------------------------------- */
dumper::Field * SolidMechanicsModel::createNodalFieldBool(
    __attribute__((unused)) const std::string & field_name,
    __attribute__((unused)) const std::string & group_name,
    __attribute__((unused)) bool padding_flag) {
  return NULL;
}

#endif
/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::dump(const std::string & dumper_name) {
  this->onDump();
  EventManager::sendEvent(SolidMechanicsModelEvent::BeforeDumpEvent());
  mesh.dump(dumper_name);
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::dump(const std::string & dumper_name, UInt step) {
  this->onDump();
  EventManager::sendEvent(SolidMechanicsModelEvent::BeforeDumpEvent());
  mesh.dump(dumper_name, step);
}

/* ------------------------------------------------------------------------- */
void SolidMechanicsModel::dump(const std::string & dumper_name, Real time,
                               UInt step) {
  this->onDump();
  EventManager::sendEvent(SolidMechanicsModelEvent::BeforeDumpEvent());
  mesh.dump(dumper_name, time, step);
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::dump() {
  this->onDump();
  EventManager::sendEvent(SolidMechanicsModelEvent::BeforeDumpEvent());
  mesh.dump();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::dump(UInt step) {
  this->onDump();
  EventManager::sendEvent(SolidMechanicsModelEvent::BeforeDumpEvent());
  mesh.dump(step);
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::dump(Real time, UInt step) {
  this->onDump();
  EventManager::sendEvent(SolidMechanicsModelEvent::BeforeDumpEvent());
  mesh.dump(time, step);
}

/* -------------------------------------------------------------------------- */
// void SolidMechanicsModel::computeCauchyStresses() {
//   AKANTU_DEBUG_IN();

//   // call compute stiffness matrix on each local elements
//   std::vector<Material *>::iterator mat_it;
//   for (mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
//     Material & mat = **mat_it;
//     if (mat.isFiniteDeformation())
//       mat.computeAllCauchyStresses(_not_ghost);
//   }

//   AKANTU_DEBUG_OUT();
// }

/* -------------------------------------------------------------------------- */
// void SolidMechanicsModel::saveStressAndStrainBeforeDamage() {
//   EventManager::sendEvent(
//       SolidMechanicsModelEvent::BeginningOfDamageIterationEvent());
// }

/* -------------------------------------------------------------------------- */
// void SolidMechanicsModel::updateEnergiesAfterDamage() {
//   EventManager::sendEvent(SolidMechanicsModelEvent::AfterDamageEvent());
// }

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::printself(std::ostream & stream, int indent) const {
  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
    ;

  stream << space << "Solid Mechanics Model [" << std::endl;
  stream << space << " + id                : " << id << std::endl;
  stream << space << " + spatial dimension : " << spatial_dimension
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
  std::vector<Material *>::const_iterator mat_it;
  for (mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    const Material & mat = *(*mat_it);
    mat.printself(stream, indent + 1);
  }
  stream << space << AKANTU_INDENT << "]" << std::endl;

  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__
