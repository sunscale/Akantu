/**
 * @file   solid_mechanics_model.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Jul 27 2010
 * @date last modification: Fri Sep 19 2014
 *
 * @brief  Implementation of the SolidMechanicsModel class
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
#include "dumpable_inline_impl.hh"
#include "integration_scheme_2nd_order.hh"
#include "element_group.hh"

#include "static_communicator.hh"

#include "dof_synchronizer.hh"
#include "element_group.hh"

#include <cmath>

#ifdef AKANTU_USE_MUMPS
#include "solver_mumps.hh"
#endif

#ifdef AKANTU_USE_IOHELPER
#  include "dumper_field.hh"
#  include "dumper_paraview.hh"
#  include "dumper_homogenizing_field.hh"
#  include "dumper_material_internal_field.hh"
#  include "dumper_elemental_field.hh"
#  include "dumper_material_padders.hh"
#  include "dumper_element_partition.hh"
#  include "dumper_iohelper.hh"
#endif

/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__

const SolidMechanicsModelOptions default_solid_mechanics_model_options(_explicit_lumped_mass, false);

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
SolidMechanicsModel::SolidMechanicsModel(Mesh & mesh,
                                         UInt dim,
                                         const ID & id,
                                         const MemoryID & memory_id) :
  Model(mesh, dim, id, memory_id), 
  BoundaryCondition<SolidMechanicsModel>(),
  time_step(NAN), f_m2a(1.0),
  mass_matrix(NULL),
  velocity_damping_matrix(NULL),
  stiffness_matrix(NULL),
  jacobian_matrix(NULL),
  element_index_by_material("element index by material", id),
  material_selector(new DefaultMaterialSelector(element_index_by_material)),
  is_default_material_selector(true),
  integrator(NULL),
  increment_flag(false), solver(NULL),
  synch_parallel(NULL),
  are_materials_instantiated(false) {

  AKANTU_DEBUG_IN();

  createSynchronizerRegistry(this);

  registerFEEngineObject<MyFEEngineType>("SolidMechanicsFEEngine", mesh, spatial_dimension);

  this->displacement = NULL;
  this->mass         = NULL;
  this->velocity     = NULL;
  this->acceleration = NULL;
  this->force        = NULL;
  this->residual     = NULL;
  this->blocked_dofs = NULL;

  this->increment    = NULL;
  this->increment_acceleration = NULL;

  this->dof_synchronizer = NULL;

  this->previous_displacement = NULL;

  materials.clear();

  mesh.registerEventHandler(*this);

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
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    delete *mat_it;
  }

  materials.clear();

  delete integrator;

  delete solver;
  delete mass_matrix;
  delete velocity_damping_matrix;

  if(stiffness_matrix && stiffness_matrix != jacobian_matrix)
    delete stiffness_matrix;

  delete jacobian_matrix;

  delete synch_parallel;

  if(is_default_material_selector) {
    delete material_selector;
    material_selector = NULL;
  }

  AKANTU_DEBUG_OUT();
}

void SolidMechanicsModel::setTimeStep(Real time_step) {
  this->time_step = time_step;

#if defined(AKANTU_USE_IOHELPER)
  this->mesh.getDumper().setTimeStep(time_step);
#endif
}

/* -------------------------------------------------------------------------- */
/* Initialisation                                                             */
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

  method = smm_options.analysis_method;

  // initialize the vectors
  initArrays();

  // set the initial condition to 0
  force->clear();
  velocity->clear();
  acceleration->clear();
  displacement->clear();

  // initialize pcb
  if(pbc_pair.size()!=0)
    initPBC();

  // initialize the time integration schemes
  switch(method) {
  case _explicit_lumped_mass:
    initExplicit();
    break;
  case _explicit_consistent_mass:
    initSolver();
    initExplicit();
    break;
  case _implicit_dynamic:
    initImplicit(true);
    break;
  case _static:
    initImplicit(false);
    break;
  default:
    AKANTU_EXCEPTION("analysis method not recognised by SolidMechanicsModel");
    break;
  }

  // initialize the materials
  if(this->parser->getLastParsedFile() != "") {
    instantiateMaterials();
  }

  if(!smm_options.no_init_materials) {
    initMaterials();
  }

  if(increment_flag)
    initBC(*this, *displacement, *increment, *force);
  else
    initBC(*this, *displacement, *force);
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initParallel(MeshPartition * partition,
                                       DataAccessor * data_accessor) {
  AKANTU_DEBUG_IN();

  if (data_accessor == NULL) data_accessor = this;
  synch_parallel = &createParallelSynch(partition,data_accessor);

  synch_registry->registerSynchronizer(*synch_parallel, _gst_material_id);
  synch_registry->registerSynchronizer(*synch_parallel, _gst_smm_mass);
  synch_registry->registerSynchronizer(*synch_parallel, _gst_smm_stress);
  synch_registry->registerSynchronizer(*synch_parallel, _gst_smm_boundary);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initFEEngineBoundary() {
  FEEngine & fem_boundary = getFEEngineBoundary();
  fem_boundary.initShapeFunctions(_not_ghost);
  fem_boundary.initShapeFunctions(_ghost);

  fem_boundary.computeNormalsOnControlPoints(_not_ghost);
  fem_boundary.computeNormalsOnControlPoints(_ghost);
}


/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initExplicit(AnalysisMethod analysis_method) {
  AKANTU_DEBUG_IN();

  //in case of switch from implicit to explicit
  if(!this->isExplicit())
    method = analysis_method;

  if (integrator) delete integrator;
  integrator = new CentralDifference();

  UInt nb_nodes = acceleration->getSize();
  UInt nb_degree_of_freedom = acceleration->getNbComponent();

  std::stringstream sstr; sstr << id << ":increment_acceleration";
  increment_acceleration = &(alloc<Real>(sstr.str(), nb_nodes, nb_degree_of_freedom, Real()));

  AKANTU_DEBUG_OUT();
}

void SolidMechanicsModel::initArraysPreviousDisplacment() {
  AKANTU_DEBUG_IN();

  SolidMechanicsModel::setIncrementFlagOn();
  UInt nb_nodes = mesh.getNbNodes();
  std::stringstream sstr_disp_t;
  sstr_disp_t << id << ":previous_displacement";
  previous_displacement = &(alloc<Real > (sstr_disp_t.str(), nb_nodes, spatial_dimension, 0.));

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Allocate all the needed vectors. By  default their are not necessarily set to
 * 0
 *
 */
void SolidMechanicsModel::initArrays() {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = mesh.getNbNodes();
  std::stringstream sstr_disp; sstr_disp << id << ":displacement";
  //  std::stringstream sstr_mass; sstr_mass << id << ":mass";
  std::stringstream sstr_velo; sstr_velo << id << ":velocity";
  std::stringstream sstr_acce; sstr_acce << id << ":acceleration";
  std::stringstream sstr_forc; sstr_forc << id << ":force";
  std::stringstream sstr_resi; sstr_resi << id << ":residual";
  std::stringstream sstr_boun; sstr_boun << id << ":blocked_dofs";

  displacement = &(alloc<Real>(sstr_disp.str(), nb_nodes, spatial_dimension, REAL_INIT_VALUE));
  //  mass         = &(alloc<Real>(sstr_mass.str(), nb_nodes, spatial_dimension, 0));
  velocity     = &(alloc<Real>(sstr_velo.str(), nb_nodes, spatial_dimension, REAL_INIT_VALUE));
  acceleration = &(alloc<Real>(sstr_acce.str(), nb_nodes, spatial_dimension, REAL_INIT_VALUE));
  force        = &(alloc<Real>(sstr_forc.str(), nb_nodes, spatial_dimension, REAL_INIT_VALUE));
  residual     = &(alloc<Real>(sstr_resi.str(), nb_nodes, spatial_dimension, REAL_INIT_VALUE));
  blocked_dofs = &(alloc<bool>(sstr_boun.str(), nb_nodes, spatial_dimension, false));

  std::stringstream sstr_curp; sstr_curp << id << ":current_position";
  current_position = &(alloc<Real>(sstr_curp.str(), 0, spatial_dimension, REAL_INIT_VALUE));

  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;
    Mesh::type_iterator it  = mesh.firstType(spatial_dimension, gt, _ek_not_defined);
    Mesh::type_iterator end = mesh.lastType(spatial_dimension, gt, _ek_not_defined);
    for(; it != end; ++it) {
      UInt nb_element = mesh.getNbElement(*it, gt);
      element_index_by_material.alloc(nb_element, 2, *it, gt);
    }
  }

  dof_synchronizer = new DOFSynchronizer(mesh, spatial_dimension);
  dof_synchronizer->initLocalDOFEquationNumbers();
  dof_synchronizer->initGlobalDOFEquationNumbers();

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

  // as long as there are ones on the diagonal of the matrix, we can put boudandary true for slaves
  std::map<UInt, UInt>::iterator it = pbc_pair.begin();
  std::map<UInt, UInt>::iterator end = pbc_pair.end();
  UInt dim = mesh.getSpatialDimension();
  while(it != end) {
    for (UInt i=0; i<dim; ++i)
      (*blocked_dofs)((*it).first,i) = true;
    ++it;
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::registerPBCSynchronizer(){
  PBCSynchronizer * synch = new PBCSynchronizer(pbc_pair);
  synch_registry->registerSynchronizer(*synch, _gst_smm_uv);
  synch_registry->registerSynchronizer(*synch, _gst_smm_mass);
  synch_registry->registerSynchronizer(*synch, _gst_smm_res);
  synch_registry->registerSynchronizer(*synch, _gst_for_dump);
  changeLocalEquationNumberForPBC(pbc_pair, mesh.getSpatialDimension());
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::updateCurrentPosition() {
  AKANTU_DEBUG_IN();
  UInt nb_nodes = mesh.getNbNodes();

  current_position->resize(nb_nodes);
  Real * current_position_val = current_position->storage();
  Real * position_val         = mesh.getNodes().storage();
  Real * displacement_val     = displacement->storage();

  /// compute current_position = initial_position + displacement
  memcpy(current_position_val, position_val, nb_nodes*spatial_dimension*sizeof(Real));
  for (UInt n = 0; n < nb_nodes*spatial_dimension; ++n) {
    *current_position_val++ += *displacement_val++;
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initializeUpdateResidualData() {
  AKANTU_DEBUG_IN();
  UInt nb_nodes = mesh.getNbNodes();
  residual->resize(nb_nodes);

  /// copy the forces in residual for boundary conditions
  memcpy(residual->storage(), force->storage(), nb_nodes*spatial_dimension*sizeof(Real));

  // start synchronization
  synch_registry->asynchronousSynchronize(_gst_smm_uv);
  synch_registry->waitEndSynchronize(_gst_smm_uv);

  updateCurrentPosition();

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
/* Explicit scheme                                                            */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/**
 * This function compute  the second member of the motion  equation.  That is to
 * say the sum of forces @f$ r = F_{ext} - F_{int} @f$. @f$ F_{ext} @f$ is given
 * by the  user in  the force  vector, and @f$  F_{int} @f$  is computed  as @f$
 * F_{int} = \int_{\Omega} N \sigma d\Omega@f$
 *
 */
void SolidMechanicsModel::updateResidual(bool need_initialize) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Assemble the internal forces");
  // f = f_ext - f_int

  // f = f_ext
  if(need_initialize) initializeUpdateResidualData();

  AKANTU_DEBUG_INFO("Compute local stresses");

  std::vector<Material *>::iterator mat_it;
  for (mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    Material & mat = **mat_it;
    mat.computeAllStresses(_not_ghost);
  }

#ifdef AKANTU_DAMAGE_NON_LOCAL
  /* ------------------------------------------------------------------------ */
  /* Computation of the non local part */
  synch_registry->asynchronousSynchronize(_gst_mnl_for_average);
  AKANTU_DEBUG_INFO("Compute non local stresses for local elements");

  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    Material & mat = **mat_it;
    mat.computeAllNonLocalStresses(_not_ghost);
  }


  AKANTU_DEBUG_INFO("Wait distant non local stresses");
  synch_registry->waitEndSynchronize(_gst_mnl_for_average);

  AKANTU_DEBUG_INFO("Compute non local stresses for ghosts elements");
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    Material & mat = **mat_it;
    mat.computeAllNonLocalStresses(_ghost);
  }
#endif

  /* ------------------------------------------------------------------------ */
  /* assembling the forces internal */
  // communicate the stress
  AKANTU_DEBUG_INFO("Send data for residual assembly");
  synch_registry->asynchronousSynchronize(_gst_smm_stress);

  AKANTU_DEBUG_INFO("Assemble residual for local elements");
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    Material & mat = **mat_it;
    mat.assembleResidual(_not_ghost);
  }


  AKANTU_DEBUG_INFO("Wait distant stresses");
  // finalize communications
  synch_registry->waitEndSynchronize(_gst_smm_stress);

  AKANTU_DEBUG_INFO("Assemble residual for ghost elements");
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    Material & mat = **mat_it;
    mat.assembleResidual(_ghost);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::computeStresses() {
  if (isExplicit()) {
    // start synchronization
    synch_registry->asynchronousSynchronize(_gst_smm_uv);
    synch_registry->waitEndSynchronize(_gst_smm_uv);

    // compute stresses on all local elements for each materials
    std::vector<Material *>::iterator mat_it;
    for (mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
      Material & mat = **mat_it;
      mat.computeAllStresses(_not_ghost);
    }

    /* ------------------------------------------------------------------------ */
#ifdef AKANTU_DAMAGE_NON_LOCAL
    /* Computation of the non local part */
    synch_registry->asynchronousSynchronize(_gst_mnl_for_average);

    for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
      Material & mat = **mat_it;
      mat.computeAllNonLocalStresses(_not_ghost);
    }

    synch_registry->waitEndSynchronize(_gst_mnl_for_average);

    for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
      Material & mat = **mat_it;
      mat.computeAllNonLocalStresses(_ghost);
    }
#endif
  } else {
    std::vector<Material *>::iterator mat_it;
    for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
      Material & mat = **mat_it;
      mat.computeAllStressesFromTangentModuli(_not_ghost);
    }
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::updateResidualInternal() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Update the residual");
  // f = f_ext - f_int - Ma - Cv = r - Ma - Cv;

  if(method != _static) {
    // f -= Ma
    if(mass_matrix) {
      // if full mass_matrix
      Array<Real> * Ma = new Array<Real>(*acceleration, true, "Ma");
      *Ma *= *mass_matrix;
      /// \todo check unit conversion for implicit dynamics
      //      *Ma /= f_m2a
      *residual -= *Ma;
      delete Ma;
    } else if (mass) {

      // else lumped mass
      UInt nb_nodes = acceleration->getSize();
      UInt nb_degree_of_freedom = acceleration->getNbComponent();

      Real * mass_val     = mass->storage();
      Real * accel_val    = acceleration->storage();
      Real * res_val      = residual->storage();
      bool * blocked_dofs_val = blocked_dofs->storage();

      for (UInt n = 0; n < nb_nodes * nb_degree_of_freedom; ++n) {
        if(!(*blocked_dofs_val)) {
          *res_val -= *accel_val * *mass_val /f_m2a;
        }
        blocked_dofs_val++;
        res_val++;
        mass_val++;
        accel_val++;
      }
    } else {
      AKANTU_DEBUG_ERROR("No function called to assemble the mass matrix.");
    }

    // f -= Cv
    if(velocity_damping_matrix) {
      Array<Real> * Cv = new Array<Real>(*velocity);
      *Cv *= *velocity_damping_matrix;
      *residual -= *Cv;
      delete Cv;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::updateAcceleration() {
  AKANTU_DEBUG_IN();

  updateResidualInternal();

  if(method == _explicit_lumped_mass) {
    /* residual = residual_{n+1} - M * acceleration_n therefore
       solution = increment acceleration not acceleration */
    solveLumped(*increment_acceleration,
                *mass,
                *residual,
                *blocked_dofs,
                f_m2a);
  } else if (method == _explicit_consistent_mass) {
    solve<NewmarkBeta::_acceleration_corrector>(*increment_acceleration);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::solveLumped(Array<Real> & x,
                                      const Array<Real> & A,
                                      const Array<Real> & b,
                                      const Array<bool> & blocked_dofs,
                                      Real alpha) {

  Real * A_val = A.storage();
  Real * b_val = b.storage();
  Real * x_val = x.storage();
  bool * blocked_dofs_val = blocked_dofs.storage();

  UInt nb_degrees_of_freedom = x.getSize() * x.getNbComponent();

  for (UInt n = 0; n < nb_degrees_of_freedom; ++n) {
    if(!(*blocked_dofs_val)) {
      *x_val = alpha * (*b_val / *A_val);
    }
    x_val++;
    A_val++;
    b_val++;
    blocked_dofs_val++;
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::explicitPred() {
  AKANTU_DEBUG_IN();

  if(increment_flag) {
    if(previous_displacement) increment->copy(*previous_displacement);
    else increment->copy(*displacement);
  }

  AKANTU_DEBUG_ASSERT(integrator,"itegrator should have been allocated: "
                      << "have called initExplicit ? "
                      << "or initImplicit ?");

  integrator->integrationSchemePred(time_step,
                                    *displacement,
                                    *velocity,
                                    *acceleration,
                                    *blocked_dofs);

  if(increment_flag) {
    Real * inc_val = increment->storage();
    Real * dis_val = displacement->storage();
    UInt nb_degree_of_freedom = displacement->getSize() * displacement->getNbComponent();

    for (UInt n = 0; n < nb_degree_of_freedom; ++n) {
      *inc_val = *dis_val - *inc_val;
      inc_val++;
      dis_val++;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::explicitCorr() {
  AKANTU_DEBUG_IN();

  integrator->integrationSchemeCorrAccel(time_step,
                                         *displacement,
                                         *velocity,
                                         *acceleration,
                                         *blocked_dofs,
                                         *increment_acceleration);

  if(previous_displacement) previous_displacement->copy(*displacement);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::solveStep() {
  AKANTU_DEBUG_IN();

  EventManager::sendEvent(SolidMechanicsModelEvent::BeforeSolveStepEvent(method));

  this->explicitPred();
  this->updateResidual();
  this->updateAcceleration();
  this->explicitCorr();

  EventManager::sendEvent(SolidMechanicsModelEvent::AfterSolveStepEvent(method));

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
/* Implicit scheme                                                            */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/**
 * Initialize the solver and create the sparse matrices needed.
 *
 */
void SolidMechanicsModel::initSolver(__attribute__((unused)) SolverOptions & options) {
#if !defined(AKANTU_USE_MUMPS) // or other solver in the future \todo add AKANTU_HAS_SOLVER in CMake
  AKANTU_DEBUG_ERROR("You should at least activate one solver.");
#else
  UInt nb_global_nodes = mesh.getNbGlobalNodes();

  delete jacobian_matrix;
  std::stringstream sstr; sstr << id << ":jacobian_matrix";
  jacobian_matrix = new SparseMatrix(nb_global_nodes * spatial_dimension, _symmetric, sstr.str(), memory_id);

  jacobian_matrix->buildProfile(mesh, *dof_synchronizer, spatial_dimension);

  if (!isExplicit()) {
    delete stiffness_matrix;
    std::stringstream sstr_sti; sstr_sti << id << ":stiffness_matrix";
    stiffness_matrix = new SparseMatrix(*jacobian_matrix, sstr_sti.str(), memory_id);
  }

#ifdef AKANTU_USE_MUMPS
  std::stringstream sstr_solv; sstr_solv << id << ":solver";
  solver = new SolverMumps(*jacobian_matrix, sstr_solv.str());

  dof_synchronizer->initScatterGatherCommunicationScheme();
#else
  AKANTU_DEBUG_ERROR("You should at least activate one solver.");
#endif //AKANTU_USE_MUMPS

  if(solver)
    solver->initialize(options);
#endif //AKANTU_HAS_SOLVER
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initJacobianMatrix() {
#ifdef AKANTU_USE_MUMPS
  
  // @todo make it more flexible: this is an ugly patch to treat the case of non
  // fix profile of the K matrix
  delete jacobian_matrix;
  std::stringstream sstr_sti; sstr_sti << id << ":jacobian_matrix";
  jacobian_matrix = new SparseMatrix(*stiffness_matrix, sstr_sti.str(), memory_id);

  std::stringstream sstr_solv; sstr_solv << id << ":solver";
  delete solver;
  solver = new SolverMumps(*jacobian_matrix, sstr_solv.str());
  if(solver)
    solver->initialize(_solver_no_options);
#else
  AKANTU_DEBUG_ERROR("You should at least activate one solver.");
#endif
}

/* -------------------------------------------------------------------------- */
/**
 * Initialize the implicit solver, either for dynamic or static cases,
 *
 * @param dynamic
 */
void SolidMechanicsModel::initImplicit(bool dynamic, SolverOptions & solver_options) {
  AKANTU_DEBUG_IN();

  method = dynamic ? _implicit_dynamic : _static;

  if (!increment) setIncrementFlagOn();

  initSolver(solver_options);

  if(method == _implicit_dynamic) {
    if(integrator) delete integrator;
    integrator = new TrapezoidalRule2();
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initialAcceleration() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Solving Ma = f");

  Solver * acc_solver = NULL;

  std::stringstream sstr; sstr << id << ":tmp_mass_matrix";
  SparseMatrix * tmp_mass = new SparseMatrix(*mass_matrix, sstr.str(), memory_id);

#ifdef AKANTU_USE_MUMPS
  std::stringstream sstr_solver; sstr << id << ":solver_mass_matrix";
  acc_solver = new SolverMumps(*mass_matrix, sstr_solver.str());

  dof_synchronizer->initScatterGatherCommunicationScheme();
#else
  AKANTU_DEBUG_ERROR("You should at least activate one solver.");
#endif //AKANTU_USE_MUMPS

  acc_solver->initialize();

  tmp_mass->applyBoundary(*blocked_dofs);

  acc_solver->setRHS(*residual);
  acc_solver->solve(*acceleration);

  delete acc_solver;
  delete tmp_mass;

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleStiffnessMatrix() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Assemble the new stiffness matrix.");

  stiffness_matrix->clear();

  // call compute stiffness matrix on each local elements
  std::vector<Material *>::iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    (*mat_it)->assembleStiffnessMatrix(_not_ghost);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SparseMatrix & SolidMechanicsModel::initVelocityDampingMatrix() {
  if(!velocity_damping_matrix)
    velocity_damping_matrix =
      new SparseMatrix(*jacobian_matrix, id + ":velocity_damping_matrix", memory_id);

  return *velocity_damping_matrix;
}

/* -------------------------------------------------------------------------- */
template<>
bool SolidMechanicsModel::testConvergence<_scc_increment>(Real tolerance, Real & error){
  AKANTU_DEBUG_IN();

  UInt nb_nodes = displacement->getSize();
  UInt nb_degree_of_freedom = displacement->getNbComponent();

  error = 0;
  Real norm[2] = {0., 0.};
  Real * increment_val    = increment->storage();
  bool * blocked_dofs_val     = blocked_dofs->storage();
  Real * displacement_val = displacement->storage();

  for (UInt n = 0; n < nb_nodes; ++n) {
    bool is_local_node = mesh.isLocalOrMasterNode(n);
    for (UInt d = 0; d < nb_degree_of_freedom; ++d) {
      if(!(*blocked_dofs_val) && is_local_node) {
        norm[0] += *increment_val * *increment_val;
        norm[1] += *displacement_val * *displacement_val;
      }
      blocked_dofs_val++;
      increment_val++;
      displacement_val++;
    }
  }

  StaticCommunicator::getStaticCommunicator().allReduce(norm, 2, _so_sum);

  norm[0] = sqrt(norm[0]);
  norm[1] = sqrt(norm[1]);

  AKANTU_DEBUG_ASSERT(!Math::isnan(norm[0]), "Something goes wrong in the solve phase");

  if (norm[1] < Math::getTolerance()) {
    error = norm[0];
    AKANTU_DEBUG_OUT();
    //    cout<<"Error 1: "<<error<<endl;
    return error < tolerance;
  }


  AKANTU_DEBUG_OUT();
  if(norm[1] > Math::getTolerance())
    error = norm[0] / norm[1];
  else
    error = norm[0]; //In case the total displacement is zero!

  //  cout<<"Error 2: "<<error<<endl;

  return (error < tolerance);
}


/* -------------------------------------------------------------------------- */
template<>
bool SolidMechanicsModel::testConvergence<_scc_residual>(Real tolerance, Real & norm) {
  AKANTU_DEBUG_IN();



  UInt nb_nodes = residual->getSize();

  norm = 0;
  Real * residual_val = residual->storage();
  bool * blocked_dofs_val = blocked_dofs->storage();

  for (UInt n = 0; n < nb_nodes; ++n) {
    bool is_local_node = mesh.isLocalOrMasterNode(n);
    if(is_local_node) {
      for (UInt d = 0; d < spatial_dimension; ++d) {
        if(!(*blocked_dofs_val)) {
          norm += *residual_val * *residual_val;
        }
        blocked_dofs_val++;
        residual_val++;
      }
    } else {
      blocked_dofs_val += spatial_dimension;
      residual_val += spatial_dimension;
    }
  }

  StaticCommunicator::getStaticCommunicator().allReduce(&norm, 1, _so_sum);

  norm = sqrt(norm);

  AKANTU_DEBUG_ASSERT(!Math::isnan(norm), "Something goes wrong in the solve phase");

  AKANTU_DEBUG_OUT();
  return (norm < tolerance);
}

/* -------------------------------------------------------------------------- */
template<>
bool SolidMechanicsModel::testConvergence<_scc_residual_mass_wgh>(Real tolerance,
                                                                  Real & norm) {
  AKANTU_DEBUG_IN();



  UInt nb_nodes = residual->getSize();

  norm = 0;
  Real * residual_val = residual->storage();
  Real * mass_val = this->mass->storage();
  bool * blocked_dofs_val = blocked_dofs->storage();

  for (UInt n = 0; n < nb_nodes; ++n) {
    bool is_local_node = mesh.isLocalOrMasterNode(n);
    if(is_local_node) {
      for (UInt d = 0; d < spatial_dimension; ++d) {
        if(!(*blocked_dofs_val)) {
          norm += *residual_val * *residual_val/(*mass_val * *mass_val);
        }
        blocked_dofs_val++;
        residual_val++;
        mass_val++;
      }
    } else {
      blocked_dofs_val += spatial_dimension;
      residual_val += spatial_dimension;
      mass_val += spatial_dimension;
    }
  }

  StaticCommunicator::getStaticCommunicator().allReduce(&norm, 1, _so_sum);

  norm = sqrt(norm);

  AKANTU_DEBUG_ASSERT(!Math::isnan(norm), "Something goes wrong in the solve phase");

  AKANTU_DEBUG_OUT();
  return (norm < tolerance);
}

/* -------------------------------------------------------------------------- */
bool SolidMechanicsModel::testConvergenceResidual(Real tolerance){
  AKANTU_DEBUG_IN();

  Real error=0;
  bool res = this->testConvergence<_scc_residual>(tolerance, error);
  AKANTU_DEBUG_OUT();
  return res;
}

/* -------------------------------------------------------------------------- */
bool SolidMechanicsModel::testConvergenceResidual(Real tolerance, Real & error){
  AKANTU_DEBUG_IN();

  bool res = this->testConvergence<_scc_residual>(tolerance, error);

  AKANTU_DEBUG_OUT();
  return res;
}

/* -------------------------------------------------------------------------- */
bool SolidMechanicsModel::testConvergenceIncrement(Real tolerance){
  AKANTU_DEBUG_IN();

  Real error=0;
  bool res = this->testConvergence<_scc_increment>(tolerance, error);

  AKANTU_DEBUG_OUT();
  return res;
}

/* -------------------------------------------------------------------------- */
bool SolidMechanicsModel::testConvergenceIncrement(Real tolerance, Real & error){
  AKANTU_DEBUG_IN();

  bool res = this->testConvergence<_scc_increment>(tolerance, error);

  AKANTU_DEBUG_OUT();
  return res;
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::implicitPred() {
  AKANTU_DEBUG_IN();

  if(previous_displacement) previous_displacement->copy(*displacement);

  if(method == _implicit_dynamic)
    integrator->integrationSchemePred(time_step,
                                      *displacement,
                                      *velocity,
                                      *acceleration,
                                      *blocked_dofs);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::implicitCorr() {
  AKANTU_DEBUG_IN();

  if(method == _implicit_dynamic) {
    integrator->integrationSchemeCorrDispl(time_step,
                                           *displacement,
                                           *velocity,
                                           *acceleration,
                                           *blocked_dofs,
                                           *increment);
  } else {
    UInt nb_nodes = displacement->getSize();
    UInt nb_degree_of_freedom = displacement->getNbComponent() * nb_nodes;

    Real * incr_val = increment->storage();
    Real * disp_val = displacement->storage();
    bool * boun_val = blocked_dofs->storage();

    for (UInt j = 0; j < nb_degree_of_freedom; ++j, ++disp_val, ++incr_val, ++boun_val){
      *incr_val *= (1. - *boun_val);
      *disp_val += *incr_val;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::updateIncrement() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(previous_displacement,"The previous displacement has to be initialized."
                      << " Are you working with Finite or Ineslactic deformations?");

  UInt nb_nodes = displacement->getSize();
  UInt nb_degree_of_freedom = displacement->getNbComponent() * nb_nodes;

  Real * incr_val = increment->storage();
  Real * disp_val = displacement->storage();
  Real * prev_disp_val = previous_displacement->storage();

  for (UInt j = 0; j < nb_degree_of_freedom; ++j, ++disp_val, ++incr_val, ++prev_disp_val)
    *incr_val = (*disp_val - *prev_disp_val);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::updatePreviousDisplacement() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(previous_displacement,"The previous displacement has to be initialized."
                      << " Are you working with Finite or Ineslactic deformations?");

  previous_displacement->copy(*displacement);

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
/* Information                                                                */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::synchronizeBoundaries() {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(synch_registry,"Synchronizer registry was not initialized."
                      << " Did you call initParallel?");
  synch_registry->synchronize(_gst_smm_boundary);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::synchronizeResidual() {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(synch_registry,"Synchronizer registry was not initialized."
                      << " Did you call initPBC?");
  synch_registry->synchronize(_gst_smm_res);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::setIncrementFlagOn() {
  AKANTU_DEBUG_IN();

  if(!increment) {
    UInt nb_nodes = mesh.getNbNodes();
    std::stringstream sstr_inc; sstr_inc << id << ":increment";
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

  updateCurrentPosition();

  Element elem;
  elem.ghost_type = ghost_type;
  elem.kind = _ek_regular;

  Mesh::type_iterator it  = mesh.firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator end = mesh.lastType(spatial_dimension, ghost_type);
  for(; it != end; ++it) {
    elem.type = *it;
    UInt nb_nodes_per_element = mesh.getNbNodesPerElement(*it);
    UInt nb_element           = mesh.getNbElement(*it);

    Array<UInt>::iterator< Vector<UInt> > eibm =
      element_index_by_material(*it, ghost_type).begin(2);

    Array<Real> X(0, nb_nodes_per_element*spatial_dimension);
    FEEngine::extractNodalToElementField(mesh, *current_position,
                                         X, *it, _not_ghost);

    Array<Real>::matrix_iterator X_el =
      X.begin(spatial_dimension, nb_nodes_per_element);

    for (UInt el = 0; el < nb_element; ++el, ++X_el, ++eibm) {
      elem.element = (*eibm)(1);
      Real el_h  = getFEEngine().getElementInradius(*X_el, *it);
      Real el_c  = mat_val[(*eibm)(0)]->getCelerity(elem);
      Real el_dt = el_h / el_c;

      min_dt = std::min(min_dt, el_dt);
    }
  }

  AKANTU_DEBUG_OUT();
  return min_dt;
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getPotentialEnergy() {
  AKANTU_DEBUG_IN();
  Real energy = 0.;

  /// call update residual on each local elements
  std::vector<Material *>::iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    energy += (*mat_it)->getPotentialEnergy();
  }

  /// reduction sum over all processors
  StaticCommunicator::getStaticCommunicator().allReduce(&energy, 1, _so_sum);

  AKANTU_DEBUG_OUT();
  return energy;
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getKineticEnergy() {
  AKANTU_DEBUG_IN();

  if (!mass && !mass_matrix)
    AKANTU_DEBUG_ERROR("No function called to assemble the mass matrix.");


  Real ekin = 0.;

  UInt nb_nodes = mesh.getNbNodes();

  Real * vel_val  = velocity->storage();
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
Real SolidMechanicsModel::getKineticEnergy(const ElementType & type, UInt index) {
  AKANTU_DEBUG_IN();

  UInt nb_quadrature_points = getFEEngine().getNbQuadraturePoints(type);

  Array<Real> vel_on_quad(nb_quadrature_points, spatial_dimension);
  Array<UInt> filter_element(1, 1, index);

  getFEEngine().interpolateOnQuadraturePoints(*velocity, vel_on_quad,
                                              spatial_dimension,
                                              type, _not_ghost,
                                              filter_element);

  Array<Real>::vector_iterator vit   = vel_on_quad.begin(spatial_dimension);
  Array<Real>::vector_iterator vend  = vel_on_quad.end(spatial_dimension);

  Vector<Real> rho_v2(nb_quadrature_points);

  Real rho = materials[element_index_by_material(type)(index, 0)]->getRho();

  for (UInt q = 0; vit != vend; ++vit, ++q) {
    rho_v2(q) = rho * vit->dot(*vit);
  }

  AKANTU_DEBUG_OUT();

  return .5*getFEEngine().integrate(rho_v2, type, index);
}


/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getExternalWork() {
  AKANTU_DEBUG_IN();

  Real * velo = velocity->storage();
  Real * forc = force->storage();
  Real * resi = residual->storage();
  bool * boun = blocked_dofs->storage();

  Real work = 0.;

  UInt nb_nodes = mesh.getNbNodes();

  for (UInt n = 0; n < nb_nodes; ++n) {
    bool is_local_node = mesh.isLocalOrMasterNode(n);
    bool is_not_pbc_slave_node = !isPBCSlaveNode(n);
    bool count_node = is_local_node && is_not_pbc_slave_node;

    for (UInt i = 0; i < spatial_dimension; ++i) {
      if (count_node) {
        if(*boun)
          work -= *resi * *velo * time_step;
        else
          work += *forc * *velo * time_step;
      }

      ++velo;
      ++forc;
      ++resi;
      ++boun;
    }
  }

  StaticCommunicator::getStaticCommunicator().allReduce(&work, 1, _so_sum);

  AKANTU_DEBUG_OUT();
  return work;
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getEnergy(const std::string & energy_id) {
  AKANTU_DEBUG_IN();

  if (energy_id == "kinetic") {
    return getKineticEnergy();
  } else if (energy_id == "external work"){
    return getExternalWork();
  }

  Real energy = 0.;
  std::vector<Material *>::iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    energy += (*mat_it)->getEnergy(energy_id);
  }

  /// reduction sum over all processors
  StaticCommunicator::getStaticCommunicator().allReduce(&energy, 1, _so_sum);

  AKANTU_DEBUG_OUT();
  return energy;
}
/* -------------------------------------------------------------------------- */
Real SolidMechanicsModel::getEnergy(const std::string & energy_id,
                                    const ElementType & type,
                                    UInt index){
  AKANTU_DEBUG_IN();

  if (energy_id == "kinetic") {
    return getKineticEnergy(type, index);
  }

  std::vector<Material *>::iterator mat_it;
  Vector<UInt> mat = element_index_by_material(type, _not_ghost).begin(2)[index];
  Real energy = materials[mat(0)]->getEnergy(energy_id, type, mat(1));

  AKANTU_DEBUG_OUT();
  return energy;
}


/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::onNodesAdded(const Array<UInt> & nodes_list,
                                       __attribute__((unused)) const NewNodesEvent & event) {
  AKANTU_DEBUG_IN();
  UInt nb_nodes = mesh.getNbNodes();

  if(displacement) displacement->resize(nb_nodes);
  if(mass        ) mass        ->resize(nb_nodes);
  if(velocity    ) velocity    ->resize(nb_nodes);
  if(acceleration) acceleration->resize(nb_nodes);
  if(force       ) force       ->resize(nb_nodes);
  if(residual    ) residual    ->resize(nb_nodes);
  if(blocked_dofs) blocked_dofs->resize(nb_nodes);

  if(previous_displacement) previous_displacement->resize(nb_nodes);
  if(increment_acceleration) increment_acceleration->resize(nb_nodes);
  if(increment) increment->resize(nb_nodes);

  if(current_position) current_position->resize(nb_nodes);

  delete dof_synchronizer;
  dof_synchronizer = new DOFSynchronizer(mesh, spatial_dimension);
  dof_synchronizer->initLocalDOFEquationNumbers();
  dof_synchronizer->initGlobalDOFEquationNumbers();

  std::vector<Material *>::iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    (*mat_it)->onNodesAdded(nodes_list, event);
  }

  if (method != _explicit_lumped_mass) {
    delete stiffness_matrix;
    delete jacobian_matrix;
    delete solver;
    SolverOptions  solver_options;
    initImplicit((method == _implicit_dynamic), solver_options);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::onElementsAdded(const Array<Element> & element_list,
                                          const NewElementsEvent & event) {
  AKANTU_DEBUG_IN();

  getFEEngine().initShapeFunctions(_not_ghost);
  getFEEngine().initShapeFunctions(_ghost);

  Array<Element>::const_iterator<Element> it  = element_list.begin();
  Array<Element>::const_iterator<Element> end = element_list.end();

  /// \todo have rules to choose the correct material
  UInt mat_id = 0;

  UInt * mat_id_vect = NULL;
  try {
    const NewMaterialElementsEvent & event_mat = dynamic_cast<const NewMaterialElementsEvent &>(event);
    mat_id_vect = event_mat.getMaterialList().storage();
  } catch(...) {  }

  for (UInt el = 0; it != end; ++it, ++el) {
    const Element & elem = *it;

    if(mat_id_vect) mat_id = mat_id_vect[el];
    else mat_id = (*material_selector)(elem);

    Material & mat = *materials[mat_id];

    UInt mat_index = mat.addElement(elem.type, elem.element, elem.ghost_type);
    Vector<UInt> id(2);
    id[0] = mat_id; id[1] = mat_index;
    if(!element_index_by_material.exists(elem.type, elem.ghost_type))
      element_index_by_material.alloc(0, 2, elem.type, elem.ghost_type);

    element_index_by_material(elem.type, elem.ghost_type).push_back(id);
  }

  std::vector<Material *>::iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    (*mat_it)->onElementsAdded(element_list, event);
  }

  if(method != _explicit_lumped_mass) AKANTU_DEBUG_TO_IMPLEMENT();

  assembleMassLumped();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::onElementsRemoved(__attribute__((unused)) const Array<Element> & element_list,
                                            const ElementTypeMapArray<UInt> & new_numbering,
                                            const RemovedElementsEvent & event) {
  //  MeshUtils::purifyMesh(mesh);

  getFEEngine().initShapeFunctions(_not_ghost);
  getFEEngine().initShapeFunctions(_ghost);

  std::vector<Material *>::iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    (*mat_it)->onElementsRemoved(element_list, new_numbering, event);
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::onNodesRemoved(__attribute__((unused)) const Array<UInt> & element_list,
                                         const Array<UInt> & new_numbering,
                                         __attribute__((unused)) const RemovedNodesEvent & event) {
  if(displacement) mesh.removeNodesFromArray(*displacement, new_numbering);
  if(mass        ) mesh.removeNodesFromArray(*mass        , new_numbering);
  if(velocity    ) mesh.removeNodesFromArray(*velocity    , new_numbering);
  if(acceleration) mesh.removeNodesFromArray(*acceleration, new_numbering);
  if(force       ) mesh.removeNodesFromArray(*force       , new_numbering);
  if(residual    ) mesh.removeNodesFromArray(*residual    , new_numbering);
  if(blocked_dofs) mesh.removeNodesFromArray(*blocked_dofs, new_numbering);

  if(increment_acceleration) mesh.removeNodesFromArray(*increment_acceleration, new_numbering);
  if(increment)              mesh.removeNodesFromArray(*increment             , new_numbering);

  delete dof_synchronizer;
  dof_synchronizer = new DOFSynchronizer(mesh, spatial_dimension);
  dof_synchronizer->initLocalDOFEquationNumbers();
  dof_synchronizer->initGlobalDOFEquationNumbers();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::reassignMaterial() {
  AKANTU_DEBUG_IN();

  std::vector< Array<Element> > element_to_add   (materials.size());
  std::vector< Array<Element> > element_to_remove(materials.size());

  Element element;
  for (ghost_type_t::iterator gt = ghost_type_t::begin(); gt != ghost_type_t::end(); ++gt) {
    GhostType ghost_type = *gt;
    element.ghost_type = ghost_type;

    Mesh::type_iterator it  = mesh.firstType(spatial_dimension, ghost_type, _ek_regular);
    Mesh::type_iterator end = mesh.lastType(spatial_dimension, ghost_type, _ek_regular);
    for(; it != end; ++it) {
      ElementType type = *it;
      element.type = type;
      element.kind = Mesh::getKind(type);

      UInt nb_element = mesh.getNbElement(type, ghost_type);
      Array<UInt> & el_index_by_mat = element_index_by_material(type, ghost_type);

      for (UInt el = 0; el < nb_element; ++el) {
        element.element = el;

        UInt old_material = el_index_by_mat(el, 0);
        UInt new_material = (*material_selector)(element);

        if(old_material != new_material) {
          element_to_add   [new_material].push_back(element);
          element_to_remove[old_material].push_back(element);
        }
      }
    }
  }

  std::vector<Material *>::iterator mat_it;
  UInt mat_index = 0;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it, ++mat_index) {
    (*mat_it)->removeElements(element_to_remove[mat_index]);
    (*mat_it)->addElements   (element_to_add[mat_index]);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
bool SolidMechanicsModel::isInternal(const std::string & field_name,
				     const ElementKind & element_kind){

  bool is_internal = false;
  
  /// check if at least one material contains field_id as an internal
  for (UInt m = 0; m < materials.size() && !is_internal; ++m) {
    is_internal |= materials[m]->isInternal(field_name, element_kind);
  }
  
  return is_internal;
}
/* -------------------------------------------------------------------------- */

ElementTypeMap<UInt> SolidMechanicsModel::getInternalDataPerElem(const std::string & field_name,
								const ElementKind & element_kind){


  if (!(this->isInternal(field_name,element_kind))) AKANTU_EXCEPTION("unknown internal " << field_name);

  for (UInt m = 0; m < materials.size() ; ++m) {
    if (materials[m]->isInternal(field_name, element_kind))
      return materials[m]->getInternalDataPerElem(field_name,element_kind);
  }
  
  return ElementTypeMap<UInt>();
}


/* -------------------------------------------------------------------------- */
ElementTypeMapArray<Real> & SolidMechanicsModel::flattenInternal(const std::string & field_name,
								const ElementKind & kind){

  std::pair<std::string,ElementKind> key(field_name,kind);
  if (this->registered_internals.count(key) == 0){
    this->registered_internals[key] =
      new ElementTypeMapArray<Real>(field_name,this->id);
  }

  ElementTypeMapArray<Real> * internal_flat = this->registered_internals[key];
  for (UInt m = 0; m < materials.size(); ++m)
    materials[m]->flattenInternal(field_name,*internal_flat,_not_ghost,kind);
  
  return  *internal_flat;
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::flattenAllRegisteredInternals(const ElementKind & kind){

  std::map<std::pair<std::string,ElementKind>,ElementTypeMapArray<Real> *> ::iterator it = this->registered_internals.begin();
  std::map<std::pair<std::string,ElementKind>,ElementTypeMapArray<Real> *>::iterator end = this->registered_internals.end();

  while (it != end){
    if (kind == it->first.second)
      this->flattenInternal(it->first.first,kind);
    ++it;
  }
}

/* -------------------------------------------------------------------------- */

void SolidMechanicsModel::onDump(){
  this->flattenAllRegisteredInternals(_ek_regular);
}
/* -------------------------------------------------------------------------- */

#ifdef AKANTU_USE_IOHELPER

dumper::Field * SolidMechanicsModel
::createElementalField(const std::string & field_name, 
		       const std::string & group_name,
		       bool padding_flag,
		       const ElementKind & kind){
  
  dumper::Field * field = NULL;

  if(field_name == "partitions") 
    field = mesh.createElementalField<UInt, dumper::ElementPartitionField>(mesh.getConnectivities(),group_name,this->spatial_dimension,kind);
  else if(field_name == "element_index_by_material") 
    field = mesh.createElementalField<UInt, Vector, dumper::ElementalField >(element_index_by_material,group_name,this->spatial_dimension,kind);
  else {

    bool is_internal = this->isInternal(field_name,kind);

    if (is_internal) {

      ElementTypeMap<UInt> nb_data_per_elem = this->getInternalDataPerElem(field_name,kind);
      ElementTypeMapArray<Real> & internal_flat = this->flattenInternal(field_name,kind);
      field = mesh.createElementalField<Real, dumper::InternalMaterialField>(internal_flat,
									     group_name,
									     this->spatial_dimension,kind,nb_data_per_elem);
      
      //treat the paddings
      if (padding_flag){
	if (field_name == "stress"){
	  if (this->spatial_dimension == 2) {
	    dumper::StressPadder<2> * foo = new dumper::StressPadder<2>(*this);
	    field = dumper::FieldComputeProxy::createFieldCompute(field,*foo);
	  }
	}
	// else if (field_name == "strain"){
	//   if (this->spatial_dimension == 2) {
	//     dumper::StrainPadder<2> * foo = new dumper::StrainPadder<2>(*this);
	//     field = dumper::FieldComputeProxy::createFieldCompute(field,*foo);
	//   }
	// }
      }

      // homogenize the field
      dumper::ComputeFunctorInterface * foo = 
	dumper::HomogenizerProxy::createHomogenizer(*field);

      field = dumper::FieldComputeProxy::createFieldCompute(field,*foo);

    } 
  }
  return field;
}

/* -------------------------------------------------------------------------- */

dumper::Field * SolidMechanicsModel::createNodalFieldReal(const std::string & field_name,
							  const std::string & group_name,
							  bool padding_flag) {
  
  std::map<std::string,Array<Real>* > real_nodal_fields;
  real_nodal_fields["displacement"             ] = displacement;
  real_nodal_fields["mass"                     ] = mass;
  real_nodal_fields["velocity"                 ] = velocity;
  real_nodal_fields["acceleration"             ] = acceleration;
  real_nodal_fields["force"                    ] = force;
  real_nodal_fields["residual"                 ] = residual;
  real_nodal_fields["increment"                ] = increment;

  dumper::Field * field = NULL;
  if (padding_flag)
    field = mesh.createNodalField(real_nodal_fields[field_name],group_name, 3);
  else
    field = mesh.createNodalField(real_nodal_fields[field_name],group_name);
  
  return field;
}
/* -------------------------------------------------------------------------- */

dumper::Field * SolidMechanicsModel::createNodalFieldBool(const std::string & field_name,
							  const std::string & group_name,
bool padding_flag) {
  
  std::map<std::string,Array<bool>* > uint_nodal_fields;
  uint_nodal_fields["blocked_dofs"             ] = blocked_dofs;

  dumper::Field * field = NULL;
  field = mesh.createNodalField(uint_nodal_fields[field_name],group_name);
  return field;

}
/* -------------------------------------------------------------------------- */
#else
/* -------------------------------------------------------------------------- */

dumper::Field * SolidMechanicsModel
::createElementalField(const std::string & field_name, 
		       const std::string & group_name,
		       bool padding_flag,
		       const ElementKind & kind){
  return NULL;
}

/* -------------------------------------------------------------------------- */

dumper::Field * SolidMechanicsModel::createNodalFieldReal(const std::string & field_name,
							  const std::string & group_name,
							  bool padding_flag) {
  return NULL;
}

/* -------------------------------------------------------------------------- */

dumper::Field * SolidMechanicsModel::createNodalFieldBool(const std::string & field_name,
							  const std::string & group_name,
bool padding_flag) {
  return NULL;
}

#endif
/* -------------------------------------------------------------------------- */

void SolidMechanicsModel::dump(const std::string & dumper_name) {
  this->onDump();
  EventManager::sendEvent(SolidMechanicsModelEvent::BeforeDumpEvent());
  synch_registry->synchronize(_gst_for_dump);
  mesh.dump(dumper_name);
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::dump(const std::string & dumper_name, UInt step) {
  this->onDump();
  EventManager::sendEvent(SolidMechanicsModelEvent::BeforeDumpEvent());
  synch_registry->synchronize(_gst_for_dump);
  mesh.dump(dumper_name, step);
}

/* ------------------------------------------------------------------------- */
void SolidMechanicsModel::dump(const std::string & dumper_name, Real time, UInt step) {
  this->onDump();
  EventManager::sendEvent(SolidMechanicsModelEvent::BeforeDumpEvent());
  synch_registry->synchronize(_gst_for_dump);
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
void SolidMechanicsModel::computeCauchyStresses() {
  AKANTU_DEBUG_IN();

  // call compute stiffness matrix on each local elements
  std::vector<Material *>::iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    Material & mat = **mat_it;
    if(mat.isFiniteDeformation())
      mat.computeAllCauchyStresses(_not_ghost);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::saveStressAndStrainBeforeDamage() {
  EventManager::sendEvent(SolidMechanicsModelEvent::BeginningOfDamageIterationEvent());
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::updateEnergiesAfterDamage() {
  EventManager::sendEvent(SolidMechanicsModelEvent::AfterDamageEvent());
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "Solid Mechanics Model [" << std::endl;
  stream << space << " + id                : " << id << std::endl;
  stream << space << " + spatial dimension : " << spatial_dimension << std::endl;
  stream << space << " + fem [" << std::endl;
  getFEEngine().printself(stream, indent + 2);
  stream << space << AKANTU_INDENT << "]" << std::endl;
  stream << space << " + nodals information [" << std::endl;
  displacement->printself(stream, indent + 2);
  mass        ->printself(stream, indent + 2);
  velocity    ->printself(stream, indent + 2);
  acceleration->printself(stream, indent + 2);
  force       ->printself(stream, indent + 2);
  residual    ->printself(stream, indent + 2);
  blocked_dofs->printself(stream, indent + 2);
  stream << space << AKANTU_INDENT << "]" << std::endl;

  stream << space << " + connectivity type information [" << std::endl;
  element_index_by_material.printself(stream, indent + 2);
  stream << space << AKANTU_INDENT << "]" << std::endl;

  stream << space << " + materials [" << std::endl;
  std::vector<Material *>::const_iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    const Material & mat = *(*mat_it);
    mat.printself(stream, indent + 1);
  }
  stream << space << AKANTU_INDENT << "]" << std::endl;

  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__
