/**
 * @file   structural_mechanics_model.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Damien Spielmann <damien.spielmann@epfl.ch>
 * @author Sébastien Hartmann <sebastien.hartmann@epfl.ch>
 * @author Fabian Barras <fabian.barras@epfl.ch>
 *
 * @date creation: Fri Jul 15 2011
 * @date last modification: Mon Sep 15 2014
 *
 * @brief  Model implementation for StucturalMechanics elements
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "structural_mechanics_model.hh"
#include "group_manager_inline_impl.cc"
#include "dumpable_inline_impl.hh"
#include "aka_math.hh"
#include "integration_scheme_2nd_order.hh"
#include "static_communicator.hh"
#include "sparse_matrix.hh"
#include "solver.hh"

#ifdef AKANTU_USE_MUMPS
#include "solver_mumps.hh"
#endif

#ifdef AKANTU_USE_IOHELPER
#  include "dumper_paraview.hh"
#  include "dumper_elemental_field.hh" 
#endif
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

const StructuralMechanicsModelOptions default_structural_mechanics_model_options(_static);

/* -------------------------------------------------------------------------- */
StructuralMechanicsModel::StructuralMechanicsModel(Mesh & mesh,
						   UInt dim,
						   const ID & id,
						   const MemoryID & memory_id) :
  Model(mesh, dim, id, memory_id), 
  time_step(NAN), f_m2a(1.0),
  stress("stress", id, memory_id),
  element_material("element_material", id, memory_id),
  set_ID("beam sets", id, memory_id),
  stiffness_matrix(NULL),
  mass_matrix(NULL),
  velocity_damping_matrix(NULL),
  jacobian_matrix(NULL),
  solver(NULL),
  rotation_matrix("rotation_matices", id, memory_id),
  increment_flag(false),
  integrator(NULL) {
  AKANTU_DEBUG_IN();

  registerFEEngineObject<MyFEEngineType>("StructuralMechanicsFEEngine", mesh, spatial_dimension);

  this->displacement_rotation = NULL;
  this->mass                  = NULL;
  this->velocity              = NULL;
  this->acceleration          = NULL;
  this->force_momentum        = NULL;
  this->residual              = NULL;
  this->blocked_dofs              = NULL;
  this->increment             = NULL;

  this->previous_displacement = NULL;

  if(spatial_dimension == 2)
    nb_degree_of_freedom = 3;
  else if (spatial_dimension == 3)
    nb_degree_of_freedom = 6;
  else {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  this->mesh.registerDumper<DumperParaview>("paraview_all", id, true);
  this->mesh.addDumpMesh(mesh, spatial_dimension, _not_ghost, _ek_structural);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
StructuralMechanicsModel::~StructuralMechanicsModel() {
  AKANTU_DEBUG_IN();

  delete integrator;
  delete solver;
  delete stiffness_matrix;
  delete jacobian_matrix;
  delete mass_matrix;
  delete velocity_damping_matrix;

  AKANTU_DEBUG_OUT();
}

void StructuralMechanicsModel::setTimeStep(Real time_step) {
  this->time_step = time_step;

#if defined(AKANTU_USE_IOHELPER)
  this->mesh.getDumper().setTimeStep(time_step);
#endif
}

/* -------------------------------------------------------------------------- */
/* Initialisation                                                             */
/* -------------------------------------------------------------------------- */

void StructuralMechanicsModel::initFull(const ModelOptions & options) {
  
  const StructuralMechanicsModelOptions & stmm_options =
    dynamic_cast<const StructuralMechanicsModelOptions &>(options);

  method = stmm_options.analysis_method;
  
  initModel();
  initArrays();

  displacement_rotation->clear();
  velocity             ->clear();
  acceleration         ->clear();
  force_momentum       ->clear();
  residual             ->clear();
  blocked_dofs             ->clear();
  increment            ->clear();


  Mesh::type_iterator it = getFEEngine().getMesh().firstType(spatial_dimension, _not_ghost, _ek_structural);
  Mesh::type_iterator end = getFEEngine().getMesh().lastType(spatial_dimension, _not_ghost, _ek_structural);
  for (; it != end; ++it) {
    computeRotationMatrix(*it);
  }

  switch(method) {
  case _implicit_dynamic:
    initImplicit();
    break;
  case _static:
    initSolver();
    break;
  default:
    AKANTU_EXCEPTION("analysis method not recognised by StructuralMechanicsModel");
      break;
  }
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::initArrays() {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = getFEEngine().getMesh().getNbNodes();
  std::stringstream sstr_disp; sstr_disp << id << ":displacement";
  std::stringstream sstr_velo; sstr_velo << id << ":velocity";
  std::stringstream sstr_acce; sstr_acce << id << ":acceleration";
  std::stringstream sstr_forc; sstr_forc << id << ":force";
  std::stringstream sstr_resi; sstr_resi << id << ":residual";
  std::stringstream sstr_boun; sstr_boun << id << ":blocked_dofs";
  std::stringstream sstr_incr; sstr_incr << id << ":increment";

  displacement_rotation = &(alloc<Real>(sstr_disp.str(), nb_nodes, nb_degree_of_freedom, REAL_INIT_VALUE));
  velocity              = &(alloc<Real>(sstr_velo.str(), nb_nodes, nb_degree_of_freedom, REAL_INIT_VALUE));
  acceleration          = &(alloc<Real>(sstr_acce.str(), nb_nodes, nb_degree_of_freedom, REAL_INIT_VALUE));
  force_momentum        = &(alloc<Real>(sstr_forc.str(), nb_nodes, nb_degree_of_freedom, REAL_INIT_VALUE));
  residual              = &(alloc<Real>(sstr_resi.str(), nb_nodes, nb_degree_of_freedom, REAL_INIT_VALUE));
  blocked_dofs              = &(alloc<bool>(sstr_boun.str(), nb_nodes, nb_degree_of_freedom, false));
  increment             = &(alloc<Real>(sstr_incr.str(), nb_nodes, nb_degree_of_freedom, REAL_INIT_VALUE));

  Mesh::type_iterator it = getFEEngine().getMesh().firstType(spatial_dimension, _not_ghost, _ek_structural);
  Mesh::type_iterator end = getFEEngine().getMesh().lastType(spatial_dimension, _not_ghost, _ek_structural);
  for (; it != end; ++it) {
    UInt nb_element = getFEEngine().getMesh().getNbElement(*it);
    UInt nb_quadrature_points       = getFEEngine().getNbQuadraturePoints(*it);

    element_material.alloc(nb_element, 1, *it, _not_ghost);
    set_ID.alloc(nb_element, 1, *it, _not_ghost);

    UInt size = getTangentStiffnessVoigtSize(*it);
    stress.alloc(nb_element * nb_quadrature_points, size , *it, _not_ghost);
  }

  dof_synchronizer = new DOFSynchronizer(getFEEngine().getMesh(), nb_degree_of_freedom);
  dof_synchronizer->initLocalDOFEquationNumbers();
  dof_synchronizer->initGlobalDOFEquationNumbers();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::initModel() {
  getFEEngine().initShapeFunctions(_not_ghost);
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::initSolver(__attribute__((unused)) SolverOptions & options) {
  AKANTU_DEBUG_IN();

  const Mesh & mesh = getFEEngine().getMesh();
#if !defined(AKANTU_USE_MUMPS) // or other solver in the future \todo add AKANTU_HAS_SOLVER in CMake
  AKANTU_DEBUG_ERROR("You should at least activate one solver.");
#else
  UInt nb_global_node = mesh.getNbGlobalNodes();

  std::stringstream sstr; sstr << id << ":jacobian_matrix";
  jacobian_matrix = new SparseMatrix(nb_global_node * nb_degree_of_freedom, _symmetric, sstr.str(), memory_id);

  jacobian_matrix->buildProfile(mesh, *dof_synchronizer, nb_degree_of_freedom);

  std::stringstream sstr_sti; sstr_sti << id << ":stiffness_matrix";
  stiffness_matrix = new SparseMatrix(*jacobian_matrix, sstr_sti.str(), memory_id);

#ifdef AKANTU_USE_MUMPS
  std::stringstream sstr_solv; sstr_solv << id << ":solver";
  solver = new SolverMumps(*jacobian_matrix, sstr_solv.str());

  dof_synchronizer->initScatterGatherCommunicationScheme();
#else
  AKANTU_DEBUG_ERROR("You should at least activate one solver.");
#endif //AKANTU_USE_MUMPS

  solver->initialize(options);
#endif //AKANTU_HAS_SOLVER

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
 void StructuralMechanicsModel::initImplicit(bool dynamic, SolverOptions & solver_options) {
   AKANTU_DEBUG_IN();
 
   if (!increment) setIncrementFlagOn();

   initSolver(solver_options);

   if(integrator) delete integrator;
   integrator = new TrapezoidalRule2();

   AKANTU_DEBUG_OUT();
 }

/* -------------------------------------------------------------------------- */
UInt StructuralMechanicsModel::getTangentStiffnessVoigtSize(const ElementType & type) {
  UInt size;
#define GET_TANGENT_STIFFNESS_VOIGT_SIZE(type)	\
  size = getTangentStiffnessVoigtSize<type>();

  AKANTU_BOOST_STRUCTURAL_ELEMENT_SWITCH(GET_TANGENT_STIFFNESS_VOIGT_SIZE);
#undef GET_TANGENT_STIFFNESS_VOIGT_SIZE

  return size;
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::assembleStiffnessMatrix() {
  AKANTU_DEBUG_IN();

  stiffness_matrix->clear();

  Mesh::type_iterator it = getFEEngine().getMesh().firstType(spatial_dimension, _not_ghost, _ek_structural);
  Mesh::type_iterator end = getFEEngine().getMesh().lastType(spatial_dimension, _not_ghost, _ek_structural);
  for (; it != end; ++it) {
    ElementType type = *it;

#define ASSEMBLE_STIFFNESS_MATRIX(type)		\
    assembleStiffnessMatrix<type>();

    AKANTU_BOOST_STRUCTURAL_ELEMENT_SWITCH(ASSEMBLE_STIFFNESS_MATRIX);
#undef ASSEMBLE_STIFFNESS_MATRIX
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<>
void StructuralMechanicsModel::computeRotationMatrix<_bernoulli_beam_2>(Array<Real> & rotations){

  ElementType type = _bernoulli_beam_2;
  Mesh & mesh = getFEEngine().getMesh();
  UInt nb_element = mesh.getNbElement(type);

  Array<UInt>::iterator< Vector<UInt> > connec_it = mesh.getConnectivity(type).begin(2);
  Array<Real>::vector_iterator nodes_it = mesh.getNodes().begin(spatial_dimension);
  Array<Real>::matrix_iterator R_it = rotations.begin(nb_degree_of_freedom, nb_degree_of_freedom);

  for (UInt e = 0; e < nb_element; ++e, ++R_it, ++connec_it) {
    Matrix<Real> & R = *R_it;
    Vector<UInt> & connec = *connec_it;

    Vector<Real> x2;
    x2 = nodes_it[connec(1)]; // X2
    Vector<Real> x1;
    x1 = nodes_it[connec(0)]; // X1

    Real le = x1.distance(x2);
    Real c = (x2(0) - x1(0)) / le;
    Real s = (x2(1) - x1(1)) / le;

    /// Definition of the rotation matrix
    R(0,0) =  c;  R(0,1) = s;  R(0,2) = 0.;
    R(1,0) = -s;  R(1,1) = c;  R(1,2) = 0.;
    R(2,0) =  0.; R(2,1) = 0.; R(2,2) = 1.;
  }
}


/* -------------------------------------------------------------------------- */
template<>
void StructuralMechanicsModel::computeRotationMatrix<_bernoulli_beam_3>(Array<Real> & rotations){
  ElementType type = _bernoulli_beam_3;
  Mesh & mesh = getFEEngine().getMesh();
  UInt nb_element = mesh.getNbElement(type);

  Array<Real>::vector_iterator n_it = mesh.getNormals(type).begin(spatial_dimension);
  Array<UInt>::iterator< Vector<UInt> > connec_it = mesh.getConnectivity(type).begin(2);
  Array<Real>::vector_iterator nodes_it = mesh.getNodes().begin(spatial_dimension);

  Matrix<Real> Pe    (spatial_dimension, spatial_dimension);
  Matrix<Real> Pg    (spatial_dimension, spatial_dimension);
  Matrix<Real> inv_Pg(spatial_dimension, spatial_dimension);
  Vector<Real> x_n(spatial_dimension); // x vect n

  Array<Real>::matrix_iterator R_it =
    rotations.begin(nb_degree_of_freedom, nb_degree_of_freedom);

  for (UInt e=0 ; e < nb_element; ++e, ++n_it, ++connec_it, ++R_it) {
    Vector<Real> & n = *n_it;
    Matrix<Real> & R = *R_it;
    Vector<UInt> & connec = *connec_it;

    Vector<Real> x;
    x = nodes_it[connec(1)]; // X2
    Vector<Real> y;
    y = nodes_it[connec(0)]; // X1

    Real l = x.distance(y);
    x -= y; // X2 - X1
    x_n.crossProduct(x, n);

    Pe.eye();
    Pe(0, 0) *=  l;
    Pe(1, 1) *= -l;

    Pg(0,0) = x(0); Pg(0,1) = x_n(0); Pg(0,2) = n(0);
    Pg(1,0) = x(1); Pg(1,1) = x_n(1); Pg(1,2) = n(1);
    Pg(2,0) = x(2); Pg(2,1) = x_n(2); Pg(2,2) = n(2);

    inv_Pg.inverse(Pg);

    Pe *= inv_Pg;
    for (UInt i = 0; i < spatial_dimension; ++i) {
      for (UInt j = 0; j < spatial_dimension; ++j) {
	R(i, j) = Pe(i, j);
	R(i + spatial_dimension,j + spatial_dimension) = Pe(i, j);
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
template<> 
void StructuralMechanicsModel::computeRotationMatrix<_kirchhoff_shell>(Array<Real> & rotations){
  ElementType type = _kirchhoff_shell;
  Mesh & mesh = getFEEngine().getMesh();
  UInt nb_element = mesh.getNbElement(type);
  
  Array<UInt>::iterator< Vector<UInt> > connec_it = mesh.getConnectivity(type).begin(3);
  Array<Real>::vector_iterator nodes_it = mesh.getNodes().begin(spatial_dimension);
  
  Matrix<Real> Pe    (spatial_dimension, spatial_dimension);
  Matrix<Real> Pg    (spatial_dimension, spatial_dimension);
  Matrix<Real> inv_Pg(spatial_dimension, spatial_dimension);
  
  Array<Real>::matrix_iterator R_it =
    rotations.begin(nb_degree_of_freedom, nb_degree_of_freedom);
  
  for (UInt e=0 ; e < nb_element; ++e, ++connec_it, ++R_it) { 
      
    Pe.eye();
    
    Matrix<Real> & R = *R_it;
    Vector<UInt> & connec = *connec_it;

    Vector<Real> x2;
    x2 = nodes_it[connec(1)]; // X2
    Vector<Real> x1;
    x1 = nodes_it[connec(0)]; // X1
    Vector<Real> x3;
    x3 = nodes_it[connec(2)]; // X3




    Vector<Real> Pg_col_1=x2-x1;
    
    Vector<Real> Pg_col_2 = x3-x1;
    
    Vector<Real> Pg_col_3(spatial_dimension);
    Pg_col_3.crossProduct(Pg_col_1, Pg_col_2);
      
      
    for (UInt i = 0; i <spatial_dimension; ++i) {
      Pg(i,0) = Pg_col_1(i);
      Pg(i,1) = Pg_col_2(i);
      Pg(i,2) = Pg_col_3(i);
    }
  

    inv_Pg.inverse(Pg);
    // Pe *= inv_Pg;
    Pe.eye(); 
    
    for (UInt i = 0; i < spatial_dimension; ++i) {
      for (UInt j = 0; j < spatial_dimension; ++j) {
	R(i, j) = Pe(i, j);
	R(i + spatial_dimension,j + spatial_dimension) = Pe(i, j);
      }
    }  

  }
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::computeRotationMatrix(const ElementType & type) {
  Mesh & mesh = getFEEngine().getMesh();

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_element = mesh.getNbElement(type);

  if(!rotation_matrix.exists(type)) {
    rotation_matrix.alloc(nb_element,
			  nb_degree_of_freedom*nb_nodes_per_element * nb_degree_of_freedom*nb_nodes_per_element,
			  type);
  } else {
    rotation_matrix(type).resize(nb_element);
  }
  rotation_matrix(type).clear();

  Array<Real>rotations(nb_element, nb_degree_of_freedom * nb_degree_of_freedom);
  rotations.clear();

#define COMPUTE_ROTATION_MATRIX(type)	\
  computeRotationMatrix<type>(rotations);

  AKANTU_BOOST_STRUCTURAL_ELEMENT_SWITCH(COMPUTE_ROTATION_MATRIX);
#undef COMPUTE_ROTATION_MATRIX


  Array<Real>::matrix_iterator R_it = rotations.begin(nb_degree_of_freedom, nb_degree_of_freedom);
  Array<Real>::matrix_iterator T_it =
    rotation_matrix(type).begin(nb_degree_of_freedom*nb_nodes_per_element,
				nb_degree_of_freedom*nb_nodes_per_element);

  for (UInt el = 0; el < nb_element; ++el, ++R_it, ++T_it) {
    Matrix<Real> & T = *T_it;
    Matrix<Real> & R = *R_it;
    T.clear();
    for (UInt k = 0; k < nb_nodes_per_element; ++k){
      for (UInt i = 0; i < nb_degree_of_freedom; ++i)
	for (UInt j = 0; j < nb_degree_of_freedom; ++j)
	  T(k*nb_degree_of_freedom + i, k*nb_degree_of_freedom + j) = R(i, j);
    }
  }
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::computeStresses() {
  AKANTU_DEBUG_IN();

  Mesh::type_iterator it = getFEEngine().getMesh().firstType(spatial_dimension, _not_ghost, _ek_structural);
  Mesh::type_iterator end = getFEEngine().getMesh().lastType(spatial_dimension, _not_ghost, _ek_structural);
  for (; it != end; ++it) {
    ElementType type = *it;

#define COMPUTE_STRESS_ON_QUAD(type)		\
    computeStressOnQuad<type>();

    AKANTU_BOOST_STRUCTURAL_ELEMENT_SWITCH(COMPUTE_STRESS_ON_QUAD);
#undef COMPUTE_STRESS_ON_QUAD
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::updateResidual() {
 AKANTU_DEBUG_IN();
 residual->copy(*force_momentum);

 Array<Real> ku(*displacement_rotation, true);
 ku *= *stiffness_matrix;
 *residual -= ku;

 this->updateResidualInternal();

 AKANTU_DEBUG_OUT();

}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::implicitPred() {
  AKANTU_DEBUG_IN();

  if(previous_displacement) previous_displacement->copy(*displacement_rotation);

  if(method == _implicit_dynamic)
    integrator->integrationSchemePred(time_step,
                                    *displacement_rotation,
                                    *velocity,
                                    *acceleration,
                                    *blocked_dofs);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::implicitCorr() {
  AKANTU_DEBUG_IN();
  
  Real * incr_val = increment->storage();
  bool * boun_val = blocked_dofs->storage();

  UInt nb_nodes = displacement_rotation->getSize();
  
  for (UInt j = 0; j < nb_nodes * nb_degree_of_freedom; ++j, ++incr_val, ++boun_val){
      *incr_val *= (1. - *boun_val);
  }

  if(method == _implicit_dynamic) {
    integrator->integrationSchemeCorrDispl(time_step,
                                           *displacement_rotation,
                                           *velocity,
                                           *acceleration,
                                           *blocked_dofs,
                                           *increment);
  } else {

    Real * disp_val = displacement_rotation->storage();
    incr_val = increment->storage();

    for (UInt j = 0; j < nb_nodes *nb_degree_of_freedom; ++j, ++disp_val){
      *disp_val += *incr_val;
    }
  }

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::updateResidualInternal() {
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
void StructuralMechanicsModel::setIncrementFlagOn() {
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
void StructuralMechanicsModel::solve() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Solving an implicit step.");

  UInt nb_nodes = displacement_rotation->getSize();

  /// todo residual = force - Kxr * d_bloq
  jacobian_matrix->copyContent(*stiffness_matrix);
  jacobian_matrix->applyBoundary(*blocked_dofs);

  jacobian_matrix->saveMatrix("Kbound.mtx");

  increment->clear();

  solver->setRHS(*residual);

  solver->factorize();
  solver->solve(*increment);

  Real * increment_val     = increment->storage();
  Real * displacement_val  = displacement_rotation->storage();
  bool * blocked_dofs_val      = blocked_dofs->storage();

  for (UInt n = 0; n < nb_nodes * nb_degree_of_freedom; ++n) {
    if(!(*blocked_dofs_val)) {
      *displacement_val += *increment_val;
    }
    else {
      *increment_val = 0.0;
    }

    displacement_val++;
    blocked_dofs_val++;
    increment_val++;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<>
bool StructuralMechanicsModel::testConvergence<_scc_increment>(Real tolerance, Real & error){
  AKANTU_DEBUG_IN();

  UInt nb_nodes = displacement_rotation->getSize();
  UInt nb_degree_of_freedom = displacement_rotation->getNbComponent();

  error = 0;
  Real norm[2] = {0., 0.};
  Real * increment_val    = increment->storage();
  bool * blocked_dofs_val     = blocked_dofs->storage();
  Real * displacement_val = displacement_rotation->storage();

  for (UInt n = 0; n < nb_nodes; ++n) {
    for (UInt d = 0; d < nb_degree_of_freedom; ++d) {
      if(!(*blocked_dofs_val)) {
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
  AKANTU_DEBUG_ASSERT(!Math::isnan(norm[0]), "Something went wrong in the solve phase");

  AKANTU_DEBUG_OUT();
  if(norm[1] > Math::getTolerance())
    error = norm[0] / norm[1];
  else
    error = norm[0]; //In case the total displacement is zero!
  return (error < tolerance);
}

/* -------------------------------------------------------------------------- */
template<>
bool StructuralMechanicsModel::testConvergence<_scc_residual>(Real tolerance, Real & norm) {
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
bool StructuralMechanicsModel::testConvergenceIncrement(Real tolerance) {
  Real error;
  bool tmp = testConvergenceIncrement(tolerance, error);

  AKANTU_DEBUG_INFO("Norm of increment : " << error);

  return tmp;
}

/* -------------------------------------------------------------------------- */
bool StructuralMechanicsModel::testConvergenceIncrement(Real tolerance, Real & error) {
  AKANTU_DEBUG_IN();
  Mesh & mesh= getFEEngine().getMesh();
  UInt nb_nodes = displacement_rotation->getSize();
  UInt nb_degree_of_freedom = displacement_rotation->getNbComponent();

  Real norm = 0;
  Real * increment_val     = increment->storage();
  bool * blocked_dofs_val      = blocked_dofs->storage();

  for (UInt n = 0; n < nb_nodes; ++n) {
    bool is_local_node = mesh.isLocalOrMasterNode(n);
    for (UInt d = 0; d < nb_degree_of_freedom; ++d) {
      if(!(*blocked_dofs_val) && is_local_node) {
	norm += *increment_val * *increment_val;
      }
      blocked_dofs_val++;
      increment_val++;
    }
  }

  StaticCommunicator::getStaticCommunicator().allReduce(&norm, 1, _so_sum);

  error = sqrt(norm);
  AKANTU_DEBUG_ASSERT(!Math::isnan(norm), "Something goes wrong in the solve phase");

  AKANTU_DEBUG_OUT();
  return (error < tolerance);
}

/* -------------------------------------------------------------------------- */
template<>
void StructuralMechanicsModel::computeTangentModuli<_bernoulli_beam_2>(Array<Real> & tangent_moduli) {
  UInt nb_element                 = getFEEngine().getMesh().getNbElement(_bernoulli_beam_2);
  UInt nb_quadrature_points       = getFEEngine().getNbQuadraturePoints(_bernoulli_beam_2);
  UInt tangent_size = 2;

  Array<Real>::matrix_iterator D_it = tangent_moduli.begin(tangent_size, tangent_size);
  Array<UInt> & el_mat = element_material(_bernoulli_beam_2, _not_ghost);

  for (UInt e = 0; e < nb_element; ++e) {
    UInt mat = el_mat(e);
    Real E = materials[mat].E;
    Real A = materials[mat].A;
    Real I = materials[mat].I;
    for (UInt q = 0; q < nb_quadrature_points; ++q, ++D_it) {
      Matrix<Real> & D = *D_it;
      D(0,0) = E * A;
      D(1,1) = E * I;
    }
  }
}
/* -------------------------------------------------------------------------- */
template<>
void StructuralMechanicsModel::computeTangentModuli<_bernoulli_beam_3>(Array<Real> & tangent_moduli) {
  UInt nb_element                 = getFEEngine().getMesh().getNbElement(_bernoulli_beam_3);
  UInt nb_quadrature_points       = getFEEngine().getNbQuadraturePoints(_bernoulli_beam_3);
  UInt tangent_size = 4;

  Array<Real>::matrix_iterator D_it = tangent_moduli.begin(tangent_size, tangent_size);

  for (UInt e = 0; e < nb_element; ++e) {
    UInt mat = element_material(_bernoulli_beam_3, _not_ghost)(e);
    Real E = materials[mat].E;
    Real A = materials[mat].A;
    Real Iz = materials[mat].Iz;
    Real Iy = materials[mat].Iy;
    Real GJ = materials[mat].GJ;
    for (UInt q = 0; q < nb_quadrature_points; ++q, ++D_it) {
      Matrix<Real> & D = *D_it;
      D(0,0) = E * A;
      D(1,1) = E * Iz;
      D(2,2) = E * Iy;
      D(3,3) = GJ;
    }
  }
}

/* -------------------------------------------------------------------------- */
template<>
void StructuralMechanicsModel::computeTangentModuli<_kirchhoff_shell>(Array<Real> & tangent_moduli) {
  UInt nb_element                 = getFEEngine().getMesh().getNbElement(_kirchhoff_shell);
  UInt nb_quadrature_points       = getFEEngine().getNbQuadraturePoints(_kirchhoff_shell);
  UInt tangent_size = 6;

  Array<Real>::matrix_iterator D_it = tangent_moduli.begin(tangent_size, tangent_size);

  for (UInt e = 0; e < nb_element; ++e) {
    UInt mat = element_material(_kirchhoff_shell, _not_ghost)(e);
    Real E = materials[mat].E;
    Real nu = materials[mat].nu;
    Real t = materials[mat].t;

    Real HH=E*t/(1-nu*nu);
    Real DD=E*t*t*t/(12*(1-nu*nu));
    
    for (UInt q = 0; q < nb_quadrature_points; ++q, ++D_it) {
      Matrix<Real> & D = *D_it;
      D(0,0) = HH;
      D(0,1) = HH*nu;
      D(1,0) = HH*nu;      
      D(1,1) = HH;
      D(2,2) = HH*(1-nu)/2;
      D(3,3) = DD;
      D(3,4) = DD*nu;
      D(4,3) = DD*nu;
      D(4,4) = DD;
      D(5,5) = DD*(1-nu)/2;
    }
  }
}

/* -------------------------------------------------------------------------- */
template<>
void StructuralMechanicsModel::transferBMatrixToSymVoigtBMatrix<_bernoulli_beam_2>(Array<Real> & b, bool local) {
  UInt nb_element                 = getFEEngine().getMesh().getNbElement(_bernoulli_beam_2);
  UInt nb_nodes_per_element       = Mesh::getNbNodesPerElement(_bernoulli_beam_2);
  UInt nb_quadrature_points       = getFEEngine().getNbQuadraturePoints(_bernoulli_beam_2);

  MyFEEngineType & fem = getFEEngineClass<MyFEEngineType>();
  Array<Real>::const_vector_iterator shape_Np  = fem.getShapesDerivatives(_bernoulli_beam_2, _not_ghost, 0).begin(nb_nodes_per_element);
  Array<Real>::const_vector_iterator shape_Mpp = fem.getShapesDerivatives(_bernoulli_beam_2, _not_ghost, 1).begin(nb_nodes_per_element);
  Array<Real>::const_vector_iterator shape_Lpp = fem.getShapesDerivatives(_bernoulli_beam_2, _not_ghost, 2).begin(nb_nodes_per_element);

  UInt tangent_size = getTangentStiffnessVoigtSize<_bernoulli_beam_2>();
  UInt bt_d_b_size  = nb_nodes_per_element * nb_degree_of_freedom;
  b.clear();
  Array<Real>::matrix_iterator B_it = b.begin(tangent_size, bt_d_b_size);

  for (UInt e = 0; e < nb_element; ++e) {
    for (UInt q = 0; q < nb_quadrature_points; ++q) {
      Matrix<Real> & B = *B_it;
      const Vector<Real> & Np  = *shape_Np;
      const Vector<Real> & Lpp = *shape_Lpp;
      const Vector<Real> & Mpp = *shape_Mpp;

      B(0,0) = Np(0);
      B(0,3) = Np(1);

      B(1,1) = Mpp(0);
      B(1,2) = Lpp(0);
      B(1,4) = Mpp(1);
      B(1,5) = Lpp(1);

      ++B_it;
      ++shape_Np;
      ++shape_Mpp;
      ++shape_Lpp;
    }

    // ++R_it;
  }
}

/* -------------------------------------------------------------------------- */
template<> 
void StructuralMechanicsModel::transferBMatrixToSymVoigtBMatrix<_bernoulli_beam_3>(Array<Real> & b, bool local) {
  MyFEEngineType & fem = getFEEngineClass<MyFEEngineType>();

  UInt nb_element                 = getFEEngine().getMesh().getNbElement(_bernoulli_beam_3);
  UInt nb_nodes_per_element       = Mesh::getNbNodesPerElement(_bernoulli_beam_3);
  UInt nb_quadrature_points       = getFEEngine().getNbQuadraturePoints(_bernoulli_beam_3);

  Array<Real>::const_vector_iterator shape_Np  = fem.getShapesDerivatives(_bernoulli_beam_3, _not_ghost, 0).begin(nb_nodes_per_element);
  Array<Real>::const_vector_iterator shape_Mpp = fem.getShapesDerivatives(_bernoulli_beam_3, _not_ghost, 1).begin(nb_nodes_per_element);
  Array<Real>::const_vector_iterator shape_Lpp = fem.getShapesDerivatives(_bernoulli_beam_3, _not_ghost, 2).begin(nb_nodes_per_element);

  UInt tangent_size = getTangentStiffnessVoigtSize<_bernoulli_beam_3>();
  UInt bt_d_b_size  = nb_nodes_per_element * nb_degree_of_freedom;

  b.clear();

  Array<Real>::matrix_iterator B_it = b.begin(tangent_size, bt_d_b_size);

  for (UInt e = 0; e < nb_element; ++e) {
    for (UInt q = 0; q < nb_quadrature_points; ++q) {
      Matrix<Real> & B  = *B_it;

      const Vector<Real> & Np  = *shape_Np;
      const Vector<Real> & Lpp = *shape_Lpp;
      const Vector<Real> & Mpp = *shape_Mpp;

      B(0,0)  =  Np(0);
      B(0,6)  =  Np(1);

      B(1,1)  =  Mpp(0);
      B(1,5)  =  Lpp(0);
      B(1,7)  =  Mpp(1);
      B(1,11) =  Lpp(1);

      B(2,2)  =  Mpp(0);
      B(2,4)  = -Lpp(0);
      B(2,8)  =  Mpp(1);
      B(2,10) = -Lpp(1);

      B(3,3)  =  Np(0);
      B(3,9)  =  Np(1);

      ++B_it;
      ++shape_Np;
      ++shape_Mpp;
      ++shape_Lpp;
    }
  }
}
/* -------------------------------------------------------------------------- */
template<> 
void StructuralMechanicsModel::transferBMatrixToSymVoigtBMatrix<_kirchhoff_shell>(Array<Real> & b, bool local) {
  MyFEEngineType & fem = getFEEngineClass<MyFEEngineType>();


  UInt nb_element                 = getFEEngine().getMesh().getNbElement(_kirchhoff_shell);
  UInt nb_nodes_per_element       = Mesh::getNbNodesPerElement(_kirchhoff_shell);
  UInt nb_quadrature_points       = getFEEngine().getNbQuadraturePoints(_kirchhoff_shell);



  Array<Real>::const_matrix_iterator shape_Np   = fem.getShapesDerivatives(_kirchhoff_shell, _not_ghost, 0).begin(2,nb_nodes_per_element);
  Array<Real>::const_matrix_iterator shape_Nx1p = fem.getShapesDerivatives(_kirchhoff_shell, _not_ghost, 1).begin(2,nb_nodes_per_element);
  Array<Real>::const_matrix_iterator shape_Nx2p = fem.getShapesDerivatives(_kirchhoff_shell, _not_ghost, 2).begin(2,nb_nodes_per_element);
  Array<Real>::const_matrix_iterator shape_Nx3p = fem.getShapesDerivatives(_kirchhoff_shell, _not_ghost, 3).begin(2,nb_nodes_per_element);
  Array<Real>::const_matrix_iterator shape_Ny1p = fem.getShapesDerivatives(_kirchhoff_shell, _not_ghost, 4).begin(2,nb_nodes_per_element);
  Array<Real>::const_matrix_iterator shape_Ny2p = fem.getShapesDerivatives(_kirchhoff_shell, _not_ghost, 5).begin(2,nb_nodes_per_element);
  Array<Real>::const_matrix_iterator shape_Ny3p = fem.getShapesDerivatives(_kirchhoff_shell, _not_ghost, 6).begin(2,nb_nodes_per_element);

  UInt tangent_size = getTangentStiffnessVoigtSize<_kirchhoff_shell>();
  UInt bt_d_b_size  = nb_nodes_per_element * nb_degree_of_freedom;

  b.clear();

  Array<Real>::matrix_iterator B_it = b.begin(tangent_size, bt_d_b_size);

  for (UInt e = 0; e < nb_element; ++e) {
    for (UInt q = 0; q < nb_quadrature_points; ++q) { 
      Matrix<Real> & B  = *B_it;

      const Matrix<Real> & Np   = *shape_Np;
      const Matrix<Real> & Nx1p = *shape_Nx1p; 
      const Matrix<Real> & Nx2p = *shape_Nx2p;
      const Matrix<Real> & Nx3p = *shape_Nx3p;
      const Matrix<Real> & Ny1p = *shape_Ny1p;
      const Matrix<Real> & Ny2p = *shape_Ny2p;
      const Matrix<Real> & Ny3p = *shape_Ny3p;

      B(0,0)  =  Np(0,0);
      B(0,6)  =  Np(0,1);  
      B(0,12) =  Np(0,2);

      B(1,1)  =  Np(1,0);
      B(1,7)  =  Np(1,1);
      B(1,13) =  Np(1,2);

      B(2,0)  =  Np(1,0);
      B(2,1)  =  Np(0,0);
      B(2,6)  =  Np(1,1);
      B(2,7)  =  Np(0,1);
      B(2,12) =  Np(1,2);
      B(2,13) =  Np(0,2);

      B(3,2)  =  Nx1p(0,0);
      B(3,3)  =  Nx2p(0,0);
      B(3,4)  =  Nx3p(0,0);
      B(3,8)  =  Nx1p(0,1);
      B(3,9)  =  Nx2p(0,1);
      B(3,10) =  Nx3p(0,1);
      B(3,14) =  Nx1p(0,2);
      B(3,15) =  Nx2p(0,2);
      B(3,16) =  Nx3p(0,2);

      B(4,2)  =  Ny1p(1,0);
      B(4,3)  =  Ny2p(1,0);
      B(4,4)  =  Ny3p(1,0);
      B(4,8)  =  Ny1p(1,1);
      B(4,9)  =  Ny2p(1,1);
      B(4,10) =  Ny3p(1,1);
      B(4,14) =  Ny1p(1,2);
      B(4,15) =  Ny2p(1,2);
      B(4,16) =  Ny3p(1,2);

      B(5,2)  =  Nx1p(1,0) + Ny1p(0,0);
      B(5,3)  =  Nx2p(1,0) + Ny2p(0,0);
      B(5,4)  =  Nx3p(1,0) + Ny3p(0,0);
      B(5,8)  =  Nx1p(1,1) + Ny1p(0,1);
      B(5,9)  =  Nx2p(1,1) + Ny2p(0,1);
      B(5,10) =  Nx3p(1,1) + Ny3p(0,1);
      B(5,14) =  Nx1p(1,2) + Ny1p(0,2);
      B(5,15) =  Nx2p(1,2) + Ny2p(0,2);
      B(5,16) =  Nx3p(1,2) + Ny3p(0,2);

      ++B_it;
      ++shape_Np;
      ++shape_Nx1p;
      ++shape_Nx2p;
      ++shape_Nx3p;
      ++shape_Ny1p;
      ++shape_Ny2p; 
      ++shape_Ny3p;
    }
  }
}


/* -------------------------------------------------------------------------- */

dumper::Field * StructuralMechanicsModel
::createNodalFieldBool(const std::string & field_name,
		       const std::string & group_name,
		       bool padding_flag) {
  
  std::map<std::string,Array<bool>* > uint_nodal_fields;
  uint_nodal_fields["blocked_dofs"             ] = blocked_dofs;

  dumper::Field * field = NULL;
  field = mesh.createNodalField(uint_nodal_fields[field_name],group_name);
  return field;

}

/* -------------------------------------------------------------------------- */


dumper::Field * StructuralMechanicsModel
::createNodalFieldReal(const std::string & field_name,
		       const std::string & group_name,
		       bool padding_flag) {

  UInt n;
  if(spatial_dimension == 2) {
    n = 2;
  } else n = 3;

  dumper::Field * field = NULL;

  if (field_name == "displacement"){
    field = mesh.createStridedNodalField(displacement_rotation,group_name,n,0,padding_flag);
  }
  if (field_name == "rotation"){
    field = mesh.createStridedNodalField(displacement_rotation,group_name,
					 nb_degree_of_freedom-n,n,padding_flag);
  }
  if (field_name == "force"){
    field = mesh.createStridedNodalField(force_momentum,group_name,n,0,padding_flag);
  }
  if (field_name == "momentum"){
    field = mesh.createStridedNodalField(force_momentum,group_name,
					 nb_degree_of_freedom-n,n,padding_flag);
  }
  if (field_name == "residual"){
    field = mesh.createNodalField(residual,group_name,padding_flag);
  }

  return field;

}

/* -------------------------------------------------------------------------- */

dumper::Field * StructuralMechanicsModel
::createElementalField(const std::string & field_name, 
		       const std::string & group_name,
		       bool padding_flag,
		       const ElementKind & kind){


  dumper::Field * field = NULL;

  if(field_name == "element_index_by_material") 
    field = mesh.createElementalField<UInt, Vector, dumper::ElementalField >(field_name,group_name,this->spatial_dimension,kind);

  return field;
}

/* -------------------------------------------------------------------------- */





__END_AKANTU__
