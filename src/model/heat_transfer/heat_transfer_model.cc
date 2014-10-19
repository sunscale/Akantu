/**
 * @file   heat_transfer_model.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Srinivasa Babu Ramisetti <srinivasa.ramisetti@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Rui Wang <rui.wang@epfl.ch>
 *
 * @date creation: Sun May 01 2011
 * @date last modification: Mon Sep 15 2014
 *
 * @brief  Implementation of HeatTransferModel class
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
#include "aka_common.hh"
#include "heat_transfer_model.hh"
#include "group_manager_inline_impl.cc"
#include "dumpable_inline_impl.hh"
#include "aka_math.hh"
#include "aka_common.hh"
#include "fe_engine_template.hh"
#include "mesh.hh"
#include "static_communicator.hh"
#include "parser.hh"
#include "generalized_trapezoidal.hh"

#ifdef AKANTU_USE_MUMPS
#include "solver_mumps.hh"
#endif

#ifdef AKANTU_USE_IOHELPER
#  include "dumper_paraview.hh"
#  include "dumper_elemental_field.hh"
#  include "dumper_element_partition.hh"
#  include "dumper_material_internal_field.hh"
#endif

/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__

const HeatTransferModelOptions default_heat_transfer_model_options(_explicit_lumped_capacity);

/* -------------------------------------------------------------------------- */
HeatTransferModel::HeatTransferModel(Mesh & mesh,
				     UInt dim,
				     const ID & id,
				     const MemoryID & memory_id) :
  Model(mesh, dim, id, memory_id), Parsable(_st_heat, id),
  integrator(new ForwardEuler()),
  stiffness_matrix(NULL),
  jacobian_matrix(NULL),
  temperature_gradient    ("temperature_gradient", id),
  temperature_on_qpoints  ("temperature_on_qpoints", id),
  conductivity_on_qpoints ("conductivity_on_qpoints", id),
  k_gradt_on_qpoints      ("k_gradt_on_qpoints", id),
  int_bt_k_gT             ("int_bt_k_gT", id),
  bt_k_gT                 ("bt_k_gT", id),
  conductivity(spatial_dimension, spatial_dimension),
  thermal_energy          ("thermal_energy", id) {
  AKANTU_DEBUG_IN();

  createSynchronizerRegistry(this);

  std::stringstream sstr; sstr << id << ":fem";
  registerFEEngineObject<MyFEEngineType>(sstr.str(), mesh,spatial_dimension);

  this->temperature= NULL;
  this->residual = NULL;
  this->blocked_dofs = NULL;

#ifdef AKANTU_USE_IOHELPER
  this->mesh.registerDumper<DumperParaview>("paraview_all", id, true);
  this->mesh.addDumpMesh(mesh, spatial_dimension, _not_ghost, _ek_regular);
#endif

  this->registerParam("conductivity"          , conductivity              , _pat_parsmod);
  this->registerParam("conductivity_variation", conductivity_variation, 0., _pat_parsmod);
  this->registerParam("temperature_reference" , T_ref                 , 0., _pat_parsmod);
  this->registerParam("capacity"              , capacity                  , _pat_parsmod);
  this->registerParam("density"               , density                   , _pat_parsmod);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::initModel() {
  getFEEngine().initShapeFunctions(_not_ghost);
  getFEEngine().initShapeFunctions(_ghost);
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::initParallel(MeshPartition * partition,
				     DataAccessor * data_accessor) {
  AKANTU_DEBUG_IN();

  if (data_accessor == NULL) data_accessor = this;
  Synchronizer & synch_parallel = createParallelSynch(partition,data_accessor);

  synch_registry->registerSynchronizer(synch_parallel, _gst_htm_capacity);
  synch_registry->registerSynchronizer(synch_parallel, _gst_htm_temperature);
  synch_registry->registerSynchronizer(synch_parallel, _gst_htm_gradient_temperature);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::initPBC() {
  AKANTU_DEBUG_IN();

  Model::initPBC();
  PBCSynchronizer * synch = new PBCSynchronizer(pbc_pair);

  synch_registry->registerSynchronizer(*synch, _gst_htm_capacity);
  synch_registry->registerSynchronizer(*synch, _gst_htm_temperature);
  changeLocalEquationNumberForPBC(pbc_pair,1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::initArrays() {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = getFEEngine().getMesh().getNbNodes();

  std::stringstream sstr_temp;      sstr_temp      << Model::id << ":temperature";
  std::stringstream sstr_temp_rate; sstr_temp_rate << Model::id << ":temperature_rate";
  std::stringstream sstr_inc;       sstr_inc       << Model::id << ":increment";
  std::stringstream sstr_ext_flx;   sstr_ext_flx   << Model::id << ":external_flux";
  std::stringstream sstr_residual;  sstr_residual  << Model::id << ":residual";
  std::stringstream sstr_lump;      sstr_lump      << Model::id << ":lumped";
  std::stringstream sstr_boun;      sstr_boun      << Model::id << ":blocked_dofs";

  temperature        = &(alloc<Real>(sstr_temp.str(),      nb_nodes, 1, REAL_INIT_VALUE));
  temperature_rate   = &(alloc<Real>(sstr_temp_rate.str(), nb_nodes, 1, REAL_INIT_VALUE));
  increment          = &(alloc<Real>(sstr_inc.str(),       nb_nodes, 1, REAL_INIT_VALUE));
  external_heat_rate = &(alloc<Real>(sstr_ext_flx.str(),   nb_nodes, 1, REAL_INIT_VALUE));
  residual           = &(alloc<Real>(sstr_residual.str(),  nb_nodes, 1, REAL_INIT_VALUE));
  capacity_lumped    = &(alloc<Real>(sstr_lump.str(),      nb_nodes, 1, REAL_INIT_VALUE));
  blocked_dofs           = &(alloc<bool>(sstr_boun.str(),      nb_nodes, 1, false));

  Mesh::ConnectivityTypeList::const_iterator it;

  /* -------------------------------------------------------------------------- */
  // byelementtype vectors
  getFEEngine().getMesh().initElementTypeMapArray(temperature_on_qpoints,
					     1,
					     spatial_dimension);

  getFEEngine().getMesh().initElementTypeMapArray(temperature_gradient,
					     spatial_dimension,
					     spatial_dimension);

  getFEEngine().getMesh().initElementTypeMapArray(conductivity_on_qpoints,
					     spatial_dimension*spatial_dimension,
					     spatial_dimension);

  getFEEngine().getMesh().initElementTypeMapArray(k_gradt_on_qpoints,
					     spatial_dimension,
					     spatial_dimension);

  getFEEngine().getMesh().initElementTypeMapArray(bt_k_gT,
					     1,
					     spatial_dimension,true);

  getFEEngine().getMesh().initElementTypeMapArray(int_bt_k_gT,
					     1,
					     spatial_dimension,true);

  getFEEngine().getMesh().initElementTypeMapArray(thermal_energy,
					     1,
					     spatial_dimension);


  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;

    const Mesh::ConnectivityTypeList & type_list =
      getFEEngine().getMesh().getConnectivityTypeList(gt);

    for(it = type_list.begin(); it != type_list.end(); ++it) {
      if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;
      UInt nb_element = getFEEngine().getMesh().getNbElement(*it, gt);
      UInt nb_quad_points = this->getFEEngine().getNbQuadraturePoints(*it, gt) * nb_element;

      temperature_on_qpoints(*it, gt).resize(nb_quad_points);
      temperature_on_qpoints(*it, gt).clear();

      temperature_gradient(*it, gt).resize(nb_quad_points);
      temperature_gradient(*it, gt).clear();

      conductivity_on_qpoints(*it, gt).resize(nb_quad_points);
      conductivity_on_qpoints(*it, gt).clear();

      k_gradt_on_qpoints(*it, gt).resize(nb_quad_points);
      k_gradt_on_qpoints(*it, gt).clear();

      bt_k_gT(*it, gt).resize(nb_quad_points);
      bt_k_gT(*it, gt).clear();

      int_bt_k_gT(*it, gt).resize(nb_element);
      int_bt_k_gT(*it, gt).clear();

      thermal_energy(*it, gt).resize(nb_element);
      thermal_energy(*it, gt).clear();
    }
  }

  /* -------------------------------------------------------------------------- */
  dof_synchronizer = new DOFSynchronizer(getFEEngine().getMesh(),1);
  dof_synchronizer->initLocalDOFEquationNumbers();
  dof_synchronizer->initGlobalDOFEquationNumbers();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::initSolver(__attribute__((unused)) SolverOptions & options) {
#if !defined(AKANTU_USE_MUMPS) // or other solver in the future \todo add AKANTU_HAS_SOLVER in CMake
  AKANTU_DEBUG_ERROR("You should at least activate one solver.");
#else
  UInt nb_global_nodes = mesh.getNbGlobalNodes();

  delete jacobian_matrix;
  std::stringstream sstr; sstr << Memory::id << ":jacobian_matrix";
  jacobian_matrix = new SparseMatrix(nb_global_nodes, _symmetric, sstr.str(), memory_id);

  jacobian_matrix->buildProfile(mesh, *dof_synchronizer, 1);

  delete stiffness_matrix;
  std::stringstream sstr_sti; sstr_sti << Memory::id << ":stiffness_matrix";
  stiffness_matrix = new SparseMatrix(*jacobian_matrix, sstr_sti.str(), memory_id);

#ifdef AKANTU_USE_MUMPS
  std::stringstream sstr_solv; sstr_solv << Memory::id << ":solver";
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
void HeatTransferModel::initImplicit(SolverOptions & solver_options) {
  method = _static;
  initSolver(solver_options);
}

/* -------------------------------------------------------------------------- */

HeatTransferModel::~HeatTransferModel()
{
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::assembleCapacityLumped(const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  FEEngine & fem = getFEEngine();

  const Mesh::ConnectivityTypeList & type_list
    = fem.getMesh().getConnectivityTypeList(ghost_type);

  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it)
  {
    if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;

    UInt nb_element = getFEEngine().getMesh().getNbElement(*it,ghost_type);
    UInt nb_quadrature_points = getFEEngine().getNbQuadraturePoints(*it, ghost_type);

    Array<Real> rho_1 (nb_element * nb_quadrature_points,1, capacity * density);
    fem.assembleFieldLumped(rho_1,1,*capacity_lumped,
			    dof_synchronizer->getLocalDOFEquationNumbers(),
			    *it, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::assembleCapacityLumped() {
  AKANTU_DEBUG_IN();

  capacity_lumped->clear();

  assembleCapacityLumped(_not_ghost);
  assembleCapacityLumped(_ghost);

  getSynchronizerRegistry().synchronize(_gst_htm_capacity);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::updateResidual() {
  AKANTU_DEBUG_IN();
  /// @f$ r = q_{ext} - q_{int} - C \dot T @f$

  // start synchronization
  synch_registry->asynchronousSynchronize(_gst_htm_temperature);
  // finalize communications
  synch_registry->waitEndSynchronize(_gst_htm_temperature);

  //clear the array
  /// first @f$ r = q_{ext} @f$
  //  residual->clear();
  residual->copy(*external_heat_rate);

  /// then @f$ r -= q_{int} @f$
  // update the not ghost ones

  updateResidual(_not_ghost);
  // update for the received ghosts
  updateResidual(_ghost);

/*  if (method == _explicit_lumped_capacity) {
    this->solveExplicitLumped();
  }*/


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::assembleStiffnessMatrix() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Assemble the new stiffness matrix.");

  stiffness_matrix->clear();

  switch(mesh.getSpatialDimension()) {
    case 1: this->assembleStiffnessMatrix<1>(_not_ghost); break;
    case 2: this->assembleStiffnessMatrix<2>(_not_ghost); break;
    case 3: this->assembleStiffnessMatrix<3>(_not_ghost); break;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void HeatTransferModel::assembleStiffnessMatrix(const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  Mesh & mesh = this->getFEEngine().getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type);
  for(; it != last_type; ++it) {
    this->assembleStiffnessMatrix<dim>(*it, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}
  
/* -------------------------------------------------------------------------- */
template <UInt dim>
void HeatTransferModel::assembleStiffnessMatrix(const ElementType & type, const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  SparseMatrix & K = *stiffness_matrix;
  const Array<Real> & shapes_derivatives = this->getFEEngine().getShapesDerivatives(type, ghost_type);

  UInt nb_element                 = mesh.getNbElement(type, ghost_type);
  UInt nb_nodes_per_element       = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points       = getFEEngine().getNbQuadraturePoints(type, ghost_type);

  /// compute @f$\mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  UInt bt_d_b_size = nb_nodes_per_element;

  Array<Real> * bt_d_b = new Array<Real>(nb_element * nb_quadrature_points,
					                               bt_d_b_size * bt_d_b_size,
					                               "B^t*D*B");

  Matrix<Real> Bt_D(nb_nodes_per_element, dim);

  Array<Real>::const_iterator< Matrix<Real> > shapes_derivatives_it = shapes_derivatives.begin(dim, nb_nodes_per_element);

  Array<Real>::iterator< Matrix<Real> > Bt_D_B_it = bt_d_b->begin(bt_d_b_size, bt_d_b_size);

  this->computeConductivityOnQuadPoints(ghost_type);
  Array<Real>::iterator< Matrix<Real> > D_it = conductivity_on_qpoints(type, ghost_type).begin(dim, dim);
  Array<Real>::iterator< Matrix<Real> > D_end = conductivity_on_qpoints(type, ghost_type).end(dim, dim);

  for (; D_it != D_end; ++D_it, ++Bt_D_B_it, ++shapes_derivatives_it) {
    Matrix<Real> & D = *D_it;
    const Matrix<Real> & B = *shapes_derivatives_it;
    Matrix<Real> & Bt_D_B = *Bt_D_B_it;

    Bt_D.mul<true, false>(B, D);
    Bt_D_B.mul<false, false>(Bt_D, B);
  }

  /// compute @f$ k_e = \int_e \mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  Array<Real> * K_e = new Array<Real>(nb_element,
				      bt_d_b_size * bt_d_b_size,
				      "K_e");

  this->getFEEngine().integrate(*bt_d_b, *K_e,
                            bt_d_b_size * bt_d_b_size,
                            type, ghost_type);

  delete bt_d_b;

  this->getFEEngine().assembleMatrix(*K_e, K, 1, type, ghost_type);
  delete K_e;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::solveStatic() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Solving Ku = f");
  AKANTU_DEBUG_ASSERT(stiffness_matrix != NULL,
                      "You should first initialize the implicit solver and assemble the stiffness matrix");

  UInt nb_nodes = temperature->getSize();
  UInt nb_degree_of_freedom = temperature->getNbComponent() * nb_nodes;

  jacobian_matrix->copyContent(*stiffness_matrix);
  jacobian_matrix->applyBoundary(*blocked_dofs);

  increment->clear();

  solver->setRHS(*residual);
  solver->factorize();
  solver->solve(*increment);

  Real * increment_val = increment->storage();
  Real * temperature_val = temperature->storage();
  bool * blocked_dofs_val = blocked_dofs->storage();

  for (UInt j = 0; j < nb_degree_of_freedom;
       ++j, ++temperature_val, ++increment_val, ++blocked_dofs_val) {
    if (!(*blocked_dofs_val)) {
      *temperature_val += *increment_val;
    }
    else {
      *increment_val = 0.0;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::computeConductivityOnQuadPoints(const GhostType & ghost_type) {
  const Mesh::ConnectivityTypeList & type_list =
    this->getFEEngine().getMesh().getConnectivityTypeList(ghost_type);
  Mesh::ConnectivityTypeList::const_iterator it;

  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;

    Array<Real> & temperature_interpolated = temperature_on_qpoints(*it, ghost_type);

    //compute the temperature on quadrature points
    this->getFEEngine().interpolateOnQuadraturePoints(*temperature,
						 temperature_interpolated,
						 1 ,*it,ghost_type);

    Array<Real>::matrix_iterator C_it =
      conductivity_on_qpoints(*it, ghost_type).begin(spatial_dimension, spatial_dimension);
    Array<Real>::matrix_iterator C_end =
      conductivity_on_qpoints(*it, ghost_type).end(spatial_dimension, spatial_dimension);

    Array<Real>::iterator<Real>  T_it = temperature_interpolated.begin();

    for (;C_it != C_end; ++C_it, ++T_it) {
      Matrix<Real> & C = *C_it;
      Real & T = *T_it;
      C = conductivity;

      Matrix<Real> variation(spatial_dimension, spatial_dimension, conductivity_variation * (T - T_ref));
      C += conductivity_variation;
    }
  }
  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
void HeatTransferModel::computeKgradT(const GhostType & ghost_type) {

  computeConductivityOnQuadPoints(ghost_type);

  const Mesh::ConnectivityTypeList & type_list =
    this->getFEEngine().getMesh().getConnectivityTypeList(ghost_type);
  Mesh::ConnectivityTypeList::const_iterator it;

  for(it = type_list.begin(); it != type_list.end(); ++it) {
    const ElementType & type = *it;
    if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;

    Array<Real> & gradient = temperature_gradient(*it, ghost_type);
    this->getFEEngine().gradientOnQuadraturePoints(*temperature,
					      gradient,
					      1 ,*it, ghost_type);

    Array<Real>::matrix_iterator C_it =
      conductivity_on_qpoints(*it, ghost_type).begin(spatial_dimension, spatial_dimension);

    Array<Real>::vector_iterator BT_it = gradient.begin(spatial_dimension);

    Array<Real>::vector_iterator k_BT_it =
      k_gradt_on_qpoints(type, ghost_type).begin(spatial_dimension);
    Array<Real>::vector_iterator k_BT_end =
      k_gradt_on_qpoints(type, ghost_type).end(spatial_dimension);

    for (;k_BT_it != k_BT_end; ++k_BT_it, ++BT_it, ++C_it) {
      Vector<Real> & k_BT = *k_BT_it;
      Vector<Real> & BT = *BT_it;
      Matrix<Real> & C = *C_it;

      k_BT.mul<false>(C, BT);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::updateResidual(const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  const Mesh::ConnectivityTypeList & type_list =
    this->getFEEngine().getMesh().getConnectivityTypeList(ghost_type);
  Mesh::ConnectivityTypeList::const_iterator it;

  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(Mesh::getSpatialDimension(*it) != spatial_dimension) continue;

    Array<Real> & shapes_derivatives =
      const_cast<Array<Real> &>(getFEEngine().getShapesDerivatives(*it,ghost_type));

    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(*it);

    // compute k \grad T
    computeKgradT(ghost_type);

    Array<Real>::vector_iterator k_BT_it =
      k_gradt_on_qpoints(*it,ghost_type).begin(spatial_dimension);

    Array<Real>::matrix_iterator B_it =
      shapes_derivatives.begin(spatial_dimension, nb_nodes_per_element);

    Array<Real>::vector_iterator Bt_k_BT_it =
      bt_k_gT(*it,ghost_type).begin(nb_nodes_per_element);

    Array<Real>::vector_iterator Bt_k_BT_end =
      bt_k_gT(*it,ghost_type).end(nb_nodes_per_element);


    for (;Bt_k_BT_it != Bt_k_BT_end; ++Bt_k_BT_it, ++B_it, ++k_BT_it) {
      Vector<Real> & k_BT = *k_BT_it;
      Vector<Real> & Bt_k_BT = *Bt_k_BT_it;
      Matrix<Real> & B = *B_it;

      Bt_k_BT.mul<true>(B, k_BT);
    }

    this->getFEEngine().integrate(bt_k_gT(*it,ghost_type),
			     int_bt_k_gT(*it,ghost_type),
			     nb_nodes_per_element, *it,ghost_type);

    this->getFEEngine().assembleArray(int_bt_k_gT(*it,ghost_type), *residual,
				  dof_synchronizer->getLocalDOFEquationNumbers(),
				 1, *it,ghost_type, empty_filter, -1);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::solveExplicitLumped() {
  AKANTU_DEBUG_IN();
  
  /// finally @f$ r -= C \dot T @f$
  // lumped C
  UInt nb_nodes            = temperature_rate->getSize();
  UInt nb_degree_of_freedom = temperature_rate->getNbComponent();

  Real * capacity_val  = capacity_lumped->storage();
  Real * temp_rate_val = temperature_rate->storage();
  Real * res_val       = residual->storage();
  bool * blocked_dofs_val  = blocked_dofs->storage();

  for (UInt n = 0; n < nb_nodes * nb_degree_of_freedom; ++n) {
    if(!(*blocked_dofs_val)) {
      *res_val -= *capacity_val * *temp_rate_val;
    }
    blocked_dofs_val++;
    res_val++;
    capacity_val++;
    temp_rate_val++;
  }

#ifndef AKANTU_NDEBUG
  getSynchronizerRegistry().synchronize(akantu::_gst_htm_gradient_temperature);
#endif

  capacity_val      = capacity_lumped->storage();
  res_val           = residual->storage();
  blocked_dofs_val      = blocked_dofs->storage();
  Real * inc           = increment->storage();

  for (UInt n = 0; n < nb_nodes * nb_degree_of_freedom; ++n) {
    if(!(*blocked_dofs_val)) {
      *inc = (*res_val / *capacity_val);
    }
    res_val++;
    blocked_dofs_val++;
    inc++;
    capacity_val++;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::explicitPred() {
  AKANTU_DEBUG_IN();

  integrator->integrationSchemePred(time_step,
				    *temperature,
				    *temperature_rate,
				    *blocked_dofs);

  UInt nb_nodes = temperature->getSize();
  UInt nb_degree_of_freedom = temperature->getNbComponent();

  Real * temp = temperature->storage();
  for (UInt n = 0; n < nb_nodes * nb_degree_of_freedom; ++n, ++temp)
    if(*temp < 0.) *temp = 0.;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::explicitCorr() {
  AKANTU_DEBUG_IN();

  integrator->integrationSchemeCorrTempRate(time_step,
					    *temperature,
					    *temperature_rate,
					    *blocked_dofs,
					    *increment);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Real HeatTransferModel::getStableTimeStep()
{
  AKANTU_DEBUG_IN();

  Real el_size;
  Real min_el_size = std::numeric_limits<Real>::max();
  Real conductivitymax = conductivity(0, 0);

  //get the biggest parameter from k11 until k33//
  for(UInt i = 0; i < spatial_dimension; i++)
    for(UInt j = 0; j < spatial_dimension; j++)
      conductivitymax = std::max(conductivity(i, j), conductivitymax);


  const Mesh::ConnectivityTypeList & type_list =
    getFEEngine().getMesh().getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(getFEEngine().getMesh().getSpatialDimension(*it) != spatial_dimension)
      continue;

    UInt nb_nodes_per_element = getFEEngine().getMesh().getNbNodesPerElement(*it);

    Array<Real> coord(0, nb_nodes_per_element*spatial_dimension);
    FEEngine::extractNodalToElementField(getFEEngine().getMesh(), getFEEngine().getMesh().getNodes(),
				    coord, *it, _not_ghost);

    Array<Real>::matrix_iterator el_coord = coord.begin(spatial_dimension, nb_nodes_per_element);
    UInt nb_element           = getFEEngine().getMesh().getNbElement(*it);

    for (UInt el = 0; el < nb_element; ++el, ++el_coord) {
	el_size    = getFEEngine().getElementInradius(*el_coord, *it);
	min_el_size = std::min(min_el_size, el_size);
    }

    AKANTU_DEBUG_INFO("The minimum element size : " << min_el_size
		      << " and the max conductivity is : "
		      << conductivitymax);
  }

  Real min_dt = 2 * min_el_size * min_el_size * density
    * capacity / conductivitymax;

  StaticCommunicator::getStaticCommunicator().allReduce(&min_dt, 1, _so_min);

  AKANTU_DEBUG_OUT();

  return min_dt;
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::readMaterials() {
  std::pair<Parser::const_section_iterator, Parser::const_section_iterator>
    sub_sect = this->parser->getSubSections(_st_heat);

  Parser::const_section_iterator it = sub_sect.first;
  const ParserSection & section = *it;

  this->parseSection(section);
}

/* -------------------------------------------------------------------------- */

void HeatTransferModel::initFull(const ModelOptions & options){
  Model::initFull(options);

  readMaterials();

  const HeatTransferModelOptions & my_options = dynamic_cast<const HeatTransferModelOptions &>(options);

  //initialize the vectors
  initArrays();
  temperature->clear();
  temperature_rate->clear();
  external_heat_rate->clear();

  method = my_options.analysis_method;

  if (method == _static) {
    initImplicit();
  }
}

/* -------------------------------------------------------------------------- */
void HeatTransferModel::initFEEngineBoundary(bool create_surface) {

  if(create_surface)
    MeshUtils::buildFacets(getFEEngine().getMesh());

  FEEngine & fem_boundary = getFEEngineBoundary();
  fem_boundary.initShapeFunctions();
  fem_boundary.computeNormalsOnControlPoints();
}

/* -------------------------------------------------------------------------- */
Real HeatTransferModel::computeThermalEnergyByNode() {
  AKANTU_DEBUG_IN();

  Real ethermal = 0.;

  Array<Real>::vector_iterator heat_rate_it =
    residual->begin(residual->getNbComponent());

  Array<Real>::vector_iterator heat_rate_end =
    residual->end(residual->getNbComponent());


  UInt n = 0;
  for(;heat_rate_it != heat_rate_end; ++heat_rate_it, ++n) {
    Real heat = 0;
    bool is_local_node = mesh.isLocalOrMasterNode(n);
    bool is_not_pbc_slave_node = !isPBCSlaveNode(n);
    bool count_node = is_local_node && is_not_pbc_slave_node;

    Vector<Real> & heat_rate = *heat_rate_it;
    for (UInt i = 0; i < heat_rate.size(); ++i) {
      if (count_node)
	heat += heat_rate[i] * time_step;
    }
    ethermal += heat;
  }

  StaticCommunicator::getStaticCommunicator().allReduce(&ethermal, 1, _so_sum);

  AKANTU_DEBUG_OUT();
  return ethermal;
}

/* -------------------------------------------------------------------------- */
template<class iterator>
void HeatTransferModel::getThermalEnergy(iterator Eth,
					 Array<Real>::const_iterator<Real> T_it,
					 Array<Real>::const_iterator<Real> T_end) const {
  for(;T_it != T_end; ++T_it, ++Eth) {
    *Eth = capacity * density * *T_it;
  }
}

/* -------------------------------------------------------------------------- */
Real HeatTransferModel::getThermalEnergy(const ElementType & type, UInt index) {
  AKANTU_DEBUG_IN();

  UInt nb_quadrature_points = getFEEngine().getNbQuadraturePoints(type);
  Vector<Real> Eth_on_quarature_points(nb_quadrature_points);

  Array<Real>::iterator<Real> T_it  = this->temperature_on_qpoints(type).begin();
  T_it  += index * nb_quadrature_points;

  Array<Real>::iterator<Real> T_end = T_it + nb_quadrature_points;

  getThermalEnergy(Eth_on_quarature_points.storage(), T_it, T_end);

  return getFEEngine().integrate(Eth_on_quarature_points, type, index);
}

/* -------------------------------------------------------------------------- */
Real HeatTransferModel::getThermalEnergy() {
  Real Eth = 0;

  Mesh & mesh = getFEEngine().getMesh();
  Mesh::type_iterator it = mesh.firstType(spatial_dimension);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension);
  for(; it != last_type; ++it) {
    UInt nb_element = getFEEngine().getMesh().getNbElement(*it, _not_ghost);
    UInt nb_quadrature_points = getFEEngine().getNbQuadraturePoints(*it, _not_ghost);
    Array<Real> Eth_per_quad(nb_element * nb_quadrature_points, 1);

    Array<Real>::iterator<Real> T_it  = this->temperature_on_qpoints(*it).begin();
    Array<Real>::iterator<Real> T_end = this->temperature_on_qpoints(*it).end();
    getThermalEnergy(Eth_per_quad.begin(), T_it, T_end);

    Eth += getFEEngine().integrate(Eth_per_quad, *it);
  }

  return Eth;
}

/* -------------------------------------------------------------------------- */
Real HeatTransferModel::getEnergy(const std::string & id) {
  AKANTU_DEBUG_IN();
  Real energy = 0;

  if("thermal") energy = getThermalEnergy();

  // reduction sum over all processors
  StaticCommunicator::getStaticCommunicator().allReduce(&energy, 1, _so_sum);

  AKANTU_DEBUG_OUT();
  return energy;
}

/* -------------------------------------------------------------------------- */
Real HeatTransferModel::getEnergy(const std::string & energy_id, const ElementType & type, UInt index) {
  AKANTU_DEBUG_IN();

  Real energy = 0.;

  if("thermal") energy = getThermalEnergy(type, index);

  AKANTU_DEBUG_OUT();
  return energy;
}

/* -------------------------------------------------------------------------- */

dumper::Field * HeatTransferModel::createNodalFieldBool(const std::string & field_name,
							const std::string & group_name,
							bool padding_flag) {
  
  std::map<std::string,Array<bool>* > uint_nodal_fields;
  uint_nodal_fields["blocked_dofs"             ] = blocked_dofs;

  dumper::Field * field = NULL;
  field = mesh.createNodalField(uint_nodal_fields[field_name],group_name);
  return field;

}
/* -------------------------------------------------------------------------- */

dumper::Field * HeatTransferModel::createNodalFieldReal(const std::string & field_name,
							const std::string & group_name,
							bool padding_flag) {
  
  std::map<std::string,Array<Real>* > real_nodal_fields;
  real_nodal_fields["temperature"              ] = temperature;
  real_nodal_fields["temperature_rate"         ] = temperature_rate;
  real_nodal_fields["external_heat_rate"       ] = external_heat_rate;
  real_nodal_fields["residual"                 ] = residual;
  real_nodal_fields["capacity_lumped"          ] = capacity_lumped;

  dumper::Field * field = 
    mesh.createNodalField(real_nodal_fields[field_name],group_name);
  
  return field;
}
/* -------------------------------------------------------------------------- */

dumper::Field * HeatTransferModel
::createElementalField(const std::string & field_name, 
		       const std::string & group_name,
		       bool padding_flag,
		       const ElementKind & element_kind){


  dumper::Field * field = NULL;

  if(field_name == "partitions") 
    field = mesh.createElementalField<UInt, dumper::ElementPartitionField>(mesh.getConnectivities(),group_name,this->spatial_dimension,element_kind);
  else if(field_name == "temperature_gradient"){
    ElementTypeMap<UInt> nb_data_per_elem = this->mesh.getNbDataPerElem(temperature_gradient,element_kind);

    field = 
      mesh.createElementalField<Real, 
				dumper::InternalMaterialField>(temperature_gradient,
							       group_name,
							       this->spatial_dimension,
							       element_kind,
							       nb_data_per_elem);
  }

  return field;
}

/* -------------------------------------------------------------------------- */
 
 


__END_AKANTU__
