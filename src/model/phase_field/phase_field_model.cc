/**
 * @file   phase_field_model.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Wed Aug 01 2018
 * @date last modification: Wed Aug 01 2018
 *
 * @brief  Implementation of PhaseFieldModel class
 *
 * @section LICENSE
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
#include "phase_field_model.hh"
#include "dumpable_inline_impl.hh"
#include "element_synchronizer.hh"
#include "fe_engine_template.hh"
#include "group_manager_inline_impl.cc"
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

  // functor can be created 
  
}
  

/* -------------------------------------------------------------------------- */
PhaseFieldModel::PhaseFieldModel(Mesh & mesh, UInt dim, const ID & id,
				 const MemoryID & memory_id)
  : Model(mesh, ModelType::_phase_field_model, dim, id, memory_id),
    damage_gradient("damage_gradient", id),
    damage_on_qpoints("damage_on_qpoints", id),
    strain_history_on_qpoints("strain_history_on_qpoints", id),
    displacement_on_qpoints("displacement_on_qpoints", id){

 AKANTU_DEBUG_IN();

 this->initDOFManager();

 this->registerDataAccessor(*this);

 if (this->mesh.isDistributed()) {
   auto & synchronizer = this->mesh.getElementSynchronizer();
   this->registerSynchronizer(synchronizer, _gst_pfm_damage);
   this->registerSynchronizer(synchronizer, _gst_pfm_gradient_damage);
 }

 registerFEEngineObject<FEEngineType>(id + ":fem", mesh, spatial_dimension);

#ifdef AKANTU_USE_IOHELPER
 this->mesh.registerDumper<DumperParaview>("phase_field", id, true);
 this->mesh.addDumpMesh(mesh, spatial_dimension, _not_ghost, _ek_regular);
#endif // AKANTU_USE_IOHELPER

 this->registerParam("l0", l_0, 0., _pat_parsmod, "length scale");
 this->registerParam("gc", g_c, _pat_parsmod, "critical local fracture energy density");
 this->registerParam("lambda", lame_lambda, _pat_parsmod, "lame parameter lambda");
 this->registerParam("mu", lame_mu, _pat_parsmod, "lame paramter mu");

 
 AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void PhaseFieldModel::initModel() {
  auto & fem = this->getFEEngine();
  fem.initShapeFunctions(_not_ghost);
  fem.initShapeFunctions(_ghost);

  damage_on_qpoints.initialize(fem, _nb_component = 1);
  damage_gradient.initialize(fem, _nb_component = spatial_dimension);
  displacement_on_qpoints.initialize(fem, _nb_component = spatial_dimension);
  strain_history_on_qpoints.initialize(fem, _nb_component = 1);
  
}


/* -------------------------------------------------------------------------- */
void PhaseFieldModel::assembleMatrix(const ID & matrix_id) {
  if (matrix_id == "K") {
    this->assembleDamageMatrix();
  } else if (matrix_id == "M" and need_to_reassemble_capacity) {
    this->assembleDamageGradMatrix();
  }
}


/* -------------------------------------------------------------------------- */
void PhaseFieldModel::initSolver(TimeStepSolverType time_step_solver_type,
                                 NonLinearSolverType) {
  DOFManager & dof_manager = this->getDOFManager();

  this->allocNodalField(this->damage, "damage");
  this->allocNodalField(this->external_force, "external_force");
  this->allocNodalField(this->internal_force, "internal_force");
  this->allocNodalField(this->blocked_dofs, "blocked_dofs");

  if (!dof_manager.hasDOFs("damage")) {
    dof_manager.registerDOFs("damage", *this->damage, _dst_nodal);
    dof_manager.registerBlockedDOFs("damage", *this->blocked_dofs);
  }

  /*if (time_step_solver_type == _tsst_dynamic ||
      time_step_solver_type == _tsst_dynamic_lumped) {
    this->allocNodalField(this->damage_rate, "damage_rate");

    if (!dof_manager.hasDOFsDerivatives("damage", 1)) {
      dof_manager.registerDOFsDerivative("damage", 1,
                                         *this->damage_rate);
    }
    }*/
}



/* -------------------------------------------------------------------------- */
void PhaseFieldModel::assembleDamageMatrix() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Assemble the new damage matrix");

  this->computeStrainHistoryOnQuadPoints(_not_ghost);
  
  if (!this->getDOFManager().hasMatrix("K")) {
    this->getDOFManager().getNewMatrix("K", getMatrixType("K"));
  }
  this->getDOFManager().clearMatrix("K");

  switch (mesh.getSpatialDimension()) {
  case 1:
    this->assembleDamageMatrix<1>(_not_ghost);
    break;
  case 2:
    this->assembleDamageMatrix<2>(_not_ghost);
    break;
  case 3:
    this->assembleDamageMatrix<3>(_not_ghost);
    break;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::assembleDamageGradMatrix() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Assemble the new damage gradient matrix");

  if (!this->getDOFManager().hasMatrix("K")) {
    this->getDOFManager().getNewMatrix("K", getMatrixType("K"));
  }
  this->getDOFManager().clearMatrix("K");

  switch (mesh.getSpatialDimension()) {
  case 1:
    this->assembleDamageGradMatrix<1>(_not_ghost);
    break;
  case 2:
    this->assembleDamageGradMatrix<2>(_not_ghost);
    break;
  case 3:
    this->assembleDamageGradMatrix<3>(_not_ghost);
    break;
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
std::tuple<ID, TimeStepSolverType>
PhaseFieldModel::getDefaultSolverID(const AnalysisMethod & method) {
  switch (method) {
  case _explicit_lumped_mass: {
    return std::make_tuple("explicit_lumped", _tsst_dynamic_lumped);
  }
  case _static: {
    return std::make_tuple("static", _tsst_static);
  }
  case _implicit_dynamic: {
    return std::make_tuple("implicit", _tsst_dynamic);
  }
  default:
    return std::make_tuple("unknown", _tsst_not_defined);
  }
}

/* -------------------------------------------------------------------------- */
ModelSolverOptions PhaseFieldModel::getDefaultSolverOptions(
    const TimeStepSolverType & type) const {
  ModelSolverOptions options;

  switch (type) {
  case _tsst_dynamic_lumped: {
    options.non_linear_solver_type = _nls_lumped;
    options.integration_scheme_type["damage"] = _ist_forward_euler;
    options.solution_type["damage"] = IntegrationScheme::_temperature_rate;
    break;
  }
  case _tsst_static: {
    options.non_linear_solver_type = _nls_newton_raphson;
    options.integration_scheme_type["damage"] = _ist_pseudo_time;
    options.solution_type["damage"] = IntegrationScheme::_not_defined;
    break;
  }
  case _tsst_dynamic: {
    if (this->method == _explicit_consistent_mass) {
      options.non_linear_solver_type = _nls_newton_raphson;
      options.integration_scheme_type["damage"] = _ist_forward_euler;
      options.solution_type["damage"] =
          IntegrationScheme::_temperature_rate;
    } else {
      options.non_linear_solver_type = _nls_newton_raphson;
      options.integration_scheme_type["damage"] = _ist_backward_euler;
      options.solution_type["damage"] = IntegrationScheme::_temperature;
    }
    break;
  }
  default:
    AKANTU_EXCEPTION(type << " is not a valid time step solver type");
  }

  return options;
}

  
/* -------------------------------------------------------------------------- */
template<UInt dim>
void PhaseFieldModel::assembleDamageMatrix(const GhostType & ghost_type) {

  AKANTU_DEBUG_IN();

 
  auto & fem = this->getFEEngine();

  for (auto && type : mesh.elementTypes(spatial_dimension, ghost_type)) {
    auto nb_element = mesh.getNbElement(type, ghost_type);
    auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    auto nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);

    auto nt_d_n = std::make_unique<Array<Real>>(
	nb_element * nb_quadrature_points,
	nb_nodes_per_element * nb_nodes_per_element, "N^t*D*N");

    // conductivity_on_qpoints = gc*l0 = scalar
    fem.computeNtDN(conductivity_on_qpoints(type, ghost_type), *nt_d_n, type,
		    ghost_type);

    /// compute @f$ K_d = \int_e \mathbf{N}^t * \mathbf{D} *
    /// \mathbf{N}@f$
    auto K_d = std::make_unique<Array<Real>>(
	nb_element, nb_nodes_per_element * nb_nodes_per_element, "K_d");

    fem.integrate(*nt_d_n, *K_d, nb_nodes_per_element * nb_nodes_per_element,
		  type, ghost_type);

    this->getDOFManager().assembleElementalMatricesToMatrix(
       "K", "phasefield", *K_d, type, ghost_type, _symmetric);
  }

  AKANTU_DEBUG_OUT();
}

  
/* -------------------------------------------------------------------------- */
template<UInt dim>
void PhaseFieldModel::assembleDamageGradMatrix(const GhostType & ghost_type) {

  AKANTU_DEBUG_IN();

 
  auto & fem = this->getFEEngine();

  for (auto && type : mesh.elementTypes(spatial_dimension, ghost_type)) {
    auto nb_element = mesh.getNbElement(type, ghost_type);
    auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    auto nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);

    auto bt_d_b = std::make_unique<Array<Real>>(
       nb_element * nb_quadrature_points,
       nb_nodes_per_element * nb_nodes_per_element, "B^t*D*B");

    // conductivity_on_qpoints = gc*l0 = scalar
    fem.computeBtDB(conductivity_on_qpoints(type, ghost_type), *bt_d_b, type,
		    ghost_type);

    /// compute @f$ K_{\grad d} = \int_e \mathbf{B}^t * \mathbf{W} *
    /// \mathbf{B}@f$
    auto K_b = std::make_unique<Array<Real>>(
	nb_element, nb_nodes_per_element * nb_nodes_per_element, "K_b");

    fem.integrate(*bt_d_b, *K_b, nb_nodes_per_element * nb_nodes_per_element,
		  type, ghost_type);

    this->getDOFManager().assembleElementalMatricesToMatrix(
       "K", "phasefield", *K_b, type, ghost_type, _symmetric);
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void PhaseFieldModel::computeStrainHistoryOnQuadPoints(
     const GhostType & ghost_type) {

  for (auto & type : mesh.elementType(spatial_dimension, ghost_type)) {
    auto & displacement_interpolated = displacement_on_qpoints(type, ghost_type);

    // compute the strain on quadrature points
    this->getFEEngine().interpolateOnIntegrationPoints(
	    *displacement, displacement_interpolated, , type, ghost_type);

    auto & strain_history = strain_history_on_qpoints(type, ghost_type);
    for (auto  && tuple :
	   zip(make_view(strain_history, spatial_dimension, spatial_dimension),
	       displacement_interpolated)) {

      // compute strain on quad from displacement_interpolated
      // commpute matrix sigma_plus and sigma_minus using
      // lame_lambda and lame_mu

      // compute phi_plus
      // updated strain_history 
      
     
    }
  }

}

/* -------------------------------------------------------------------------- */
void PhaseFieldModel::assembleResidual() {

  AKANTU_DEBUG_IN();

  this->assembleInternalForces();

  this->getDOFManager().assembleToResidual("damage",
					   *this->external_force, 1);
  this->getDOFManager().assembleToResidual("damage",
					   *this->internal_force, 1);


  AKANTU_DEBUG_OUT();
}



/* -------------------------------------------------------------------------- */
void PhaseFieldModel::assembleInternalForces() {
  AKANTU_DEBUG_IN();

  this->internal_force->clear();

  this->synchronize(_gst_pfm_damage);
  auto & fem = this->getFFEngine();

  for (auto ghost_type: ghost_types) {

    for (auto type: mesh.elementTypes(spatial_dimension, ghost_type)) {
      UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

      auto & strain_history_on_qpoints_vect = strain_history_on_qpoints(type, ghost_type);

      UInt nb_quad_points = strain_history_on_qpoints_vect.size();
      Array<Real> nt_strain_history(nb_quad_points, nb_nodes_per_element);
      fem.computeNtb(strain_history_on_qpoints_vect, nt_strain_history, type, ghost_type);

      UInt nb_elements = mesh.getNbElement(type, ghost_type);
      Array<Real> int_nt_srain_history(nb_elements, nb_nodes_per_element);

      fem.integrate(nt_strain_history, int_nt_srain_history, nb_nodes_per_element, type,
		    ghost_type);

      this->getDOFManager().assembleElementalArrayLocalArray(
	  int_nt_srain_history, *this->internal_force, type, ghost_type, -1);
    }
  }

  AKANTU_DEBUG_OUT();
}
  
/* -------------------------------------------------------------------------- */
void PhaseFieldModel::readMaterials() {
  auto sect = this->getParserSection();

  if (not std::get<1>(sect)) {
    const auto & section = std::get<0>(sect);
    this->parseSection(section);
  }

  gc_on_qpoints.set(g_c);
}
  
/* -------------------------------------------------------------------------- */
void PhaseFieldModel::initFullImpl(const ModelOptions & options) {
  Model::initFullImpl(options);

  readMaterials();
}

  
}
