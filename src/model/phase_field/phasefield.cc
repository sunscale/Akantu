/**
 * @file   pahsefield.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Mon Mar 2 2020
 * @date last modification: Mon Mar 2 2020
 *
 * @brief  Implementation of the common part of the phasefield class
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
#include "phasefield.hh"
#include "phase_field_model.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
PhaseField::PhaseField(PhaseFieldModel & model, const ID & id)
  : Memory(id, model.getMemoryID()), Parsable(ParserType::_phasefield, id),
    fem(model.getFEEngine()), name(""), model(model),
    spatial_dimension(this->model.getSpatialDimension()),
    element_filter("element_filter", id, this->memory_id),
    damage("damage", *this), phi("phi", *this),
    strain("strain", *this), driving_force("driving_force", *this),
    damage_energy("damage_energy", *this),
    damage_energy_density("damage_energy_density", *this) {

  AKANTU_DEBUG_IN();

  /// for each connectivity types allocate the element filer array of the
  /// material
  element_filter.initialize(model.getMesh(),
                            _spatial_dimension = spatial_dimension,
                            _element_kind = _ek_regular);
  this->initialize();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
  PhaseField::PhaseField(PhaseFieldModel & model, UInt dim, const Mesh & mesh,
			 FEEngine & fe_engine, const ID & id)
  : Memory(id, model.getMemoryID()), Parsable(ParserType::_phasefield, id),
    fem(fe_engine), name(""), model(model),
    spatial_dimension(this->model.getSpatialDimension()),
    element_filter("element_filter", id, this->memory_id),
    damage("damage", *this, dim, fe_engine, this->element_filter),
    phi("phi", *this, dim, fe_engine, this->element_filter),
    strain("strain", *this, dim, fe_engine, this->element_filter),
    driving_force("driving_force", *this, dim, fe_engine, this->element_filter),
    damage_energy("damage_energy", *this, dim, fe_engine, this->element_filter),
    damage_energy_density("damage_energy_density", *this, dim, fe_engine,
			  this->element_filter){

  AKANTU_DEBUG_IN();

  /// for each connectivity types allocate the element filer array of the
  /// material
  element_filter.initialize(mesh,
                            _spatial_dimension = spatial_dimension,
                            _element_kind = _ek_regular);
  this->initialize();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */  
PhaseField::~PhaseField() = default;

/* -------------------------------------------------------------------------- */
void PhaseField::initialize() {
  registerParam("name", name, std::string(), _pat_parsable | _pat_readable);
  registerParam("l0", l0, Real(0.), _pat_parsmod,
		"length scale parameter");
  registerParam("gc", g_c, _pat_parsmod,
		"critical local fracture energy density");
  registerParam("E", E, _pat_parsmod, "Young's modulus");
  registerParam("nu", nu, _pat_parsmod, "Poisson ratio");

  damage.initialize(0);
  
  phi.initialize(1);
  driving_force.initialize(1);

  strain.initialize(spatial_dimension * spatial_dimension);
  damage_energy_density.initialize(1);
  damage_energy.initialize(1);
}

/* -------------------------------------------------------------------------- */
void PhaseField::initPhaseField() {
  AKANTU_DEBUG_IN();

  this->phi.initializeHistory();

  this->resizeInternals();

  updateInternalParameters();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void PhaseField::resizeInternals() {
  AKANTU_DEBUG_IN();
  for (auto it = internal_vectors_real.begin();
       it != internal_vectors_real.end(); ++it)
    it->second->resize();

  for (auto it = internal_vectors_uint.begin();
       it != internal_vectors_uint.end(); ++it)
    it->second->resize();

  for (auto it = internal_vectors_bool.begin();
       it != internal_vectors_bool.end(); ++it)
    it->second->resize();
  AKANTU_DEBUG_OUT();
}
  

/* -------------------------------------------------------------------------- */
void PhaseField::updateInternalParameters() {
  this->lambda = this->nu * this->E / ((1 + this->nu) * (1 - 2 * this->nu));
  this->mu = this->E / (2 * (1 + this->nu));
}

/* -------------------------------------------------------------------------- */
void PhaseField::computeAllDrivingForces(GhostType ghost_type) {

  AKANTU_DEBUG_IN();

  UInt spatial_dimension = model.getSpatialDimension();

  for (const auto & type :
       element_filter.elementTypes(spatial_dimension, ghost_type)) {
    Array<UInt> & elem_filter = element_filter(type, ghost_type);

    if (elem_filter.size() == 0)
      continue;
    
    computeDrivingForce(type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
  
}
  
/* -------------------------------------------------------------------------- */
void PhaseField::assembleInternalForces(GhostType ghost_type) {

  AKANTU_DEBUG_IN();

  auto & internal_forces = model.getInternalForce();

  for (auto type : element_filter.elementTypes(_ghost_type = ghost_type)) {
  
    Array<UInt> & elem_filter = element_filter(type, ghost_type);
    if (elem_filter.size() == 0)
      continue;
 
    UInt nb_element = elem_filter.size();
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    UInt nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);

    Array<Real> nt_driving_force(nb_quadrature_points, nb_nodes_per_element);
    fem.computeNtb(driving_force(type, ghost_type), nt_driving_force, type,
		   ghost_type, elem_filter);
    
    Array<Real> int_nt_driving_force(nb_element, nb_nodes_per_element);
    fem.integrate(nt_driving_force, int_nt_driving_force,
		  nb_nodes_per_element, type, ghost_type, elem_filter);

    model.getDOFManager().assembleElementalArrayLocalArray(
	 int_nt_driving_force, internal_forces, type, ghost_type, 1,
	 elem_filter);
  }

  AKANTU_DEBUG_OUT();
}
 
/* -------------------------------------------------------------------------- */
void PhaseField::assembleStiffnessMatrix(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Assemble the new stiffness matrix");

  for (auto type : element_filter.elementTypes(spatial_dimension, ghost_type)) {

    Array<UInt> & elem_filter = element_filter(type, ghost_type);
    if (elem_filter.size() == 0) {
      AKANTU_DEBUG_OUT();
      return;
    }
    
    auto nb_element = elem_filter.size();
    auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    auto nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);


    auto nt_b_n = std::make_unique<Array<Real>>(
        nb_element * nb_quadrature_points,
        nb_nodes_per_element * nb_nodes_per_element, "N^t*b*N");

    auto bt_d_b = std::make_unique<Array<Real>>(
        nb_element * nb_quadrature_points,
        nb_nodes_per_element * nb_nodes_per_element, "B^t*D*B");

    // damage_energy_on_qpoints = gc/l0 + phi = scalar
    auto & damage_energy_density_vect =
      damage_energy_density(type, ghost_type);

    // damage_energy_on_qpoints = gc*l0 = scalar
    auto & damage_energy_vect =
      damage_energy(type, ghost_type);

    fem.computeBtDB(damage_energy_vect, *bt_d_b, 2, type, ghost_type,
		    elem_filter);
    
    fem.computeNtbN(damage_energy_density_vect, *nt_b_n, 2, type, ghost_type,
		    elem_filter);

    /// compute @f$ K_{\grad d} = \int_e \mathbf{N}^t * \mathbf{w} *
    /// \mathbf{N}@f$
    auto K_n = std::make_unique<Array<Real>>(
        nb_element, nb_nodes_per_element * nb_nodes_per_element, "K_n");

    fem.integrate(*nt_b_n, *K_n, nb_nodes_per_element * nb_nodes_per_element,
                  type, ghost_type, elem_filter);

    model.getDOFManager().assembleElementalMatricesToMatrix(
		  "K", "damage", *K_n, type, _not_ghost, _symmetric, elem_filter);

    /// compute @f$ K_{\grad d} = \int_e \mathbf{B}^t * \mathbf{W} *
    /// \mathbf{B}@f$
    auto K_b = std::make_unique<Array<Real>>(
        nb_element, nb_nodes_per_element * nb_nodes_per_element, "K_b");

    fem.integrate(*bt_d_b, *K_b, nb_nodes_per_element * nb_nodes_per_element,
                  type, ghost_type, elem_filter);

    model.getDOFManager().assembleElementalMatricesToMatrix(
		 "K", "damage", *K_b, type, _not_ghost, _symmetric, elem_filter);
  }

  AKANTU_DEBUG_OUT();
}  

/* -------------------------------------------------------------------------- */
void PhaseField::beforeSolveStep() {
  this->savePreviousState();
}

/* -------------------------------------------------------------------------- */
void PhaseField::afterSolveStep() {

}
  
/* -------------------------------------------------------------------------- */
void PhaseField::savePreviousState() {
  AKANTU_DEBUG_IN();

  for (auto pair : internal_vectors_real)
    if (pair.second->hasHistory())
      pair.second->saveCurrentValues();

  AKANTU_DEBUG_OUT();
}
  
  
/* -------------------------------------------------------------------------- */
void PhaseField::printself(std::ostream & stream, int indent) const {
 std::string space(indent, AKANTU_INDENT);
  std::string type = getID().substr(getID().find_last_of(':') + 1);

  stream << space << "PhaseField Material " << type << " [" << std::endl;
  Parsable::printself(stream, indent);
  stream << space << "]" << std::endl;
}

  
}
