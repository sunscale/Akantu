/**
 * @file   material_phase_field.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Tue Aug 21 2018
 * @date last modification: Tue Aug 21 2018
 *
 * @brief  Implementation of the common part of the phasefield material class
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
#include "material_phase_field.hh"
#include "phase_field_model.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
MaterialPhaseField::MaterialPhaseField(PhaseFieldModel & model, const ID & id)
    : Memory(id, model.getMemoryID()), Parsable(ParserType::_material, id),
      is_init(false), fem(model.getFEEngine()),
      name(""), model(model),
      spatial_dimension(this->model.getSpatialDimension()),
      element_filter("element_filter", id, this->memory_id),
      damage("damage", *this), 
      stress("stress", *this), eigengradu("eigen_grad_u", *this),
      gradu("grad_u", *this), green_strain("green_strain", *this),
      piola_kirchhoff_2("piola_kirchhoff_2", *this),
      dissipated_energy("dissipated_energy", *this), is_non_local(false),
      use_previous_stress(false), use_previous_gradu(false),
      interpolation_inverse_coordinates("interpolation inverse coordinates",
                                        *this),
      interpolation_points_matrices("interpolation points matrices", *this) {
  AKANTU_DEBUG_IN();

  element_filter.initialize(model.getMesh(),
			    _spatial_dimension = spatial_dimension);

  this->initialize();
  AKANTU_DEBUG_OUT();
}

MaterialPhaseField::MaterialPhaseField(PhaseFieldModel & model, UInt dim, const Mesh & mesh,
                   FEEngine & fe_engine, const ID & id)
    : Memory(id, model.getMemoryID()), Parsable(ParserType::_material, id),
      is_init(false), fem(fe_engine), finite_deformation(false), name(""),
      model(model), spatial_dimension(dim),
      element_filter("element_filter", id, this->memory_id),
      damage("damage", *this), 
      stress("stress", *this, dim, fe_engine, this->element_filter),
      eigengradu("eigen_grad_u", *this, dim, fe_engine, this->element_filter),
      gradu("gradu", *this, dim, fe_engine, this->element_filter),
      green_strain("green_strain", *this, dim, fe_engine, this->element_filter),
      piola_kirchhoff_2("piola_kirchhoff_2", *this, dim, fe_engine,
                        this->element_filter),
      potential_energy("potential_energy", *this, dim, fe_engine,
                       this->element_filter),
      is_non_local(false), use_previous_stress(false),
      use_previous_gradu(false),
      interpolation_inverse_coordinates("interpolation inverse_coordinates",
                                        *this, dim, fe_engine,
                                        this->element_filter),
      interpolation_points_matrices("interpolation points matrices", *this, dim,
                                    fe_engine, this->element_filter) {

  AKANTU_DEBUG_IN();
  element_filter.initialize(mesh, _spatial_dimension = spatial_dimension);
  
  this->initialize();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
MaterialPhaseField::~MaterialPhaseField() = default;

/* -------------------------------------------------------------------------- */
void MaterialPhaseField::initialize() {
  registerParam("E", E, Real(0.), _pat_parsable | _pat_modifiable,
                "Young's modulus");
  registerParam("nu", nu, Real(0.5), _pat_parsable | _pat_modifiable,
		"Poisson's ratio");
  registerParam("Gc", gc, Real(0.), _pat_parsable | _pat_modifiable,
		"Griffith's fracture energy");
  registerParam("l0", l0, Real(0.01), _pat_parsable | _pat_modifiable,
		"length scale");

  /// allocate gradu stress for local elements
  eigengradu.initialize(spatial_dimension * spatial_dimension);
  gradu.initialize(spatial_dimension * spatial_dimension);
  stress.initialize(spatial_dimension * spatial_dimension);

  dissipated_energy.initialize(1);
  damage.initialize(0);
  
  this->model.registerEventHandler(*this);
}

/* -------------------------------------------------------------------------- */
void MaterialPhaseField::initMaterial() {
  AKANTU_DEBUG_IN();

  if (finite_deformation) {
    this->piola_kirchhoff_2.initialize(spatial_dimension * spatial_dimension);
    if (use_previous_stress)
      this->piola_kirchhoff_2.initializeHistory();
    this->green_strain.initialize(spatial_dimension * spatial_dimension);
  }

  if (use_previous_stress)
    this->stress.initializeHistory();
  if (use_previous_gradu)
    this->gradu.initializeHistory();

  for (auto it = internal_vectors_real.begin();
       it != internal_vectors_real.end(); ++it)
    it->second->resize();

  for (auto it = internal_vectors_uint.begin();
       it != internal_vectors_uint.end(); ++it)
    it->second->resize();

  for (auto it = internal_vectors_bool.begin();
       it != internal_vectors_bool.end(); ++it)
    it->second->resize();

  is_init = true;

  updateInternalParameters();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialPhaseField::savePreviousState() {
  AKANTU_DEBUG_IN();

  for (auto pair : internal_vectors_real)
    if (pair.second->hasHistory())
      pair.second->saveCurrentValues();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialPhaseField::restorePreviousState() {
  AKANTU_DEBUG_IN();

  for (auto pair : internal_vectors_real)
    if (pair.second->hasHistory())
      pair.second->restorePreviousValues();

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
/**
 * Compute  the damage  matrix by  assembling @f$\int_{\omega}  N^t  \times D
 * \times N d\omega @f$
 *
 * @param[in] current_position nodes postition + displacements
 * @param[in] ghost_type compute the residual for _ghost or _not_ghost element
 */
void MaterialPhaseField::assembleDamageMatrix(GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  
  UInt spatial_dimension = model.getSpatialDimension();

  for (auto type : element_filter.elementTypes(spatial_dimension, ghost_type)) {
    switch (spatial_dimension) {
    case 1: {
      assembleDamageMatrix<1>(type, ghost_type);
      break;
    }
    case 2: {
      assembleDamageMatrix<2>(type, ghost_type);
      break;
    }
    case 3: {
      assembleDamageMatrix<3>(type, ghost_type);
      break;
    }  
    }
  }
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void MaterialPhaseField::assembleDamageMatrix(const ElementType & type,
					      GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Array<UInt> & elem_filter = element_filter(type, ghost_type);
  if (elem_filter.size == 0) {
    AKANTU_DEBUG_OUT();
    return;
  }

  auto nb_element = elem_filter.size();
  auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  auto nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);

  auto nt_d_n = std::make_unique<Array<Real>>(
	nb_element * nb_quadrature_points,
	nb_nodes_per_element * nb_nodes_per_element, "N^t*D*N");

  fem.computeNtDN(conductivity_on_qpoints(tyep, ghost_type), *nt_d_n, type,
		  ghost_type);

  auto K_d = std::make_unique<Array<Real> >(
      nb_element, nb_nodes_per_element * nb_nodes_per_element, "K_d" );

  fem.integrate(*nt_d_n, *K_d, nb_nodes_per_element * nb_nodes_per_element, "K_d");

  model.getDOFManager().assembleElementalMatricesToMatrix(
	       "K", "damage", *K_d, type, ghost_type, _symmetric, elem_filter);
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void MaterialPhaseField::assembleDamageGradMatrix(const ElementType & type,
						  GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Array<UInt> & elem_filter = element_filter(type, ghost_type);
  if (elem_filter.size == 0) {
    AKANTU_DEBUG_OUT();
    return;
  }

  auto nb_element = elem_filter.size();
  auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  auto nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);

  auto bt_d_b = std::make_unique<Array<Real>>(nb_element * nb_quadrature_points,
					      nb_nodes_per_element * nb_nodes_per_element, "B^t*D*B");

  /// compute @f$ K_{\grad d} = \int_e \mathbf{B}^t * \mathbf{W} * \mathbf{B}@f$
  fem.computeBtDB(conductivity_on_qpoints(type, ghost_type), *bt_d_b, type,
		  ghost_type);

  auto K_b = std::make_unique<Array<Real> >(
	nb_element, nb_nodes_per_element * nb_nodes_per_element, "K_b");

  fem.integrate(*bt_d_b, *K_b, nb_nodes_per_element * nb_nodes_per_element,
		type, ghost_type);
  model.getDOFManager().assembleElementalMatricsToMatrix(
	      "K", "damage", *K_b, type, ghost_type, _symmetric, elem_filter);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialPhaseField::computeHistoryField(const ElementType & el_type,
					     Array<Real> & history_array,
					     GhostType ghost_type) {
  model.getDisplacement();
}
  
/* -------------------------------------------------------------------------- */
void MaterialPhaseField::computeHistoryFieldOnQuadPoints(
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

