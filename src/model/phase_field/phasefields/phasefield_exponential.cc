/**
 * @file   phasefield_exponential.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Mon Mar 3 2020
 * @date last modification: Mon Mar 3 2020
 *
 * @brief  Specialization of the phasefield class for exponential field
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
#include "phasefield_exponential.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
PhaseFieldExponential::PhaseFieldExponential(PhaseFieldModel & model,
					      const ID & id)
  : PhaseField(model, id) {
  AKANTU_DEBUG_IN();
  this->initialize();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void PhaseFieldExponential::initialize() {
  Matrix<double> d(spatial_dimension, spatial_dimension);
  d.eye(this->g_c * this->l0);
  damage_energy.set(d);
}

/* -------------------------------------------------------------------------- */
void PhaseFieldExponential::computeDrivingForce(const ElementType & el_type,
						GhostType ghost_type) {

  AKANTU_DEBUG_IN();

  for (auto && tuple : zip(this->phi(el_type, ghost_type),
			   this->driving_force(el_type, ghost_type))) {
    computeDrivingForceOnQuad(std::get<0>(tuple), std::get<1>(tuple));
  }
  
  AKANTU_DEBUG_OUT();
}

  
/* -------------------------------------------------------------------------- */
/*void PhaseFieldExponential::computeAllDrivingEnergy(GhostType ghost_type) {

  AKANTU_DEBUG_IN();

  for (auto type : element_filter.elementTypes(_ghost_type = ghost_type)) {
    computeDamageEnergyDensity(type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
  }*/
  
  
/* -------------------------------------------------------------------------- */
/*void PhaseFieldExponential::assembleStiffnessMatrix(const ElementType & type,
						     GhostType ghost_type) {

  auto nb_element = mesh.getNbElement(type);
  auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  auto nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);

  auto nt_b_n = std::make_unique<Array<Real>>(
      nb_element * nb_quadrature_points,
      nb_nodes_per_element * nb_nodes_per_element, "N^t*b*N");

  // damage_energy_on_qpoints = gc/l0 + phi = scalar
  fem.computeNtbN(damage_energy_density(type, ghost_type), *nt_b_n, 2,
		  type, ghost_type, elem_filter);

  /// compute @f$ K_{\grad d} = \int_e \mathbf{N}^t * \mathbf{w} *
  /// \mathbf{N}@f$
  auto K_n = std::make_unique<Array<Real>>(
	nb_element, nb_nodes_per_element * nb_nodes_per_element, "K_n");

  fem.integrate(*nt_b_n, *K_n, nb_nodes_per_element * nb_nodes_per_element,
		type, ghost_type, elem_filter);

  model.getDOFManager().assembleElementalMatricesToMatrix(
	"K", "damage", *K_n, type, ghost_type, _symmetric, elem_filter);
  
  auto bt_d_b = std::make_unique<Array<Real>>(
       nb_element * nb_quadrature_points,
       nb_nodes_per_element * nb_nodes_per_element, "B^t*D*B");

  // damage_energy_on_qpoints = gc*l0 = scalar
  fem.computeBtDB(damage_energy(type, ghost_type), *bt_d_b, 2, type,
		  ghost_type, elem_filter);

  /// compute @f$ K_{\grad d} = \int_e \mathbf{B}^t * \mathbf{W} *
  /// \mathbf{B}@f$
  auto K_b = std::make_unique<Array<Real>>(
	nb_element, nb_nodes_per_element * nb_nodes_per_element, "K_b");

  fem.integrate(*bt_d_b, *K_b, nb_nodes_per_element * nb_nodes_per_element,
		type, ghost_type, elem_filter);

  model.getDOFManager().assembleElementalMatricesToMatrix(
	"K", "damage", *K_b, type, ghost_type, _symmetric, elem_filter);
}*/


/* -------------------------------------------------------------------------- */
/*void PhaseFieldExponential::computeDamageEnergyDensity(const ElementType & type,
						     GhostType ghost_type) {

  AKANTU_DEBUG_IN();

  for (auto && values :
         zip(make_view(damage_energy(type, ghost_type)),
             make_view(phi(type, ghost_type)))) {

    auto & dam_energy = std::get<0>(values);
    auto & phi_history = std::get<1>(values);
    dam_energy = g_c / l_0 + 2.0 * phi_history;
  }
 

  AKANTU_DEBUG_OUT();
}*/

/* -------------------------------------------------------------------------- */
/*void PhaseFieldExponential::computePhiHistory(const ElementType & type,
					       GhostType ghost_type) {

  AKANTU_DEBUG_IN();
  
  Matrix<Real> strain_plus(spatial_dimension, spatial_dimension);
  Matrix<Real> strain_minus(spatial_dimension, spatial_dimension);
  Matrix<Real> strain_dir(spatial_dimension, spatial_dimension);
  Matrix<Real> strain_diag_plus(spatial_dimension, spatial_dimension);
  Matrix<Real> strain_diag_minus(spatial_dimension, spatial_dimension);

  Vector<Real> strain_values(spatial_dimension);

  Real trace_plus, trace_minus, phi_plus;
  
  for (auto && values :
         zip(make_view(strain(type, ghost_type), spatial_dimension,
                       spatial_dimension),
             make_view(phi(type, ghost_type)))) {

    auto & strain_quad = std::get<0>(values);
    auto & phi_history = std::get<1>(values);

    strain_plus.clear();
    strain_minus.clear();
    strain_dir.clear();
    strain_values.clear();
    strain_diag_plus.clear();
    strain_diag_minus.clear();
    
    strain_quad.eig(strain_values, strain_dir);

    for (UInt i = 0; i < spatial_dimension; i++) {
      strain_diag_plus(i, i) = std::max(Real(0.), strain_values(i));
      strain_diag_minus(i, i) = std::min(Real(0.), strain_values(i));
    }

    Matrix<Real> mat_tmp(spatial_dimension, spatial_dimension);
    Matrix<Real> sigma_plus(spatial_dimension, spatial_dimension);
    Matrix<Real> sigma_minus(spatial_dimension, spatial_dimension);

    mat_tmp.mul<false, true>(strain_diag_plus, strain_dir);
    strain_plus.mul<false, false>(strain_dir, mat_tmp);
    mat_tmp.mul<false, true>(strain_diag_minus, strain_dir);
    strain_minus.mul<false, true>(strain_dir, mat_tmp);

    trace_plus = std::max(Real(0.), strain_quad.trace());
    trace_minus = std::min(Real(0.), strain_quad.trace());
    
    for (UInt i = 0; i < spatial_dimension; i++) {
      for (UInt j = 0; j < spatial_dimension; j++) {
	sigma_plus(i, j) =
	  (i == j) * lambda * trace_plus + 2 * mu * strain_plus(i, j);
	sigma_minus(i, j) =
	  (i == j) * lambda * trace_minus + 2 * mu * strain_minus(i, j);
      }
    }

    phi_plus = 1. / 2 * sigma_plus.doubleDot(strain_quad);
    
    if (phi_plus > phi_history) {
      phi_history = phi_plus;
    }
  }
  }*/


  
/* -------------------------------------------------------------------------- */
void PhaseFieldExponential::updateInternalParameters() {
  PhaseField::updateInternalParameters();
}
  
}
