/**
 * @file   material_igfem_iterative_stiffness_reduction.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Thu Mar 10 08:37:43 2016
 *
 * @brief  Implementation of igfem material iterative stiffness reduction
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_igfem_iterative_stiffness_reduction.hh"
#include <math.h>

__BEGIN_AKANTU__
template<UInt spatial_dimension>

/* -------------------------------------------------------------------------- */
MaterialIGFEMIterativeStiffnessReduction<spatial_dimension>::MaterialIGFEMIterativeStiffnessReduction(SolidMechanicsModel & model,
											  const ID & id)  :
  Material(model, id),
  MaterialIGFEMSawToothDamage<spatial_dimension>(model, id),
  eps_u("ultimate_strain", *this),
  reduction_step("damage_step", *this),
  D("tangent", *this),
  Gf(0.),
  el_h(0.),
  max_reductions(0),
  reduction_constant(0.) {
  AKANTU_DEBUG_IN();

  this->registerParam("Gf",                  Gf,                  _pat_parsable, "fracture energy");
  this->registerParam("element_size",                  el_h,                  _pat_parsable, "element size");
  this->registerParam("max_reductions",                  max_reductions, UInt(10),                  _pat_parsable | _pat_modifiable, "max reductions");
  this->registerParam("reduction_constant",                  reduction_constant, 2.,                  _pat_parsable | _pat_modifiable, "reduction constant");

  this->eps_u.initialize(1);
  this->D.initialize(1);
  this->reduction_step.initialize(1);

  this->internals_to_transfer.push_back("ultimate_strain");
  this->internals_to_transfer.push_back("tangent");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialIGFEMIterativeStiffnessReduction<spatial_dimension>::computeNormalizedEquivalentStress(const Array<Real> & grad_u,
											      ElementType el_type,
											      GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  /// storage for the current stress
  Matrix<Real> sigma(spatial_dimension, spatial_dimension);
  /// Vector to store eigenvalues of current stress tensor
  Vector<Real> eigenvalues(spatial_dimension);

  /// iterators on the needed internal fields
  Array<Real>::const_scalar_iterator Sc_it = this->Sc(el_type, ghost_type).begin();
  Array<Real>::scalar_iterator dam_it = this->damage(el_type, ghost_type).begin();
  Array<Real>::scalar_iterator equivalent_stress_it = this->equivalent_stress(el_type, ghost_type).begin();
  Array<Real>::const_matrix_iterator grad_u_it = grad_u.begin(spatial_dimension,
                                                              spatial_dimension);
  Array<Real>::const_matrix_iterator grad_u_end = grad_u.end(spatial_dimension,
                                                             spatial_dimension);
  Real * mu_ptr = this->mu(el_type, ghost_type).storage();
  Real * lambda_ptr = this->lambda(el_type, ghost_type).storage();
  
  /// loop over all the quadrature points and compute the equivalent stress
  for(;grad_u_it != grad_u_end; ++grad_u_it) {
    /// compute the stress
    sigma.clear();
    MaterialIGFEMElastic<spatial_dimension>::computeStressOnQuad(*grad_u_it, sigma, *lambda_ptr, *mu_ptr);
    MaterialIGFEMSawToothDamage<spatial_dimension>::computeDamageAndStressOnQuad(sigma, *dam_it);
    /// compute eigenvalues
    sigma.eig(eigenvalues);

    /// find max eigenvalue and normalize by tensile strength
    *equivalent_stress_it = *(std::max_element(eigenvalues.storage(),
                                               eigenvalues.storage() + spatial_dimension)) / (*Sc_it);
    ++Sc_it;
    ++equivalent_stress_it;
    ++dam_it;
    ++lambda_ptr;
    ++mu_ptr;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
UInt MaterialIGFEMIterativeStiffnessReduction<spatial_dimension>::updateDamage() {
  UInt nb_damaged_elements = 0;

  if (this->norm_max_equivalent_stress >= 1.) {

    AKANTU_DEBUG_IN();

    /// update the damage only on non-ghosts elements! Doesn't make sense to update on ghost.
    GhostType ghost_type = _not_ghost;;

    Mesh::type_iterator it = this->model->getFEEngine().getMesh().firstType(spatial_dimension, ghost_type);
    Mesh::type_iterator last_type = this->model->getFEEngine().getMesh().lastType(spatial_dimension, ghost_type);

    /// get the Young's modulus of the damageable sub-material 
    ID mat_name = this->sub_material_names[1];
    Real E = dynamic_cast<MaterialElastic<spatial_dimension> & >(this->model->getMaterial(mat_name)).getYoungsModulus();
    
    /// loop over all the elements
    for(; it != last_type; ++it) {
      ElementType el_type = *it;

      /// get iterators on the needed internal fields
      Array<UInt> & sub_mat = this->sub_material(el_type, ghost_type);
      Array<UInt>::iterator<UInt> sub_mat_it = sub_mat.begin();
      Array<Real>::const_scalar_iterator equivalent_stress_it = this->equivalent_stress(el_type, ghost_type).begin();
      Array<Real>::const_scalar_iterator equivalent_stress_end = this->equivalent_stress(el_type, ghost_type).end();
      Array<Real>::scalar_iterator dam_it = this->damage(el_type, ghost_type).begin();
      Array<UInt>::scalar_iterator reduction_it = this->reduction_step(el_type, ghost_type).begin();
      Array<Real>::const_scalar_iterator eps_u_it = this->eps_u(el_type, ghost_type).begin();
      Array<Real>::scalar_iterator Sc_it = this->Sc(el_type, ghost_type).begin();
      Array<Real>::const_scalar_iterator D_it = this->D(el_type, ghost_type).begin();

      /// loop over all the elements of the given type
      UInt nb_element = this->element_filter(el_type, ghost_type).getSize();
      UInt nb_quads = this->fem->getNbIntegrationPoints(el_type, ghost_type);
      bool damage_element = false;
      for (UInt e = 0; e < nb_element; ++e) {
	damage_element = false;
	/// check if damage occurs in the element
	for (UInt q = 0; q < nb_quads; ++q, ++reduction_it, ++sub_mat_it, ++equivalent_stress_it) {
	  if (*equivalent_stress_it >= (1 - this->dam_tolerance) * this->norm_max_equivalent_stress && *sub_mat_it != 0) {
	    /// check if this element can still be damaged
	    if (*reduction_it == this->max_reductions)
	      continue;	    
	    damage_element = true;
	  }
	}
	
	if (damage_element) {
	  /// damage the element
	  nb_damaged_elements += 1;
	  sub_mat_it -= nb_quads;
	  reduction_it -= nb_quads;
	  for (UInt q = 0; q < nb_quads; ++q) {
	    if (*sub_mat_it) {
	      /// increment the counter of stiffness reduction steps
	      *reduction_it += 1;
	      if (*reduction_it == this->max_reductions)
		*dam_it =  this->max_damage;
	      else {
		/// update the damage on this quad
		*dam_it = 1. - (1./std::pow(this->reduction_constant, *reduction_it));
		/// update the stiffness on this quad
		*Sc_it = (*eps_u_it) * (1. - (*dam_it) ) * E  * (*D_it)/ ((1. - (*dam_it) ) * E + (*D_it));
	      }
	    }	      
	    ++sub_mat_it;
	    ++dam_it;
	    ++reduction_it;
	    ++eps_u_it;
	    ++Sc_it;
	    ++D_it;
	  }	    
	}
	else {
	  dam_it += nb_quads;
	  eps_u_it += nb_quads;
	  Sc_it += nb_quads;
	  D_it +=nb_quads;
	}
      }
    }
  }
  StaticCommunicator & comm = akantu::StaticCommunicator::getStaticCommunicator();
  comm.allReduce(&nb_damaged_elements, 1, _so_sum);
  AKANTU_DEBUG_OUT();
  return nb_damaged_elements;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialIGFEMIterativeStiffnessReduction<spatial_dimension>::onElementsAdded(__attribute__((unused)) const Array<Element> & element_list,
										  __attribute__((unused)) const NewElementsEvent & event) {

  MaterialIGFEMSawToothDamage<spatial_dimension>::onElementsAdded(element_list, event);
  /// set the correct damage iteration step (is UInt->cannot be interpolated)
  Real val = 0.;
  for (ghost_type_t::iterator g = ghost_type_t::begin(); g != ghost_type_t::end(); ++g) {
    GhostType ghost_type = *g;
    /// loop over all types in the material
    typedef ElementTypeMapArray<UInt>:: type_iterator iterator;
    iterator it = this->element_filter.firstType(spatial_dimension, ghost_type, _ek_igfem);
    iterator last_type = this->element_filter.lastType(spatial_dimension, ghost_type, _ek_igfem);
    /// loop over all types in the filter
    for(; it != last_type; ++it) {
      const ElementType el_type = *it;
      Array<Real>::scalar_iterator dam_it = this->damage(el_type, ghost_type).begin();
      Array<UInt>::scalar_iterator reduction_it = this->reduction_step(el_type, ghost_type).begin();
      UInt nb_element = this->element_filter(el_type, ghost_type).getSize();
      UInt nb_quads = this->fem->getNbIntegrationPoints(el_type);
      UInt * sub_mat_ptr = this->sub_material(el_type, ghost_type).storage();
      for (UInt q = 0; q < nb_element * nb_quads; ++q, ++sub_mat_ptr) {
	if (*sub_mat_ptr)
	  for (UInt i = 0; i < this->max_reductions; ++i) {
	    val = 1 - (1./std::pow(this->reduction_constant, i));
	    if (Math::are_float_equal(val, *dam_it))
	      *reduction_it = i;
	  }
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
  
INSTANTIATE_MATERIAL(MaterialIGFEMIterativeStiffnessReduction);


__END_AKANTU__
