/**
 * @file   material_iterative_stiffness_reduction.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Thu Feb 18 16:03:56 2016
 *
 * @brief  Implementation of material iterative stiffness reduction
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
#include "material_iterative_stiffness_reduction.hh"
#include "solid_mechanics_model_RVE.hh"
#include <math.h>

__BEGIN_AKANTU__
template<UInt spatial_dimension>

/* -------------------------------------------------------------------------- */
MaterialIterativeStiffnessReduction<spatial_dimension>::MaterialIterativeStiffnessReduction(SolidMechanicsModel & model,
											  const ID & id)  :
  Material(model, id),
  MaterialDamageIterative<spatial_dimension>(model, id),
  eps_u("ultimate_strain", *this),
  D("tangent", *this),
  Gf(0.),
  crack_band_width(0.),
  reduction_constant(0.) {
  AKANTU_DEBUG_IN();

  this->registerParam("Gf",                  Gf,                  _pat_parsable | _pat_modifiable, "fracture energy");
  this->registerParam("crack_band_width",                  crack_band_width,                  _pat_parsable | _pat_modifiable, "crack_band_width");
  this->registerParam("reduction_constant",                  reduction_constant, 2.,                  _pat_parsable | _pat_modifiable, "reduction constant");

  this->eps_u.initialize(1);
  this->D.initialize(1);


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension> 
void MaterialIterativeStiffnessReduction<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialDamageIterative<spatial_dimension>::initMaterial();

  for (ghost_type_t::iterator g = ghost_type_t::begin(); g != ghost_type_t::end(); ++g) {
    GhostType ghost_type = *g;
    /// loop over all types in the material
    typedef ElementTypeMapArray<UInt>:: type_iterator iterator;
    iterator it = this->element_filter.firstType(spatial_dimension, ghost_type, _ek_regular);
    iterator last_type = this->element_filter.lastType(spatial_dimension, ghost_type, _ek_regular);
    /// loop over all types in the filter
    for(; it != last_type; ++it) {
      ElementType el_type = *it;
      /// get the stiffness on each quad point
      Array<Real>::const_iterator<Real> Sc_it = this->Sc(el_type, ghost_type).begin();
      /// get the tangent of the tensile softening on each quad point
      Array<Real>::iterator<Real> D_it = this->D(el_type, ghost_type).begin();
      Array<Real>::iterator<Real> D_end = this->D(el_type, ghost_type).end();
      /// get the ultimate strain on each quad 
      Array<Real>::iterator<Real> eps_u_it = this->eps_u(el_type, ghost_type).begin();
      // compute the tangent and the ultimate strain for each quad
      for (; D_it != D_end; ++Sc_it, ++D_it, ++eps_u_it) {
	*eps_u_it = ( (2. * this->Gf) / (*Sc_it * this->crack_band_width) );
	*D_it = *(Sc_it) / ((*eps_u_it) - ( (*Sc_it) / this->E ) );
      }
    }
  }
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialIterativeStiffnessReduction<spatial_dimension>::computeNormalizedEquivalentStress(const Array<Real> & grad_u,
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

  /// loop over all the quadrature points and compute the equivalent stress
  for(;grad_u_it != grad_u_end; ++ grad_u_it) {
    /// compute the stress
    sigma.clear();
    MaterialElastic<spatial_dimension>::computeStressOnQuad(*grad_u_it, sigma, 0.);
    MaterialDamageIterative<spatial_dimension>::computeDamageAndStressOnQuad(sigma, *dam_it);
    /// compute eigenvalues
    sigma.eig(eigenvalues);

    /// find max eigenvalue and normalize by tensile strength
    *equivalent_stress_it = *(std::max_element(eigenvalues.storage(),
                                               eigenvalues.storage() + spatial_dimension)) / (*Sc_it);
    ++Sc_it;
    ++equivalent_stress_it;
    ++dam_it;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
UInt MaterialIterativeStiffnessReduction<spatial_dimension>::updateDamage() {
  UInt nb_damaged_elements = 0;

  if (this->norm_max_equivalent_stress >= 1.) {

    AKANTU_DEBUG_IN();

    /// update the damage only on non-ghosts elements! Doesn't make sense to update on ghost.
    GhostType ghost_type = _not_ghost;;

    Mesh::type_iterator it = this->model->getFEEngine().getMesh().firstType(spatial_dimension, ghost_type);
    Mesh::type_iterator last_type = this->model->getFEEngine().getMesh().lastType(spatial_dimension, ghost_type);

    /// loop over all the elements
    for(; it != last_type; ++it) {
      ElementType el_type = *it;

      /// get iterators on the needed internal fields
      Array<Real>::const_scalar_iterator equivalent_stress_it = this->equivalent_stress(el_type, ghost_type).begin();
      Array<Real>::const_scalar_iterator equivalent_stress_end = this->equivalent_stress(el_type, ghost_type).end();
      Array<Real>::scalar_iterator dam_it = this->damage(el_type, ghost_type).begin();
      Array<UInt>::scalar_iterator reduction_it = this->reduction_step(el_type, ghost_type).begin();
      Array<Real>::const_scalar_iterator eps_u_it = this->eps_u(el_type, ghost_type).begin();
      Array<Real>::scalar_iterator Sc_it = this->Sc(el_type, ghost_type).begin();
      Array<Real>::const_scalar_iterator D_it = this->D(el_type, ghost_type).begin();

      /// loop over all the quads of the given element type
      for (; equivalent_stress_it != equivalent_stress_end; ++equivalent_stress_it, ++dam_it, ++reduction_it,
	     ++eps_u_it, ++Sc_it, ++D_it) {

	/// check if damage occurs
	if (*equivalent_stress_it >= (1 - this->dam_tolerance) * this->norm_max_equivalent_stress) {

	  /// check if this element can still be damaged
	  if (*reduction_it == this->max_reductions)
	    continue;

	  /// increment the counter of stiffness reduction steps
	  *reduction_it += 1;
	  if (*reduction_it == this->max_reductions)
	    *dam_it =  this->max_damage;
	  else {
	    /// update the damage on this quad
	    *dam_it = 1. - (1./std::pow(this->reduction_constant, *reduction_it));
	    /// update the stiffness on this quad
	    *Sc_it = (*eps_u_it) * (1. - (*dam_it) ) * this->E  * (*D_it)/ ((1. - (*dam_it) ) * this->E + (*D_it));
	  }
	  nb_damaged_elements += 1;
	}
      }
    }
  }

  const SolidMechanicsModelRVE * rve_model = dynamic_cast<const SolidMechanicsModelRVE *>(this->model);
  if (rve_model == NULL) {
    StaticCommunicator & comm = akantu::StaticCommunicator::getStaticCommunicator();
    comm.allReduce(&nb_damaged_elements, 1, _so_sum);
  }
  AKANTU_DEBUG_OUT();
  return nb_damaged_elements;
}

/* -------------------------------------------------------------------------- */



INSTANTIATE_MATERIAL(MaterialIterativeStiffnessReduction);


__END_AKANTU__
