/**
 * @file   material_damage_iterative.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 * @date   Sun Mar  8 12:54:30 2015
 *
 * @brief  Specialization of the class material damage to damage only one gauss
 * point at a time and propagate damage in a linear way. Max principal stress
 * criterion is used as a failure criterion.
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "material_orthotropic_damage_iterative.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialOrthotropicDamageIterative<spatial_dimension>::MaterialOrthotropicDamageIterative(SolidMechanicsModel & model,
								    const ID & id)  :
  Material(model, id),
  MaterialOrthotropicDamage<spatial_dimension>(model, id),
  Sc("Sc", *this),
  equivalent_stress("equivalent_stress", *this),
  equiv_stress_dir("equiv_stress_dir", *this),
  norm_max_equivalent_stress(0) {
  AKANTU_DEBUG_IN();

  this->registerParam("Sc",                  Sc,                  _pat_parsable, "critical stress threshold");
  this->registerParam("prescribed_dam",      prescribed_dam, 0.1, _pat_parsable | _pat_modifiable, "increase of damage in every step" );
  this->registerParam("dam_threshold",       dam_threshold,  0.8,  _pat_parsable | _pat_modifiable, "damage threshold at which damage damage will be set to 1" );
  this->registerParam("dam_tolerance",       dam_tolerance,  0.01,  _pat_parsable | _pat_modifiable, "damage tolerance to decide if quadrature point will be damageed" );
  this->registerParam("max_damage",       max_damage,  0.99999,  _pat_parsable | _pat_modifiable, "maximum damage value" );


  this->use_previous_stress          = true;
  this->use_previous_gradu           = true;
  this->Sc.initialize(1);
  this->equivalent_stress.initialize(1);
  this->equiv_stress_dir.initialize(1);


  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialOrthotropicDamageIterative<spatial_dimension>::computeNormalizedEquivalentStress(
                                                                                   ElementType el_type,
                                                                                   GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  /// Vector to store eigenvalues of current stress tensor
  Vector<Real> eigenvalues(spatial_dimension);

  Array<Real>::const_iterator<Real> Sc_it = Sc(el_type).begin();
  Array<Real>::iterator<Real> equivalent_stress_it
    = equivalent_stress(el_type).begin();

  Array<UInt>::scalar_iterator e_stress_dir_it
    = this->equiv_stress_dir(el_type).begin();

  Array<Real>::const_matrix_iterator sigma_it
    = this->stress(el_type, ghost_type).begin(spatial_dimension,
					      spatial_dimension);

  Array<Real>::const_matrix_iterator sigma_end
    = this->stress(el_type, ghost_type).end(spatial_dimension,
					    spatial_dimension);

  Array<Real>::const_matrix_iterator dam_it
    = this->damage(el_type, ghost_type).begin(spatial_dimension,
					      spatial_dimension);

  Array<Real>::const_matrix_iterator dam_dir_it
    = this->damage_dir_vecs(el_type, ghost_type).begin(spatial_dimension,
						       spatial_dimension);

  Matrix<Real> sigma_rotated(this->spatial_dimension, this->spatial_dimension);
  Matrix<Real> rotation_tmp(this->spatial_dimension, this->spatial_dimension);

  for(;sigma_it != sigma_end; ++sigma_it,
	++e_stress_dir_it, ++dam_it, ++dam_dir_it,
	++Sc_it, ++equivalent_stress_it) {

    /// if damage directions are not defined yet
    if (Math::are_float_equal((*dam_it).trace(), 0 )) {
      /// compute eigenvalues and returns them sorted
      (*sigma_it).eig(eigenvalues);
      *equivalent_stress_it = eigenvalues(0) / *(Sc_it);
    } else {
      /// Rotate the stress in the coordinate system of the damage
      this->rotateIntoDamageFrame(*sigma_it,
				  sigma_rotated,
				  *dam_dir_it,
				  rotation_tmp);

      /// find the largest normal stress
      Real max_stress = 0;
      for (UInt i  = 0; i < this->spatial_dimension; ++i) {
	if (sigma_rotated(i,i) > max_stress) { 
	  max_stress = sigma_rotated(i,i); 
	  *e_stress_dir_it = i; 
	}
      }
      //(*sigma_it).eig(eigenvalues);
      //*equivalent_stress_it = eigenvalues(0) / *(Sc_it);
      *equivalent_stress_it = sigma_rotated(*e_stress_dir_it, *e_stress_dir_it) / *(Sc_it);
    }
  }

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialOrthotropicDamageIterative<spatial_dimension>::computeAllStresses(GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  /// reset normalized maximum equivalent stress
  if(ghost_type==_not_ghost)
    norm_max_equivalent_stress = 0;
 
  MaterialOrthotropicDamage<spatial_dimension>::computeAllStresses(ghost_type);

  /// find global Gauss point with highest stress
  StaticCommunicator & comm = akantu::StaticCommunicator::getStaticCommunicator();
  comm.allReduce(&norm_max_equivalent_stress, 1, _so_max);

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialOrthotropicDamageIterative<spatial_dimension>::findMaxNormalizedEquivalentStress(ElementType el_type,
										   GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  if(ghost_type==_not_ghost) {

    // const Array<Real> & e_stress = equivalent_stress(el_type);

    // if (e_stress.begin() != e_stress.end() ) {
    //   Array<Real>::const_iterator<Real> equivalent_stress_it_max = std::max_element(e_stress.begin(),e_stress.end());
    //   /// check if max equivalent stress for this element type is greater than the current norm_max_eq_stress
    //   if (*equivalent_stress_it_max > norm_max_equivalent_stress) 
    // 	norm_max_equivalent_stress = *equivalent_stress_it_max;
    // }
    const Array<Real> & e_stress = equivalent_stress(el_type);
    Array<Real>::const_iterator<Real> equivalent_stress_it = e_stress.begin();
    Array<Real>::const_iterator<Real> equivalent_stress_end = e_stress.end();
    Array<Real> & dam = this->damage(el_type);
    Array<Real>::const_matrix_iterator dam_it = dam.begin(this->spatial_dimension, this->spatial_dimension);


    for (; equivalent_stress_it != equivalent_stress_end; ++equivalent_stress_it, ++dam_it ) {
      /// check if max equivalent stress for this element type is greater than the current norm_max_eq_stress and if the element is not already fully damaged
      if (*equivalent_stress_it > norm_max_equivalent_stress && (*dam_it).trace() < std::min(max_damage, this->spatial_dimension / this->eta * max_damage)) {
	norm_max_equivalent_stress = *equivalent_stress_it;
      }

    }



  }
  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialOrthotropicDamageIterative<spatial_dimension>::computeStress(ElementType el_type,
									  GhostType ghost_type) {
  AKANTU_DEBUG_IN();


  MaterialOrthotropicDamage<spatial_dimension>::computeStress(el_type, ghost_type);

  Array<Real>::matrix_iterator damage_iterator = this->damage(el_type, ghost_type).begin(this->spatial_dimension, this->spatial_dimension);
  Array<Real>::matrix_iterator damage_dir_it = this->damage_dir_vecs(el_type, ghost_type).begin(this->spatial_dimension, this->spatial_dimension);

  /// for the computation of the Cauchy stress the matrices (1-D) and
  /// (1-D)^(1/2) are needed. For the formulation see Engineering
  /// Damage Mechanics by Lemaitre and Desmorat.

  Matrix<Real> one_minus_D(this->spatial_dimension, this->spatial_dimension);
  Matrix<Real> sqrt_one_minus_D(this->spatial_dimension, this->spatial_dimension);
  Matrix<Real> one_minus_D_rotated(this->spatial_dimension, this->spatial_dimension);
  Matrix<Real> sqrt_one_minus_D_rotated(this->spatial_dimension, this->spatial_dimension);
  Matrix<Real> rotation_tmp(this->spatial_dimension, this->spatial_dimension);

  /// create matrix to store the first term of the computation of the
  /// Cauchy stress
  Matrix<Real> first_term(this->spatial_dimension, this->spatial_dimension);
  Matrix<Real> third_term(this->spatial_dimension, this->spatial_dimension);
  
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  /// rotate the tensors from the damage principal coordinate system to the CS of the computation
  if ( !(Math::are_float_equal((*damage_iterator).trace(), 0)) ) {
    /// compute (1-D) and (1-D)^1/2
    this->computeOneMinusD(one_minus_D, *damage_iterator);
    this->computeSqrtOneMinusD(one_minus_D, sqrt_one_minus_D);

    this->rotateIntoComputationFrame(one_minus_D,
				     one_minus_D_rotated,
				     *damage_dir_it,
				     rotation_tmp);

    this->rotateIntoComputationFrame(sqrt_one_minus_D,
				     sqrt_one_minus_D_rotated,
				     *damage_dir_it,
				     rotation_tmp);
  } else {
    this->computeOneMinusD(one_minus_D_rotated, *damage_iterator);
    this->computeSqrtOneMinusD(one_minus_D_rotated, sqrt_one_minus_D_rotated);
  }

  computeDamageAndStressOnQuad(sigma,
			       one_minus_D_rotated,
			       sqrt_one_minus_D_rotated,
			       *damage_iterator,
			       first_term,
			       third_term);

  ++damage_dir_it;
  ++damage_iterator;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  computeNormalizedEquivalentStress(el_type, ghost_type);
  norm_max_equivalent_stress = 0;
  findMaxNormalizedEquivalentStress(el_type, ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
UInt MaterialOrthotropicDamageIterative<spatial_dimension>::updateDamage() {
  UInt nb_damaged_elements = 0;
  AKANTU_DEBUG_ASSERT(prescribed_dam > 0.,
		      "Your prescribed damage must be greater than zero");

  if (norm_max_equivalent_stress >= 1.) {

    AKANTU_DEBUG_IN();
    GhostType ghost_type = _not_ghost;;

    Mesh::type_iterator it = this->model->getFEEngine().getMesh().firstType(spatial_dimension, ghost_type);
    Mesh::type_iterator last_type = this->model->getFEEngine().getMesh().lastType(spatial_dimension, ghost_type);

    for(; it != last_type; ++it) {
      ElementType el_type = *it;

      const Array<Real> & e_stress = equivalent_stress(el_type);
      Array<Real>::const_iterator<Real> equivalent_stress_it = e_stress.begin();
      Array<Real>::const_iterator<Real> equivalent_stress_end = e_stress.end();
      Array<Real> & dam = this->damage(el_type);
      Array<Real>::matrix_iterator dam_it = dam.begin(this->spatial_dimension, this->spatial_dimension);
      Array<Real>::const_matrix_iterator stress_it =				
	this->stress(el_type, ghost_type).begin(this->spatial_dimension,     
					       this->spatial_dimension);
      Array<Real>::matrix_iterator damage_directions_it = 
	this->damage_dir_vecs(el_type, ghost_type).begin(this->spatial_dimension,
							 this->spatial_dimension);

      Array<UInt>::scalar_iterator e_stress_dir_it = 
	this->equiv_stress_dir(el_type, ghost_type).begin();

      for (; equivalent_stress_it != equivalent_stress_end; ++equivalent_stress_it, ++dam_it, ++stress_it, ++damage_directions_it, ++e_stress_dir_it) {

	/// check if damage occurs
	if (*equivalent_stress_it >= (1-dam_tolerance)*norm_max_equivalent_stress) {

	  /// check if damage coordinate system has been defined already

	  if (Math::are_float_equal((*dam_it).trace(), 0)) {
	    /// first time damage occurs at this Gauss point -> define
	    /// directions of damage -> directions of current strain
	    Vector<Real> eigenvalues(this->spatial_dimension);
	    /// eigenvalues are returned sorted, with the first entry
	    /// corresponding to the largest eigenvalue
      
	    (*stress_it).eig(eigenvalues, *damage_directions_it);

	    /// damage occurs in the direction of the maximum principale strain
	    (*dam_it)(0, 0) += prescribed_dam;
	  }
	  
	  else {
	    /// damage coordinate system has already been
	    /// defined. 
	    /// Direction of damage increase:
	    UInt dir = *e_stress_dir_it;
	    
	    /// margin for damage increase: ensure that trace(D) <= dim/eta * Dmax as discussed with Fabrice

	    Real margin = std::min(max_damage, this->spatial_dimension / this->eta * max_damage) - (*dam_it).trace();
	    if ((*dam_it)(dir,dir) < dam_threshold) {
	      if (margin < prescribed_dam)
		(*dam_it)(dir, dir) += margin; 
	      else 
		(*dam_it)(dir, dir) += prescribed_dam;
	    }
	    else 
	      (*dam_it)(dir, dir) += margin;
	  }
	  nb_damaged_elements += 1;
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
void MaterialOrthotropicDamageIterative<spatial_dimension>::updateEnergiesAfterDamage(ElementType el_type, GhostType ghost_type) {
  MaterialOrthotropicDamage<spatial_dimension>::updateEnergies(el_type, ghost_type);
}

/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialOrthotropicDamageIterative);


__END_AKANTU__
