/**
 * @file   material_igfem_elastic.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  Specializaton of material class for the igfem elastic material
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */
/* -------------------------------------------------------------------------- */
#include "material_igfem_elastic.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt dim>
MaterialIGFEMElastic<dim>::MaterialIGFEMElastic(SolidMechanicsModel & model, const ID & id)  :
  MaterialIGFEM(model, id),
  E(2),
  nu(2),
  lambda("lambda", *this),
  mu("mu", *this),
  kpa("kappa", *this),
  lambda_default(0),
  mu_default(0),
  kpa_default(0) {
  AKANTU_DEBUG_IN();
  this->initialize();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
void MaterialIGFEMElastic<dim>::initialize() {
  this->registerParam("E"      , E      , _pat_parsable | _pat_modifiable, "Young's modulus"        );
  this->registerParam("nu"     , nu     , _pat_parsable | _pat_modifiable, "Poisson's ratio"        );
  this->registerParam("Plane_Stress", plane_stress, false, _pat_parsmod, "Is plane stress");

  this->lambda           .initialize(1);
  this->mu           .initialize(1);
  this->kpa           .initialize(1);
}


/* -------------------------------------------------------------------------- */
template<UInt dim>
void MaterialIGFEMElastic<dim>::initMaterial() {
  AKANTU_DEBUG_IN();

  Material::initMaterial();

  if (dim == 1) this->nu.clear(); /// set the Poisson ratios to zero
  this->updateDefaultInternals();

  /// set the Lamé constants at all quad points to the constants of the first sub-material
  this->lambda.setDefaultValue(lambda_default);
  this->mu.setDefaultValue(mu_default);
  this->kpa.setDefaultValue(kpa_default);


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
void MaterialIGFEMElastic<dim>::updateDefaultInternals(const UInt default_idx) {

  this->lambda_default = this->nu(default_idx) * this->E(default_idx) / ((1 + this->nu(default_idx)) * (1 - 2*this->nu(default_idx)));
  this->mu_default     = this->E(default_idx) / (2 * (1 + this->nu(default_idx)));
  
  this->kpa_default    = this->lambda_default + 2./3. * this->mu_default;
}

/* -------------------------------------------------------------------------- */
template<>
void MaterialIGFEMElastic<2>::updateDefaultInternals(const UInt default_idx) {

  this->lambda_default = this->nu(default_idx) * this->E(default_idx) / ((1 + this->nu(default_idx)) * (1 - 2*this->nu(default_idx)));
  this->mu_default     = this->E(default_idx) / (2 * (1 + this->nu(default_idx)));

  if(this->plane_stress) this->lambda_default = this->nu(default_idx) * this->E(default_idx) / ((1 + this->nu(default_idx))*(1 - this->nu(default_idx)));

  this->kpa_default    = this->lambda_default + 2./3. * this->mu_default;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialIGFEMElastic<spatial_dimension>::computeStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Parent::computeStress(el_type, ghost_type);

  if (!this->finite_deformation) {
    /// get pointer to internals
    Real * lambda_ptr = this->lambda(el_type, ghost_type).storage();
    Real * mu_ptr = this->mu(el_type, ghost_type).storage();
    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);
    this->computeStressOnQuad(grad_u, sigma, *lambda_ptr, *mu_ptr);
    ++lambda_ptr;
    ++mu_ptr;
    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;   
  } else {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialIGFEMElastic<spatial_dimension>::computeTangentModuli(__attribute__((unused)) const ElementType & el_type,
                                                              Array<Real> & tangent_matrix,
                                                              __attribute__((unused)) GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  /// get pointer to internals
  Real * lambda_ptr = this->lambda(el_type, ghost_type).storage();
  Real * mu_ptr = this->mu(el_type, ghost_type).storage();
  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_matrix);
  this->computeTangentModuliOnQuad(tangent, *lambda_ptr, *mu_ptr);
  ++lambda_ptr;
  ++mu_ptr;
  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialIGFEMElastic<spatial_dimension>::computePotentialEnergy(ElementType el_type,
                                                                GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  // MaterialThermal<spatial_dimension>::computePotentialEnergy(el_type, ghost_type);

  // if(ghost_type != _not_ghost) return;
  // Array<Real>::scalar_iterator epot = this->potential_energy(el_type, ghost_type).begin();

  // if (!this->finite_deformation) {
  //   MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  //   this->computePotentialEnergyOnQuad(grad_u, sigma, *epot);
  //   ++epot;

  //   MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
  // } else {
  //   Matrix<Real> E(spatial_dimension, spatial_dimension);

  //   MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);
  //   this->template gradUToGreenStrain<spatial_dimension>(grad_u, E);

  //   this->computePotentialEnergyOnQuad(E, sigma, *epot);
  //   ++epot;

  //   MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
  // }

  AKANTU_DEBUG_TO_IMPLEMENT();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialIGFEMElastic<spatial_dimension>::computePotentialEnergyByElement(ElementType type, UInt index,
                                                                         Vector<Real> & epot_on_quad_points) {
  // Array<Real>::matrix_iterator gradu_it =
  //   this->gradu(type).begin(spatial_dimension,
  //                            spatial_dimension);
  // Array<Real>::matrix_iterator gradu_end =
  //   this->gradu(type).begin(spatial_dimension,
  //                            spatial_dimension);
  // Array<Real>::matrix_iterator stress_it =
  //   this->stress(type).begin(spatial_dimension,
  //                            spatial_dimension);

  // if (this->finite_deformation)
  //   stress_it = this->piola_kirchhoff_2(type).begin(spatial_dimension,
  //                                                        spatial_dimension);

  // UInt nb_quadrature_points = this->model->getFEEngine().getNbQuadraturePoints(type);

  // gradu_it  += index*nb_quadrature_points;
  // gradu_end += (index+1)*nb_quadrature_points;
  // stress_it  += index*nb_quadrature_points;

  // Real * epot_quad = epot_on_quad_points.storage();

  // Matrix<Real> grad_u(spatial_dimension, spatial_dimension);

  // for(;gradu_it != gradu_end; ++gradu_it, ++stress_it, ++epot_quad) {

  //   if (this->finite_deformation)
  //     this->template gradUToGreenStrain<spatial_dimension>(*gradu_it, grad_u);
  //   else
  //     grad_u.copy(*gradu_it);

  //   this->computePotentialEnergyOnQuad(grad_u, *stress_it, *epot_quad);
  // }
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */


INSTANTIATE_MATERIAL(MaterialIGFEMElastic);

__END_AKANTU__
