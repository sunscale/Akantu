/**
 * @file   material_vreepeerlings_non_local_inline_impl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Cyprien Wolff <cyprien.wolff@epfl.ch>
 *
 *
 * @brief  Specialization of the material class for the non-local Vree-Peerlings
 * material
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension, template <UInt> class MatParent>
MaterialVreePeerlingsNonLocal<spatial_dimension, MatParent>::
    MaterialVreePeerlingsNonLocal(SolidMechanicsModel & model, const ID & id)
    : Material(model, id), MaterialVreePeerlingsNonLocalParent(model, id),
      equi_strain_non_local("equi-strain_non_local", *this),
      equi_strain_rate_non_local("equi-strain-rate_non_local", *this) {
  AKANTU_DEBUG_IN();

  this->is_non_local = true;

  this->equi_strain_non_local.initialize(1);
  this->equi_strain_rate_non_local.initialize(1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension, template <UInt> class MatParent>
void MaterialVreePeerlingsNonLocal<spatial_dimension,
                                   MatParent>::initMaterial() {
  AKANTU_DEBUG_IN();

  this->registerNonLocalVariable(this->equi_strain, this->equi_strain_non_local,
                                 1);
  this->registerNonLocalVariable(this->equi_strain_rate,
                                 this->equi_strain_rate_non_local, 1);

  MaterialVreePeerlingsNonLocalParent::initMaterial();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

// template<UInt spatial_dimension, class WeigthFunction, template <UInt> class
// MatParent>
// void MaterialVreePeerlingsNonLocal<spatial_dimension, WeigthFunction,
// MatParent>::computeStress(ElementType el_type,
//										     GhostType ghost_type) {
//   AKANTU_DEBUG_IN();
//
//  Real * dam = this->damage(el_type, ghost_type).storage();
//  Real * equi_straint = equi_strain(el_type, ghost_type).storage();
//  Real * equi_straint_rate = equi_strain_rate(el_type, ghost_type).storage();
//  Real * Kapaq = this->Kapa(el_type, ghost_type).storage();
//  Real * crit_strain = this->critical_strain(el_type, ghost_type).storage();
//  Real * crit_strain_rate = this->critical_strain_rate(el_type,
//  ghost_type).storage();
//  Real * rdr_damage = this->recorder_damage(el_type, ghost_type).storage();
//  Real  * nb_damage = this->number_damage(el_type, ghost_type).storage();
//  Real dt = this->model.getTimeStep();
//
//  Vector<UInt> & elem_filter = this->element_filter(el_type, ghost_type);
//  Vector<Real> & velocity = this->model.getVelocity();
//  Vector<Real> & strain_rate_vrplgs = this->strain_rate_vreepeerlings(el_type,
//  ghost_type);
//
//
//  this->model.getFEEngine().gradientOnQuadraturePoints(velocity,
//  strain_rate_vrplgs,
//						   spatial_dimension,
//						   el_type, ghost_type, &elem_filter);
//
//  Vector<Real>::iterator<types::RMatrix> strain_rate_vrplgs_it =
//    strain_rate_vrplgs.begin(spatial_dimension, spatial_dimension);
//
//
//  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);
//
//  types::RMatrix & strain_rate = *strain_rate_vrplgs_it;
//
//
//
//  MaterialVreePeerlings<spatial_dimension>::computeStressOnQuad(grad_u, sigma,
//								*dam,
//								*equi_straint,
//								*equi_straint_rate,
//								*Kapaq,
//								dt,
//								strain_rate,
//								*crit_strain,
//								*crit_strain_rate,
//								*rdr_damage,
//								*nb_damage);
//  ++dam;
//  ++equi_straint;
//  ++equi_straint_rate;
//  ++Kapaq;
//  ++strain_rate_vrplgs_it;
//  ++crit_strain;
//  ++crit_strain_rate;
//  ++rdr_damage;
//  ++nb_damage;
//
//  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
//
//  AKANTU_DEBUG_OUT();
//}
//

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension, template <UInt> class MatParent>
void MaterialVreePeerlingsNonLocal<
    spatial_dimension, MatParent>::computeNonLocalStress(ElementType el_type,
                                                         GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real * dam = this->damage(el_type, ghost_type).storage();
  Real * Kapaq = this->Kapa(el_type, ghost_type).storage();
  Real * equi_strain_nl =
      this->equi_strain_non_local(el_type, ghost_type).storage();
  Real * equi_strain_rate_nl =
      this->equi_strain_rate_non_local(el_type, ghost_type).storage();
  // Real * equi_strain_rate_nl = this->equi_strain_rate(el_type,
  // ghost_type).storage();

  Real dt = this->model.getTimeStep();
  Real * FullDam_Valstrain =
      this->Full_dam_value_strain(el_type, ghost_type).storage();
  Real * FullDam_Valstrain_rate =
      this->Full_dam_value_strain_rate(el_type, ghost_type).storage();
  Real * Nb_damage = this->Number_damage(el_type, ghost_type).storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  this->computeDamageAndStressOnQuad(
      sigma, *dam, *equi_strain_nl, *equi_strain_rate_nl, *Kapaq, dt,
      *FullDam_Valstrain, *FullDam_Valstrain_rate, *Nb_damage);
  ++dam;
  ++equi_strain_nl;
  ++equi_strain_rate_nl;
  ++Kapaq;
  ++FullDam_Valstrain;
  ++FullDam_Valstrain_rate;
  ++Nb_damage;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}
