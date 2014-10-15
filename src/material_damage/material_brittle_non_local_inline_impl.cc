/**
 * @file   material_brittle_non_local_inline_impl.cc
 *
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 *
 * @date   Tue Mar 05 18:18:29 2014
 *
 * @brief  MaterialBrittleNonLocal inline function implementation
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
__END_AKANTU__

#if defined(AKANTU_DEBUG_TOOLS)
#include "aka_debug_tools.hh"
#include <string>
#endif

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeigthFunction>
MaterialBrittleNonLocal<spatial_dimension, WeigthFunction>::MaterialBrittleNonLocal(SolidMechanicsModel & model, const ID & id)  :
  Material(model, id),
  MaterialBrittleNonLocalParent(model, id),
  Sigma_max("Sigma max", *this),
  Sigma_maxnl("Sigma_max non local", *this),
  Sigma_fracture("Sigma_fracture", *this){
  AKANTU_DEBUG_IN();
  this->is_non_local = true;
  this->Sigma_max.initialize(1);
  this->Sigma_maxnl.initialize(1);
  this->Sigma_fracture.initialize(1);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeigthFunction>
void MaterialBrittleNonLocal<spatial_dimension, WeigthFunction>::initMaterial() {
  AKANTU_DEBUG_IN();
  this->registerNonLocalVariable(Sigma_max, Sigma_maxnl, 1);
  MaterialBrittleNonLocalParent::initMaterial();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeigthFunction>
void MaterialBrittleNonLocal<spatial_dimension, WeigthFunction>::computeStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real * dam = this->damage(el_type, ghost_type).storage();
  Real * Sigma_maxt  = this->Sigma_max(el_type, ghost_type).storage();
  Real * fracture_stress  = this->Sigma_fracture(el_type, ghost_type).storage();

  Array<Real> & velocity = this->model->getVelocity();
  Array<Real> & strain_rate_brittle = this->strain_rate_brittle(el_type, ghost_type);
  Array<UInt> & elem_filter = this->element_filter(el_type, ghost_type);

  this->model->getFEEngine().gradientOnQuadraturePoints(velocity, strain_rate_brittle,
						   spatial_dimension,
						   el_type, ghost_type, elem_filter);

  Array<Real>::iterator< Matrix<Real> > strain_rate_it =
    this->strain_rate_brittle(el_type, ghost_type).begin(spatial_dimension, spatial_dimension);

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  Matrix<Real> & grad_v = *strain_rate_it;

  MaterialBrittle<spatial_dimension>::computeStressOnQuad(grad_u, grad_v, sigma, *dam, *Sigma_maxt, *fracture_stress);
  ++dam;
  ++strain_rate_it;
  ++Sigma_maxt;
  ++fracture_stress;
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeigthFunction>
void MaterialBrittleNonLocal<spatial_dimension, WeigthFunction>::computeNonLocalStress(ElementType type,
										      GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real * dam  = this->damage(type, ghost_type).storage();
  Real * Sigma_maxnlt = this->Sigma_maxnl(type, ghost_type).storage();
  Real * fracture_stress = this->Sigma_fracture(type, ghost_type).storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(type, ghost_type);
  this->computeDamageAndStressOnQuad(sigma, *dam, *Sigma_maxnlt, *fracture_stress);

  ++dam;
  ++Sigma_maxnlt;
  ++fracture_stress;
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}
