/**
 * @file   material_damage_iterative_non_local_inline_impl.cc
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Thu Mar 27 12:13:18 2014
 *
 * @brief  MaterialDamageIterativeNonLocal inline function implementation
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
MaterialDamageIterativeNonLocal<spatial_dimension, WeigthFunction>::MaterialDamageIterativeNonLocal(SolidMechanicsModel & model, const ID & id)  :
  Material(model, id),
  MaterialDamageIterativeNonLocalParent(model, id),
  grad_u_nl("grad_u non local", *this) {
  AKANTU_DEBUG_IN();
  this->is_non_local = true;
  this->grad_u_nl.initialize(spatial_dimension*spatial_dimension);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeigthFunction>
void MaterialDamageIterativeNonLocal<spatial_dimension, WeigthFunction>::initMaterial() {
  AKANTU_DEBUG_IN();
  this->registerNonLocalVariable(this->gradu, grad_u_nl, spatial_dimension*spatial_dimension);
  MaterialDamageIterativeNonLocalParent::initMaterial();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeigthFunction>
void MaterialDamageIterativeNonLocal<spatial_dimension, WeigthFunction>::computeStress(ElementType type,
										       GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension, template <UInt> class WeigthFunction>
void MaterialDamageIterativeNonLocal<spatial_dimension, WeigthFunction>::computeNonLocalStress(ElementType el_type,
											       GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  
  MaterialDamage<spatial_dimension>::computeStress(el_type, ghost_type);

  Real * dam = this->damage(el_type, ghost_type).storage();
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  this->computeDamageAndStressOnQuad(sigma,*dam);
  ++dam;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  this->computeNormalizedEquivalentStress(this->grad_u_nl(el_type, ghost_type), el_type, ghost_type);
  this->norm_max_equivalent_stress = 0;
  this->findMaxNormalizedEquivalentStress(el_type, ghost_type);

  AKANTU_DEBUG_OUT();
}
