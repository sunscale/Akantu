/**
 * @file   material_marigo.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Sun Jul 09 2017
 *
 * @brief  Specialization of the material class for the marigo material
 *
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
#include "material_marigo.hh"
#include "solid_mechanics_model.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
MaterialMarigo<spatial_dimension>::MaterialMarigo(SolidMechanicsModel & model,
                                                  const ID & id)
    : MaterialDamage<spatial_dimension>(model, id), Yd("Yd", *this),
      damage_in_y(false), yc_limit(false) {
  AKANTU_DEBUG_IN();

  this->registerParam("Sd", Sd, Real(5000.), _pat_parsable | _pat_modifiable);
  this->registerParam("epsilon_c", epsilon_c, Real(0.), _pat_parsable,
                      "Critical strain");
  this->registerParam("Yc limit", yc_limit, false, _pat_internal,
                      "As the material a critical Y");
  this->registerParam("damage_in_y", damage_in_y, false, _pat_parsable,
                      "Use threshold (1-D)Y");
  this->registerParam("Yd", Yd, _pat_parsable, "Damaging energy threshold");

  this->Yd.initialize(1);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialMarigo<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialDamage<spatial_dimension>::initMaterial();

  updateInternalParameters();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialMarigo<spatial_dimension>::updateInternalParameters() {
  MaterialDamage<spatial_dimension>::updateInternalParameters();
  Yc = .5 * epsilon_c * this->E * epsilon_c;
  if (std::abs(epsilon_c) > std::numeric_limits<Real>::epsilon())
    yc_limit = true;
  else
    yc_limit = false;
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialMarigo<spatial_dimension>::computeStress(ElementType el_type,
                                                      GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Array<Real>::scalar_iterator dam = this->damage(el_type, ghost_type).begin();
  Array<Real>::scalar_iterator Yd_q = this->Yd(el_type, ghost_type).begin();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  Real Y = 0.;
  computeStressOnQuad(grad_u, sigma, *dam, Y, *Yd_q);

  ++dam;
  ++Yd_q;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

INSTANTIATE_MATERIAL(marigo, MaterialMarigo);

} // namespace akantu
