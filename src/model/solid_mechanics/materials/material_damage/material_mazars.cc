/**
 * @file   material_mazars.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Sun Jul 09 2017
 *
 * @brief  Specialization of the material class for the damage material
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
#include "material_mazars.hh"
#include "solid_mechanics_model.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
MaterialMazars<spatial_dimension>::MaterialMazars(SolidMechanicsModel & model,
                                                  const ID & id)
    : MaterialDamage<spatial_dimension>(model, id), K0("K0", *this),
      damage_in_compute_stress(true) {
  AKANTU_DEBUG_IN();

  this->registerParam("K0", K0, _pat_parsable, "K0");
  this->registerParam("At", At, Real(0.8), _pat_parsable, "At");
  this->registerParam("Ac", Ac, Real(1.4), _pat_parsable, "Ac");
  this->registerParam("Bc", Bc, Real(1900.), _pat_parsable, "Bc");
  this->registerParam("Bt", Bt, Real(12000.), _pat_parsable, "Bt");
  this->registerParam("beta", beta, Real(1.06), _pat_parsable, "beta");

  this->K0.initialize(1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialMazars<spatial_dimension>::computeStress(ElementType el_type,
                                                      GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real * dam = this->damage(el_type, ghost_type).storage();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  Real Ehat = 0;
  computeStressOnQuad(grad_u, sigma, *dam, Ehat);
  ++dam;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */

INSTANTIATE_MATERIAL(mazars, MaterialMazars);

} // namespace akantu
