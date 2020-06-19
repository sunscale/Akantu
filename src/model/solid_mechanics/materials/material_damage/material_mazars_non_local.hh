/**
 * @file   material_mazars_non_local.hh
 *
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Mon Sep 11 2017
 *
 * @brief  Mazars non-local description
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
#include "aka_common.hh"
#include "material_damage_non_local.hh"
#include "material_mazars.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_MAZARS_NON_LOCAL_HH__
#define __AKANTU_MATERIAL_MAZARS_NON_LOCAL_HH__

namespace akantu {

/**
 * Material Mazars Non local
 *
 * parameters in the material files :
 */
template <UInt spatial_dimension>
class MaterialMazarsNonLocal
    : public MaterialDamageNonLocal<spatial_dimension,
                                    MaterialMazars<spatial_dimension>> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using MaterialNonLocalParent =
      MaterialDamageNonLocal<spatial_dimension,
                             MaterialMazars<spatial_dimension>>;

  MaterialMazarsNonLocal(SolidMechanicsModel & model, const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  /// constitutive law for all element of a type
  void computeStress(ElementType el_type,
                     GhostType ghost_type = _not_ghost) override;

  void computeNonLocalStress(ElementType el_type,
                             GhostType ghost_type = _not_ghost) override;

  void registerNonLocalVariables() override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// the ehat per quadrature points to perform the averaging
  InternalField<Real> Ehat;

  InternalField<Real> non_local_variable;
};

} // namespace akantu

#endif /* __AKANTU_MATERIAL_MAZARS_NON_LOCAL_HH__ */
