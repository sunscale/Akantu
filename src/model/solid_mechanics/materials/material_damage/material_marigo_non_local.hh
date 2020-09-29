/**
 * @file   material_marigo_non_local.hh
 *
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  Marigo non-local description
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
#include "material_marigo.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_MARIGO_NON_LOCAL_HH_
#define AKANTU_MATERIAL_MARIGO_NON_LOCAL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */

/**
 * Material Marigo
 *
 * parameters in the material files :
 */
template <UInt spatial_dimension>
class MaterialMarigoNonLocal
    : public MaterialDamageNonLocal<spatial_dimension,
                                    MaterialMarigo<spatial_dimension>> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using MaterialMarigoNonLocalParent = MaterialDamageNonLocal<spatial_dimension,
                                 MaterialMarigo<spatial_dimension>>;
  MaterialMarigoNonLocal(SolidMechanicsModel & model, const ID & id = "");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  void registerNonLocalVariables() override;

  /// constitutive law
  void computeStress(ElementType el_type,
                     GhostType ghost_type = _not_ghost) override;

  void computeNonLocalStress(ElementType type,
                             GhostType ghost_type = _not_ghost) override;

private:
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Y, Y, Real);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  InternalField<Real> Y;
  InternalField<Real> Ynl;
};

} // namespace akantu

#endif /* AKANTU_MATERIAL_MARIGO_NON_LOCAL_HH_ */
