/**
 * @file   material_anisotropic_damage.cc
 *
 * @author Nicolas Richart
 *
 * @date creation  mer jun 26 2019
 *
 * @brief A Documented file.
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
#include "material_anisotropic_damage.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
/* -------------------------------------------------------------------------- */

namespace {
  template <UInt dim>
  std::unique_ptr<Material>
  materialAnisotropicDamage(std::integral_constant<UInt, dim>,
                            const ID & option, SolidMechanicsModel & model,
                            const ID & id) {
    if (option == "" or option == "mazars") {
      return std::make_unique<MaterialAnisotropicDamage<
          dim, EquivalentStrainMazars, DamageThresholdTan>>(model, id);
    } else if (option == "mazars-drucker-prager") {
      return std::make_unique<MaterialAnisotropicDamage<
          dim, EquivalentStrainMazarsDruckerPrager, DamageThresholdTan>>(model,
                                                                         id);
    } else {
      AKANTU_EXCEPTION("The option "
                       << option << " is not valid for the material " << id);
    }
  }

  template <class... Args>
  decltype(auto) dimensionDispatch(UInt dim, Args &&... args) {
    switch (dim) {
    case 1:
      return materialAnisotropicDamage(std::integral_constant<UInt, 1>{},
                                       std::forward<Args>(args)...);
    case 2:
      return materialAnisotropicDamage(std::integral_constant<UInt, 2>{},
                                       std::forward<Args>(args)...);
    case 3:
      return materialAnisotropicDamage(std::integral_constant<UInt, 3>{},
                                       std::forward<Args>(args)...);
    default: { AKANTU_EXCEPTION("In what dimension are you leaving ?"); }
    }
  }
} // namespace
static bool material_is_alocated_anisotropic_damage [[gnu::unused]] =
    MaterialFactory::getInstance().registerAllocator(
        "anisotropic_damage",
        [](UInt dim, const ID & option, SolidMechanicsModel & model,
           const ID & id) -> std::unique_ptr<Material> {
          return dimensionDispatch(dim, option, model, id);
        });

} // namespace akantu
