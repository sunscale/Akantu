/**
 * @file   material_anisotropic_damage.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  mar jun 25 2019
 *
 * @brief A Documented file.
 *
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
#include "material_elastic.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_ANISOTROPIC_DAMAGE_HH__
#define __AKANTU_MATERIAL_ANISOTROPIC_DAMAGE_HH__

namespace akantu {

template <UInt dim, template <UInt> class EquivalentStrain,
          template <UInt> class DamageThreshold,
          template <UInt> class Parent = MaterialElastic>
class MaterialAnisotropicDamage : public Parent<dim> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialAnisotropicDamage(SolidMechanicsModel & model, const ID & id = "");
  ~MaterialAnisotropicDamage() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void computeStress(ElementType el_type, GhostType ghost_type) override;

private:
  void damageStress(Matrix<double> & sigma, const Matrix<double> & sigma_el,
                    const Matrix<double> & D, Real TrD);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  Real Dc{0.99};

  /// damage internal variable
  InternalField<Real> damage;

  /// elastic stress
  InternalField<Real> elastic_stress;

  /// equivalent strain
  InternalField<Real> equivalent_strain;

  /// trace of the damageThreshold
  InternalField<Real> trace_damage;

  /// damage criteria
  EquivalentStrain<dim> equivalent_strain_function;

  /// damage evolution
  DamageThreshold<dim> damage_threshold_function;
};

} // namespace akantu

#include "material_anisotropic_damage_tmpl.hh"

#endif /* __AKANTU_MATERIAL_ANISOTROPIC_DAMAGE_HH__ */
