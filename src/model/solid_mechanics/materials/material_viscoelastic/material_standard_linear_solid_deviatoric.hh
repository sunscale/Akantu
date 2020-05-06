/**
 * @file   material_standard_linear_solid_deviatoric.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Vladislav Yastrebov <vladislav.yastrebov@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Sun Dec 03 2017
 *
 * @brief  Material Visco-elastic, based on Standard Solid rheological model,
 * see
 * [] J.C.  Simo, T.J.R. Hughes, "Computational  Inelasticity", Springer (1998),
 * see Sections 10.2 and 10.3
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
#include "material_elastic.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_STANDARD_LINEAR_SOLID_DEVIATORIC_HH__
#define __AKANTU_MATERIAL_STANDARD_LINEAR_SOLID_DEVIATORIC_HH__

namespace akantu {

/**
 * Material standard linear solid deviatoric
 *
 *
 * @verbatim

             E_\inf
      ------|\/\/\|------
      |                 |
   ---|                 |---
      |                 |
      ----|\/\/\|--[|----
            E_v   \eta

 @endverbatim
 *
 * keyword : sls_deviatoric
 *
 * parameters in the material files :
 *   - E   : Initial Young's modulus @f$ E = E_i + E_v @f$
 *   - eta : viscosity
 *   - Ev  : stiffness of the viscous element
 */

template <UInt spatial_dimension>
class MaterialStandardLinearSolidDeviatoric
    : public MaterialElastic<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialStandardLinearSolidDeviatoric(SolidMechanicsModel & model,
                                        const ID & id = "");
  ~MaterialStandardLinearSolidDeviatoric() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the material computed parameter
  void initMaterial() override;

  /// update the internal parameters (for modifiable parameters)
  void updateInternalParameters() override;

  /// set material to steady state
  void setToSteadyState(ElementType el_type,
                        GhostType ghost_type = _not_ghost) override;

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type,
                     GhostType ghost_type = _not_ghost) override;

protected:
  /// update the dissipated energy, is called after the stress have been
  /// computed
  void updateDissipatedEnergy(ElementType el_type, GhostType ghost_type);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// give the dissipated energy for the time step
  Real getDissipatedEnergy() const;
  Real getDissipatedEnergy(ElementType type, UInt index) const;

  /// get the energy using an energy type string for the time step
  Real getEnergy(const std::string & type) override;
  Real getEnergy(const std::string & energy_id, ElementType type,
                 UInt index) override;
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// viscosity, viscous elastic modulus
  Real eta, Ev, E_inf;

  Vector<Real> etas;

  /// history of deviatoric stress
  InternalField<Real> stress_dev;

  /// Internal variable: history integral
  InternalField<Real> history_integral;

  /// Dissipated energy
  InternalField<Real> dissipated_energy;
};

} // namespace akantu

#endif /* __AKANTU_MATERIAL_STANDARD_LINEAR_SOLID_DEVIATORIC_HH__ */
