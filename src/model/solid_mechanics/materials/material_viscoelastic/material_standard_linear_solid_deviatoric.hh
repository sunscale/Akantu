/**
 * @file   material_standard_linear_solid_deviatoric.hh
 *
 * @author Vladislav Yastrebov <vladislav.yastrebov@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Wed Feb 08 2012
 * @date last modification: Wed Nov 13 2013
 *
 * @brief  Material Visco-elastic, based on Standard Solid rheological model, see
 * [] J.C.  Simo, T.J.R. Hughes, "Computational  Inelasticity", Springer (1998),
 * see Sections 10.2 and 10.3
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_common.hh"
#include "material_elastic.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_STANDARD_LINEAR_SOLID_DEVIATORIC_HH__
#define __AKANTU_MATERIAL_STANDARD_LINEAR_SOLID_DEVIATORIC_HH__

__BEGIN_AKANTU__

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

template<UInt spatial_dimension>
class MaterialStandardLinearSolidDeviatoric : public MaterialElastic<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialStandardLinearSolidDeviatoric(SolidMechanicsModel & model, const ID & id = "");

  virtual ~MaterialStandardLinearSolidDeviatoric() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  void initMaterial();

  virtual void updateInternalParameters();

  void setToSteadyState(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);


protected:

  /// update the dissipated energy, is called after the stress have been computed
  void updateDissipatedEnergy(ElementType el_type, GhostType ghost_type);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// give the dissipated energy for the time step
  Real getDissipatedEnergy() const;
  Real getDissipatedEnergy(ElementType type, UInt index) const;

  virtual Real getEnergy(std::string type);
  virtual Real getEnergy(std::string energy_id, ElementType type, UInt index);
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// viscosity, viscous elastic modulus
  Real eta, Ev, E_inf;

  /// history of deviatoric stress
  InternalField<Real> stress_dev;

  /// Internal variable: history integral
  InternalField<Real> history_integral;

  /// Dissipated energy
  InternalField<Real> dissipated_energy;
};

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_STANDARD_LINEAR_SOLID_DEVIATORIC_HH__ */


