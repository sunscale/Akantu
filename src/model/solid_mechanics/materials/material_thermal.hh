/**
 * @file   material_thermal.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Aug 18 2015
 *
 * @brief  Material isotropic thermo-elastic
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "material.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_THERMAL_HH__
#define __AKANTU_MATERIAL_THERMAL_HH__

namespace akantu {
template<UInt spatial_dimension>
class MaterialThermal : public virtual Material {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialThermal(SolidMechanicsModel & model, const ID & id = "");
  MaterialThermal(SolidMechanicsModel & model,
                  UInt dim,
                  const Mesh & mesh,
                  FEEngine & fe_engine,
                  const ID & id = "");

  virtual ~MaterialThermal() {};

protected:
  void initialize();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void initMaterial();

  /// constitutive law for all element of a type
  virtual void computeStress(ElementType el_type, GhostType ghost_type);

  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:


  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// Young modulus
  Real E;

  /// Poisson ratio
  Real nu;

  /// Thermal expansion coefficient
  /// TODO : implement alpha as a matrix
  Real alpha;

  /// Temperature field
  InternalField<Real> delta_T;

  /// Current thermal stress
  InternalField<Real> sigma_th;

  /// Tell if we need to use the previous thermal stress
  bool use_previous_stress_thermal;
};

} // akantu

#endif /* __AKANTU_MATERIAL_THERMAL_HH__ */
