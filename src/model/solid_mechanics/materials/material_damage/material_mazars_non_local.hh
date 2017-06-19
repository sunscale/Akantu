/**
 * @file   material_mazars_non_local.hh
 *
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Thu Oct 08 2015
 *
 * @brief  Mazars non-local description
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
#include "material_mazars.hh"
#include "material_non_local.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_MAZARS_NON_LOCAL_HH__
#define __AKANTU_MATERIAL_MAZARS_NON_LOCAL_HH__

namespace akantu {

/**
 * Material Mazars Non local
 *
 * parameters in the material files :
 */
template<UInt spatial_dimension>
class MaterialMazarsNonLocal : public MaterialMazars<spatial_dimension>,
			       public MaterialNonLocal<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  typedef MaterialNonLocal<spatial_dimension> MaterialNonLocalParent;

  MaterialMazarsNonLocal(SolidMechanicsModel & model, const ID & id = "");

  virtual ~MaterialMazarsNonLocal() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  void initMaterial();

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// constitutive law
  void computeNonLocalStresses(GhostType ghost_type = _not_ghost);
  void computeNonLocalStress(Array<Real> & Ehatnl,
			     ElementType el_type,
			     GhostType ghost_type = _not_ghost);

protected:
  /// associate the non-local variables of the material to their neighborhoods
  virtual void nonLocalVariableToNeighborhood();

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
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "material_mazars_non_local_inline_impl.cc"

} // akantu

#endif /* __AKANTU_MATERIAL_MAZARS_NON_LOCAL_HH__ */
