
/**
 * @file   material_damage_iterative_non_local.hh
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Thu Mar 27 11:08:37 2014
 *
 * @brief  MaterialDamageIterativeNonLocal header for non-local damage
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
#include "aka_common.hh"
#include "material_damage_iterative.hh"
#include "material_damage_non_local.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_DAMAGE_ITERATIVE_NON_LOCAL_HH__
#define __AKANTU_MATERIAL_DAMAGE_ITERATIVE_NON_LOCAL_HH__

__BEGIN_AKANTU__

/**
 * Material Damage Iterative Non local
 *
 * parameters in the material files :
 */
template<UInt spatial_dimension, template <UInt> class WeightFunction = BaseWeightFunction>
class MaterialDamageIterativeNonLocal : public MaterialDamageNonLocal<spatial_dimension, MaterialDamageIterative<spatial_dimension>, WeightFunction> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  typedef MaterialDamageNonLocal<spatial_dimension, MaterialDamageIterative<spatial_dimension>, WeightFunction> MaterialDamageIterativeNonLocalParent;
  MaterialDamageIterativeNonLocal(SolidMechanicsModel & model, const ID & id = "");

  virtual ~MaterialDamageIterativeNonLocal() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  void initMaterial();

protected:
  void computeStress(ElementType type, GhostType ghost_type);

  void computeNonLocalStress(ElementType type, GhostType ghost_type = _not_ghost);
private:

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  InternalField<Real> grad_u_nl;  
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_damage_iterative_non_local_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_DAMAGE_ITERATIVE_NON_LOCAL_HH__ */
