/**
 * @file   test_material_damage.hh
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Wed Sep 23 17:16:30 2015
 *
 * @brief  test material damage for the non-local remove damage test
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
#include "material_damage.hh"
#include "material_non_local.hh"

#ifndef __TEST_MATERIAL_DAMAGE_HH__
#define __TEST_MATERIAL_DAMAGE_HH__

__BEGIN_AKANTU__

template<UInt dim>
class TestMaterialDamage : public  MaterialDamage<dim, MaterialElastic>,
									 public MaterialNonLocal<dim> {

/* -------------------------------------------------------------------------- */
/* Constructor/Destructor                                                     */
/* -------------------------------------------------------------------------- */

public:
  
  TestMaterialDamage(SolidMechanicsModel & model, const ID & id);
  virtual ~TestMaterialDamage() {};
  typedef MaterialNonLocal<dim> MyNonLocalParent;

/* -------------------------------------------------------------------------- */
/* Methods                                                                    */
/* -------------------------------------------------------------------------- */
public:
  void initMaterial();

  //void computeNonLocalStress(ElementType type, GhostType ghost_type = _not_ghost);

  void computeNonLocalStresses(GhostType ghost_type) {};

  void insertQuadsInNeighborhoods(GhostType ghost_type);

protected:
  /// associate the non-local variables of the material to their neighborhoods
  virtual void nonLocalVariableToNeighborhood();

/* -------------------------------------------------------------------------- */
/* Members                                                                   */
/* -------------------------------------------------------------------------- */
private:
  InternalField<Real> grad_u_nl; 

};

__END_AKANTU__

#endif /* __TEST_MATERIAL_DAMAGE_HH__ */
