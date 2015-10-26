/**
 * @file   test_material.hh
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Wed Sep 23 17:16:30 2015
 *
 * @brief  test material for the non-local neighborhood base test
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
#include "material_elastic.hh"
#include "material_non_local.hh"

#ifndef __TEST_MATERIAL_HH__
#define __TEST_MATERIAL_HH__

__BEGIN_AKANTU__

template<UInt dim>
class TestMaterial : public MaterialElastic<dim>,
		     public MaterialNonLocal<dim>{

/* -------------------------------------------------------------------------- */
/* Constructor/Destructor                                                     */
/* -------------------------------------------------------------------------- */

public:
  
  TestMaterial(SolidMechanicsModel & model, const ID & id);
  virtual ~TestMaterial() {};
  typedef MaterialNonLocal<dim> MyNonLocalParent;
  typedef MaterialElastic<dim> MyElasticParent;
/* -------------------------------------------------------------------------- */
/* Methods                                                                    */
/* -------------------------------------------------------------------------- */
public:
  void initMaterial();

  void computeNonLocalStresses(GhostType ghost_type) {};

  void insertQuadsInNeighborhoods(GhostType ghost_type);

  virtual void registerNeighborhood();

protected:

/* -------------------------------------------------------------------------- */
/* Members                                                                   */
/* -------------------------------------------------------------------------- */
private:
  InternalField<Real> grad_u_nl; 

};

__END_AKANTU__

#endif /* __TEST_MATERIAL_HH__ */
