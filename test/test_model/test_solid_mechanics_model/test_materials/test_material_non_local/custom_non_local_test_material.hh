/**
 * @file   custom_non_local_test_material.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Aug 23 2012
 * @date last modification: Thu Oct 15 2015
 *
 * @brief  Custom material to test the non local implementation
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

#include "material_non_local.hh"
#include "material_elastic.hh"

#ifndef __CUSTOM_NON_LOCAL_TEST_MATERIAL_HH__
#define __CUSTOM_NON_LOCAL_TEST_MATERIAL_HH__

__BEGIN_AKANTU__

template<UInt dim>
class CustomNonLocalTestMaterial : public MaterialElastic<dim>,
                                   public MaterialNonLocal<dim> {
public:
  typedef MaterialNonLocal<dim> MyNonLocalParent;
  typedef MaterialElastic<dim> MyElasticParent;

  CustomNonLocalTestMaterial(SolidMechanicsModel & model, const ID & id);

  /* ------------------------------------------------------------------------ */
  virtual void initMaterial();

  void computeNonLocalStress(ElementType el_type, GhostType ghost_type);
  void computeStress(ElementType el_type, GhostType ghost_type);

protected:

  /* ------------------------------------------------------------------------ */
  void computeNonLocalStresses(GhostType ghost_type) {
    AKANTU_DEBUG_IN();

    Mesh::type_iterator it = this->model->getFEEngine().getMesh().firstType(dim, ghost_type);
    Mesh::type_iterator last_type = this->model->getFEEngine().getMesh().lastType(dim, ghost_type);
    for(; it != last_type; ++it) {
      computeNonLocalStress(*it, ghost_type);
    }

    AKANTU_DEBUG_OUT();
  }

public:
  /* ------------------------------------------------------------------------ */
  virtual inline UInt getNbDataForElements(const Array<Element> & elements,
					   SynchronizationTag tag) const {
    return MyNonLocalParent::getNbDataForElements(elements, tag) +
      MyElasticParent::getNbDataForElements(elements, tag);
  }
  virtual inline void packElementData(CommunicationBuffer & buffer,
				      const Array<Element> & elements,
				      SynchronizationTag tag) const {
    MyNonLocalParent::packElementData(buffer, elements, tag);
    MyElasticParent::packElementData(buffer, elements, tag);
  }

  virtual inline void unpackElementData(CommunicationBuffer & buffer,
				 const Array<Element> & elements,
				 SynchronizationTag tag) {
    MyNonLocalParent::unpackElementData(buffer, elements, tag);
    MyElasticParent::unpackElementData(buffer, elements, tag);
  }


  void setDamage(Real dam) {
    this->local_damage.setDefaultValue(dam);
  }
protected:
  InternalField<Real> local_damage;
  InternalField<Real> damage;
};

__END_AKANTU__

#endif /* __CUSTOM_NON_LOCAL_TEST_MATERIAL_HH__ */
