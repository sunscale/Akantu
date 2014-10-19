/**
 * @file   material_damage_non_local.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Aug 23 2012
 * @date last modification: Thu Jun 05 2014
 *
 * @brief  interface for non local damage material
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
#include "material_non_local.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_DAMAGE_NON_LOCAL_HH__
#define __AKANTU_MATERIAL_DAMAGE_NON_LOCAL_HH__

__BEGIN_AKANTU__

template<UInt spatial_dimension,
         class MaterialDamageLocal,
	 template <UInt> class WeightFunction = BaseWeightFunction>
class MaterialDamageNonLocal : public MaterialDamageLocal,
			       public MaterialNonLocal<spatial_dimension, WeightFunction> {
public:
  typedef MaterialNonLocal<spatial_dimension, WeightFunction> MaterialNonLocalParent;
  typedef MaterialDamageLocal MaterialDamageParent;

  MaterialDamageNonLocal(SolidMechanicsModel & model, const ID & id)  :
    Material(model, id),
    MaterialDamageParent(model, id), MaterialNonLocalParent(model, id) { };

  /* ------------------------------------------------------------------------ */
  virtual void initMaterial() {
    MaterialDamageParent::initMaterial();
    MaterialNonLocalParent::initMaterial();
  }

protected:
  /* -------------------------------------------------------------------------- */
  virtual void computeNonLocalStress(ElementType type, GhostType ghost_type = _not_ghost) = 0;

  /* ------------------------------------------------------------------------ */
  void computeNonLocalStresses(GhostType ghost_type) {
    AKANTU_DEBUG_IN();

    Mesh::type_iterator it = this->model->getFEEngine().getMesh().firstType(spatial_dimension, ghost_type);
    Mesh::type_iterator last_type = this->model->getFEEngine().getMesh().lastType(spatial_dimension, ghost_type);
    for(; it != last_type; ++it) {
      computeNonLocalStress(*it, ghost_type);
    }

    AKANTU_DEBUG_OUT();
  }

public:
  /* ------------------------------------------------------------------------ */
  virtual inline UInt getNbDataForElements(const Array<Element> & elements,
					   SynchronizationTag tag) const {
    return MaterialNonLocalParent::getNbDataForElements(elements, tag) +
      MaterialDamageParent::getNbDataForElements(elements, tag);
  }
  virtual inline void packElementData(CommunicationBuffer & buffer,
				      const Array<Element> & elements,
				      SynchronizationTag tag) const {
    MaterialNonLocalParent::packElementData(buffer, elements, tag);
    MaterialDamageParent::packElementData(buffer, elements, tag);
  }

  virtual inline void unpackElementData(CommunicationBuffer & buffer,
				 const Array<Element> & elements,
				 SynchronizationTag tag) {
    MaterialNonLocalParent::unpackElementData(buffer, elements, tag);
    MaterialDamageParent::unpackElementData(buffer, elements, tag);
  }

};

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_DAMAGE_NON_LOCAL_HH__ */
