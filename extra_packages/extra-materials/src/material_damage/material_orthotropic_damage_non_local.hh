/**
 * @file   material_orthotropic_damage_non_local.hh
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Sun Mar 22 21:10:27 2015
 *
 * @brief  interface for non local orthotropic damage material
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

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "material_non_local.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_ORTHOTROPIC_DAMAGE_NON_LOCAL_HH_
#define AKANTU_MATERIAL_ORTHOTROPIC_DAMAGE_NON_LOCAL_HH_

namespace akantu {

template <UInt spatial_dimension, class MaterialOrthotropicDamageLocal>
class MaterialOrthotropicDamageNonLocal
    : public MaterialOrthotropicDamageLocal,
      public MaterialNonLocal<spatial_dimension> {
public:
  typedef MaterialNonLocal<spatial_dimension> MaterialNonLocalParent;
  typedef MaterialOrthotropicDamageLocal MaterialOrthotropicDamageParent;

  MaterialOrthotropicDamageNonLocal(SolidMechanicsModel & model, const ID & id)
      : Material(model, id), MaterialOrthotropicDamageParent(model, id),
        MaterialNonLocalParent(model, id){};

  /* ------------------------------------------------------------------------ */
  virtual void initMaterial() {
    MaterialOrthotropicDamageParent::initMaterial();
    MaterialNonLocalParent::initMaterial();
  }

protected:
  /* --------------------------------------------------------------------------
   */
  virtual void computeNonLocalStress(ElementType type,
                                     GhostType ghost_type = _not_ghost) = 0;

  /* ------------------------------------------------------------------------ */
  void computeNonLocalStresses(GhostType ghost_type) {
    AKANTU_DEBUG_IN();

    Mesh::type_iterator it = this->model.getFEEngine().getMesh().firstType(
        spatial_dimension, ghost_type);
    Mesh::type_iterator last_type =
        this->model.getFEEngine().getMesh().lastType(spatial_dimension,
                                                     ghost_type);
    for (; it != last_type; ++it) {
      computeNonLocalStress(*it, ghost_type);
    }

    AKANTU_DEBUG_OUT();
  }

public:
  /* ------------------------------------------------------------------------ */
  virtual inline UInt getNbDataForElements(const Array<Element> & elements,
                                           SynchronizationTag tag) const {
    return MaterialNonLocalParent::getNbDataForElements(elements, tag) +
           MaterialOrthotropicDamageParent::getNbDataForElements(elements, tag);
  }
  virtual inline void packElementData(CommunicationBuffer & buffer,
                                      const Array<Element> & elements,
                                      SynchronizationTag tag) const {
    MaterialNonLocalParent::packElementData(buffer, elements, tag);
    MaterialOrthotropicDamageParent::packElementData(buffer, elements, tag);
  }

  virtual inline void unpackElementData(CommunicationBuffer & buffer,
                                        const Array<Element> & elements,
                                        SynchronizationTag tag) {
    MaterialNonLocalParent::unpackElementData(buffer, elements, tag);
    MaterialOrthotropicDamageParent::unpackElementData(buffer, elements, tag);
  }
};

} // namespace akantu

#endif /* AKANTU_MATERIAL_ORTHOTROPIC_DAMAGE_NON_LOCAL_HH_ */
