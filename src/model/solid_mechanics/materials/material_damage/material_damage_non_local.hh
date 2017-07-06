/**
 * @file   material_damage_non_local.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Aug 23 2012
 * @date last modification: Wed Oct 21 2015
 *
 * @brief  interface for non local damage material
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
#include "material_non_local.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_DAMAGE_NON_LOCAL_HH__
#define __AKANTU_MATERIAL_DAMAGE_NON_LOCAL_HH__

namespace akantu {

template <UInt dim, class MaterialDamageLocal>
class MaterialDamageNonLocal : public MaterialNonLocal<dim, MaterialDamageLocal> {
public:
  using MaterialParent = MaterialNonLocal<dim, MaterialDamageLocal>;

  MaterialDamageNonLocal(SolidMechanicsModel & model, const ID & id)
      : MaterialParent(model, id) {};

  /* ------------------------------------------------------------------------ */
  virtual void initMaterial() {
    MaterialParent::initMaterial();
  }

protected:
  /* ------------------------------------------------------------------------ */
  virtual void computeNonLocalStress(ElementType type,
                                     GhostType ghost_type = _not_ghost) = 0;

  /* ------------------------------------------------------------------------ */
  void computeNonLocalStresses(GhostType ghost_type) override {
    AKANTU_DEBUG_IN();

    for (auto type :
         this->element_filter.elementTypes(dim, ghost_type)) {
      auto & elem_filter = this->element_filter(type, ghost_type);
      if (elem_filter.getSize() == 0)
        continue;

      computeNonLocalStress(type, ghost_type);
    }

    AKANTU_DEBUG_OUT();
  }

public:
  /* ------------------------------------------------------------------------ */
  inline UInt getNbData(const Array<Element> & elements,
                        const SynchronizationTag & tag) const override {
    return MaterialParent::getNbData(elements, tag);
  }
  inline void packData(CommunicationBuffer & buffer,
                       const Array<Element> & elements,
                       const SynchronizationTag & tag) const override {
    MaterialParent::packData(buffer, elements, tag);
  }

  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<Element> & elements,
                         const SynchronizationTag & tag) {
    MaterialParent::unpackData(buffer, elements, tag);
  }
};

} // namespace akantu

#endif /* __AKANTU_MATERIAL_DAMAGE_NON_LOCAL_HH__ */
