/**
 * @file   material_selector.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Nov 13 2013
 * @date last modification: Thu Dec 17 2015
 *
 * @brief  class describing how to choose a material for a given element
 *
 * @section LICENSE
 *
 * Copyright  (©)  2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de Lausanne)
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
#include "mesh.hh"

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_SELECTOR_HH__
#define __AKANTU_MATERIAL_SELECTOR_HH__

/* -------------------------------------------------------------------------- */

namespace akantu {

class SolidMechanicsModel;

/**
 * main class to assign same or different materials for different
 * elements
 */
class MaterialSelector {
public:
  MaterialSelector() : fallback_value(0) {}
  virtual ~MaterialSelector() {}
  virtual UInt operator()(__attribute__((unused)) const Element & element) {
    return fallback_value;
  }

  void setFallback(UInt f) { fallback_value = f; }

protected:
  UInt fallback_value;
};

/* -------------------------------------------------------------------------- */
/**
 * class that assigns the first material to regular elements by default
 */
class DefaultMaterialSelector : public MaterialSelector {
public:
  DefaultMaterialSelector(const ElementTypeMapArray<UInt> & material_index)
      : material_index(material_index) {}

  UInt operator()(const Element & element) {
    try {
      DebugLevel dbl = debug::getDebugLevel();
      debug::setDebugLevel(dblError);

      const Array<UInt> & mat_indexes =
          material_index(element.type, element.ghost_type);
      UInt mat = this->fallback_value;

      if (element.element < mat_indexes.getSize())
        mat = mat_indexes(element.element);

      debug::setDebugLevel(dbl);
      return mat;
    } catch (...) {
      return MaterialSelector::operator()(element);
    }
  }

private:
  const ElementTypeMapArray<UInt> & material_index;
};

/* -------------------------------------------------------------------------- */
/**
 * Use elemental data to assign materials
 */
template <typename T>
class ElementDataMaterialSelector : public MaterialSelector {
public:
  ElementDataMaterialSelector(const ElementTypeMapArray<T> & element_data,
                              const SolidMechanicsModel & model,
                              UInt first_index = 1)
      : element_data(element_data), model(model), first_index(first_index) {}

  inline T elementData(const Element & element) {
    DebugLevel dbl = debug::getDebugLevel();
    debug::setDebugLevel(dblError);
    T data = element_data(element.type, element.ghost_type)(element.element);
    debug::setDebugLevel(dbl);
    return data;
  }

  inline UInt operator()(const Element & element) {
    return MaterialSelector::operator()(element);
  }

protected:
  /// list of element with the specified data (i.e. tag value)
  const ElementTypeMapArray<T> & element_data;

  /// the model that the materials belong
  const SolidMechanicsModel & model;

  /// first material index: equal to 1 if none specified
  UInt first_index;
};

/* -------------------------------------------------------------------------- */
/**
 * class to use mesh data information to assign different materials
 * where name is the tag value: tag_0, tag_1
 */
template <typename T>
class MeshDataMaterialSelector : public ElementDataMaterialSelector<T> {
public:
  MeshDataMaterialSelector(const std::string & name,
                           const SolidMechanicsModel & model,
                           UInt first_index = 1);
};

} // akantu

#endif /* __AKANTU_MATERIAL_SELECTOR_HH__ */
