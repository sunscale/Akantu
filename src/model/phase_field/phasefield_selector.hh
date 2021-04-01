/**
 * @file   phasefield_selector.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Sun Mar 1 2020
 * @date last modification: Sun Mar 1 2020
 *
 * @brief class describing how to choose a phasefield variable
 * function for a given element
 *
 * @section LICENSE
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "element.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */
#include <memory>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_PHASEFIELD_SELECTOR_HH__
#define __AKANTU_PHASEFIELD_SELECTOR_HH__

/* -------------------------------------------------------------------------- */
namespace akantu {

class PhaseFieldModel;

/**
 * main class to assign same or different phasefields for different
 * elements
 */
class PhaseFieldSelector : public std::enable_shared_from_this<PhaseFieldSelector> {
public:
  PhaseFieldSelector() = default;
  virtual ~PhaseFieldSelector() = default;
  virtual inline UInt operator()(const Element & element) {
    if (fallback_selector)
      return (*fallback_selector)(element);

    return fallback_value;
  }

  inline void setFallback(UInt f) { fallback_value = f; }
  inline void
  setFallback(const std::shared_ptr<PhaseFieldSelector> & fallback_selector) {
    this->fallback_selector = fallback_selector;
  }

  inline void setFallback(PhaseFieldSelector & fallback_selector) {
    this->fallback_selector = fallback_selector.shared_from_this();
  }

  inline std::shared_ptr<PhaseFieldSelector> & getFallbackSelector() {
    return this->fallback_selector;
  }

  inline UInt getFallbackValue() { return this->fallback_value; }

protected:
  UInt fallback_value{0};
  std::shared_ptr<PhaseFieldSelector> fallback_selector;
};

/* -------------------------------------------------------------------------- */
/**
 * class that assigns the first phasefield to regular elements by default
 */
class DefaultPhaseFieldSelector : public PhaseFieldSelector {
public:
  explicit DefaultPhaseFieldSelector(
      const ElementTypeMapArray<UInt> & phasefield_index)
      : phasefield_index(phasefield_index) {}

  UInt operator()(const Element & element) override {
    if (not phasefield_index.exists(element.type, element.ghost_type))
      return PhaseFieldSelector::operator()(element);

    const auto & phase_indexes = phasefield_index(element.type, element.ghost_type);
    if (element.element < phase_indexes.size()) {
      auto && tmp_phase = phase_indexes(element.element);
      if (tmp_phase != UInt(-1))
        return tmp_phase;
    }

    return PhaseFieldSelector::operator()(element);
  }

private:
  const ElementTypeMapArray<UInt> & phasefield_index;
};

/* -------------------------------------------------------------------------- */
/**
 * Use elemental data to assign phasefields
 */
template <typename T>
class ElementDataPhaseFieldSelector : public PhaseFieldSelector {
public:
  ElementDataPhaseFieldSelector(const ElementTypeMapArray<T> & element_data,
                              const PhaseFieldModel & model,
                              UInt first_index = 1)
      : element_data(element_data), model(model), first_index(first_index) {}

  inline T elementData(const Element & element) {
    DebugLevel dbl = debug::getDebugLevel();
    debug::setDebugLevel(dblError);
    T data = element_data(element.type, element.ghost_type)(element.element);
    debug::setDebugLevel(dbl);
    return data;
  }

  inline UInt operator()(const Element & element) override {
    return PhaseFieldSelector::operator()(element);
  }

protected:
  /// list of element with the specified data (i.e. tag value)
  const ElementTypeMapArray<T> & element_data;

  /// the model that the materials belong
  const PhaseFieldModel & model;

  /// first phasefield index: equal to 1 if none specified
  UInt first_index;
};

/* -------------------------------------------------------------------------- */
/**
 * class to use mesh data information to assign different phasefields
 * where name is the tag value: tag_0, tag_1
 */
template <typename T>
class MeshDataPhaseFieldSelector : public ElementDataPhaseFieldSelector<T> {
public:
  MeshDataPhaseFieldSelector(const std::string & name,
                           const PhaseFieldModel & model,
                           UInt first_index = 1);
};

} // namespace akantu

#endif /* __AKANTU_PHASEFIELD_SELECTOR_HH__ */
