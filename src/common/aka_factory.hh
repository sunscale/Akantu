/**
 * @file   aka_factory.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sun Jul 09 2017
 * @date last modification: Fri Dec 08 2017
 *
 * @brief  This is a generic factory
 *
 *
 * Copyright (©) 2016-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */
#include <functional>
#include <map>
#include <memory>
#include <string>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_AKA_FACTORY_HH_
#define AKANTU_AKA_FACTORY_HH_

namespace akantu {

template <class Base, class T = ID, class... Args> class Factory {
  using allocator_t = std::function<std::unique_ptr<Base>(Args...)>;

private:
  Factory() = default;

public:
  Factory(const Factory &) = delete;
  Factory & operator=(const Factory &) = delete;

  static Factory & getInstance() {
    static Factory instance;
    return instance;
  }
  /* ------------------------------------------------------------------------ */
  bool registerAllocator(const T & id, const allocator_t & allocator) {
    if (allocators.find(id) != allocators.end()) {
      AKANTU_EXCEPTION("The id "
                       << " is already registered in the "
                       << debug::demangle(typeid(Base).name()) << " factory");
    }
    allocators[id] = allocator;
    return true;
  }

  template <typename... AArgs>
  std::unique_ptr<Base> allocate(const T & id, AArgs &&... args) const {
    if (allocators.find(id) == allocators.end()) {
      AKANTU_EXCEPTION("The id "
                       << " is not registered in the "
                       << debug::demangle(typeid(Base).name()) << " factory.");
    }
    return std::forward<std::unique_ptr<Base>>(
        allocators.at(id)(std::forward<AArgs>(args)...));
  }

  std::vector<T> getPossibleAllocators() {
    std::vector<T> keys;
    for (auto & e : allocators) {
      keys.push_back(e.first);
    }
    return keys;
  }

private:
  std::map<T, allocator_t> allocators;
};

} // namespace akantu

#endif /* AKANTU_AKA_FACTORY_HH_ */
