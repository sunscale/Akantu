/**
 * @file   synchronized_array_inline_impl.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Tue Dec 02 2014
 * @date last modification: Fri Feb 23 2018
 *
 * @brief  inlined methods for the synchronized array
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "synchronized_array.hh"
/* -------------------------------------------------------------------------- */
namespace akantu {
/* -------------------------------------------------------------------------- */
template <typename T>
template <typename P>
inline void SynchronizedArray<T>::push_back(P && value) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(deleted_elements.size() == 0,
                      "Cannot push_back element if SynchronizedArray"
                          << " is already modified without synchronization");

  Array<T>::push_back(std::forward<P>(value));
  this->nb_added_elements++;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <typename T> inline void SynchronizedArray<T>::erase(UInt i) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(nb_added_elements == 0,
                      "Cannot erase element if SynchronizedArray"
                          << " is already modified without synchronization");

  for (UInt j = 0; j < this->nb_component; ++j)
    this->values[i * this->nb_component + j] =
        this->values[(this->size_ - 1) * this->nb_component + j];
  this->size_--;

  this->deleted_elements.push_back(i);

  AKANTU_DEBUG_OUT();
}

}
