/**
 * @file   dumper_quadrature_point_iterator.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  Description of quadrature point iterator
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef __AKANTU_DUMPER_QUADRATURE_POINT_ITERATOR_HH__
#define __AKANTU_DUMPER_QUADRATURE_POINT_ITERATOR_HH__

/* -------------------------------------------------------------------------- */
#include "dumper_elemental_field.hh"

namespace akantu {
namespace dumpers {

/* -------------------------------------------------------------------------- */
template <typename types>
class quadrature_point_iterator
    : public element_iterator<types, quadrature_point_iterator> {
  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using parent = element_iterator<types, dumpers::quadrature_point_iterator>;
  using data_type = typename types::data_type;
  using return_type = typename types::return_type;
  using field_type = typename types::field_type;
  using array_iterator = typename types::array_iterator;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  quadrature_point_iterator(const field_type & field,
                            const typename field_type::type_iterator & t_it,
                            const typename field_type::type_iterator & t_it_end,
                            const array_iterator & array_it,
                            const array_iterator & array_it_end,
                            const GhostType ghost_type = _not_ghost)
      : parent(field, t_it, t_it_end, array_it, array_it_end, ghost_type) {}

  return_type operator*() { return *this->array_it; }
};

/* -------------------------------------------------------------------------- */

} // namespace dumpers
} // namespace akantu

#endif /* __AKANTU_DUMPER_QUADRATURE_POINT_ITERATOR_HH__ */
