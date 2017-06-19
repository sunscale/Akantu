/**
 * @file   dumper_quadrature_point_iterator.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Sep 02 2014
 * @date last modification: Thu Sep 17 2015
 *
 * @brief  Description of quadrature point iterator
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

#ifndef __AKANTU_DUMPER_QUADRATURE_POINT_ITERATOR_HH__
#define __AKANTU_DUMPER_QUADRATURE_POINT_ITERATOR_HH__

/* -------------------------------------------------------------------------- */
#include "dumper_elemental_field.hh"

namespace akantu {
__BEGIN_AKANTU_DUMPER__

/* -------------------------------------------------------------------------- */
template<typename types>
class quadrature_point_iterator
   : public element_iterator<types, quadrature_point_iterator> {
  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */
public:

  typedef element_iterator<types, dumper::quadrature_point_iterator> parent;
  typedef typename types::data_type   data_type;
  typedef typename types::return_type return_type;
  typedef typename types::field_type  field_type;
  typedef typename types::array_iterator array_iterator;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  quadrature_point_iterator(const field_type & field,
			    const typename field_type::type_iterator & t_it,
			    const typename field_type::type_iterator & t_it_end,
			    const array_iterator & array_it,
			    const array_iterator & array_it_end,
			    const GhostType ghost_type = _not_ghost) :
    parent(field, t_it, t_it_end, array_it, array_it_end, ghost_type) { }

  return_type operator*() {
    return *this->array_it;
  }
};

/* -------------------------------------------------------------------------- */

__END_AKANTU_DUMPER__
} // akantu

#endif /* __AKANTU_DUMPER_QUADRATURE_POINT_ITERATOR_HH__ */
