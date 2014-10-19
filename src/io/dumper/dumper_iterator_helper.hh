/**
 * @file   dumper_iterator_helper.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Sep 02 2014
 * @date last modification: Tue Sep 02 2014
 *
 * @brief  Helper to write field iterators
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

#ifndef __AKANTU_DUMPER_ITERATOR_HELPER_HH__
#define __AKANTU_DUMPER_ITERATOR_HELPER_HH__
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__
__BEGIN_AKANTU_DUMPER__

template<class T, class R>
class iterator_helper {

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */
public:


  typedef typename Array<T>::template const_iterator< R > internal_iterator;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

public:

  static internal_iterator begin(const Array<T> & vect, UInt m, UInt n, UInt size) {
    return vect.begin_reinterpret(n*m, size);
  }

  static internal_iterator end(const Array<T> & vect, UInt m, UInt n, UInt size) {
    return vect.end_reinterpret(n*m, size);
  }
};

/* -------------------------------------------------------------------------- */

template<class T>
class iterator_helper<T, Matrix<T> > {

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */

public:

  typedef typename Array<T>::template const_iterator< Matrix<T> > internal_iterator;

public:

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */


  static internal_iterator begin(const Array<T> & vect, UInt m, UInt n, UInt size) {
    return vect.begin_reinterpret(m, n, size);
  }

  static internal_iterator end(const Array<T> & vect, UInt m, UInt n, UInt size) {
    return vect.end_reinterpret(m, n, size);
  }
};


/* -------------------------------------------------------------------------- */

__END_AKANTU_DUMPER__
__END_AKANTU__


#endif /* __AKANTU_DUMPER_ITERATOR_HELPER_HH__ */
