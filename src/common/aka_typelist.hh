/**
 * @file   aka_typelist.hh
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date creation: Fri Jan 04 2013
 * @date last modification: Mon Jun 19 2017
 *
 * @brief  Objects that support the visitor design pattern
 *
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

#ifndef AKANTU_TYPELIST_HH_
#define AKANTU_TYPELIST_HH_

#include "aka_common.hh"

namespace akantu {

struct Empty_type {};
class Null_type {};

template <class T, class U> struct Typelist {
  typedef T Head;
  typedef U Tail;
};

template <typename T1 = Null_type, typename T2 = Null_type,
          typename T3 = Null_type, typename T4 = Null_type,
          typename T5 = Null_type, typename T6 = Null_type,
          typename T7 = Null_type, typename T8 = Null_type,
          typename T9 = Null_type, typename T10 = Null_type,
          typename T11 = Null_type, typename T12 = Null_type,
          typename T13 = Null_type, typename T14 = Null_type,
          typename T15 = Null_type, typename T16 = Null_type,
          typename T17 = Null_type, typename T18 = Null_type,
          typename T19 = Null_type, typename T20 = Null_type>
struct MakeTypelist {
private:
  typedef typename MakeTypelist<T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12,
                                T13, T14, T15, T16, T17, T18, T19, T20>::Result
      TailResult;

public:
  typedef Typelist<T1, TailResult> Result;
};

template <> struct MakeTypelist<> { typedef Null_type Result; };

////////////////////////////////////////////////////////////////////////////////
// class template Length
// Computes the length of a typelist
// Invocation (TList is a typelist):
// Length<TList>::value
// returns a compile-time constant containing the length of TList, not counting
//     the end terminator (which by convention is Null_type)
////////////////////////////////////////////////////////////////////////////////

template <class TList> struct Length;
template <> struct Length<Null_type> {
  enum { value = 0 };
};

template <class T, class U> struct Length<Typelist<T, U>> {
  enum { value = 1 + Length<U>::value };
};

////////////////////////////////////////////////////////////////////////////////
// class template TypeAt
// Finds the type at a given index in a typelist
// Invocation (TList is a typelist and index is a compile-time integral
//     constant):
// TypeAt<TList, index>::Result
// returns the type in position 'index' in TList
// If you pass an out-of-bounds index, the result is a compile-time error
////////////////////////////////////////////////////////////////////////////////

template <class TList, unsigned int index> struct TypeAt;

template <class Head, class Tail> struct TypeAt<Typelist<Head, Tail>, 0> {
  typedef Head Result;
};

template <class Head, class Tail, unsigned int i>
struct TypeAt<Typelist<Head, Tail>, i> {
  typedef typename TypeAt<Tail, i - 1>::Result Result;
};

////////////////////////////////////////////////////////////////////////////////
// class template Erase
// Erases the first occurence, if any, of a type in a typelist
// Invocation (TList is a typelist and T is a type):
// Erase<TList, T>::Result
// returns a typelist that is TList without the first occurence of T
////////////////////////////////////////////////////////////////////////////////

template <class TList, class T> struct Erase;

template <class T> // Specialization 1
struct Erase<Null_type, T> {
  typedef Null_type Result;
};

template <class T, class Tail> // Specialization 2
struct Erase<Typelist<T, Tail>, T> {
  typedef Tail Result;
};

template <class Head, class Tail, class T> // Specialization 3
struct Erase<Typelist<Head, Tail>, T> {
  typedef Typelist<Head, typename Erase<Tail, T>::Result> Result;
};

template <class TList, class T> struct IndexOf;

template <class T> struct IndexOf<Null_type, T> {
  enum { value = -1 };
};

template <class T, class Tail> struct IndexOf<Typelist<T, Tail>, T> {
  enum { value = 0 };
};

template <class Head, class Tail, class T>
struct IndexOf<Typelist<Head, Tail>, T> {
private:
  enum { temp = IndexOf<Tail, T>::value };

public:
  enum { value = (temp == -1 ? -1 : 1 + temp) };
};

} // namespace akantu

#endif /* AKANTU_TYPELIST_HH_ */
