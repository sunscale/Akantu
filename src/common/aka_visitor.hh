/**
 * @file   aka_visitor.hh
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

#ifndef __AKANTU_VISITOR_HH__
#define __AKANTU_VISITOR_HH__

#include "aka_typelist.hh"

namespace akantu {

///////////////////////////////////////////////////////////////////////////
// visitor class template, adapted from the Andrei Alexandrescu's
// "Modern C++ Design"

enum Visit_type { Mutable, Immutable };

template <class T, typename R = void, Visit_type = Mutable> class StrictVisitor;

template <class T, typename R> class StrictVisitor<T, R, Mutable> {
public:
  typedef R ReturnType;
  typedef T ParamType;
  virtual ~StrictVisitor() {}
  virtual ReturnType Visit(ParamType &) = 0;
};

template <class T, typename R> class StrictVisitor<T, R, Immutable> {
public:
  typedef R ReturnType;
  typedef const T ParamType;
  virtual ~StrictVisitor() {}
  virtual ReturnType Visit(ParamType &) = 0;
};

/// class template StrictVisitor (specialization)

template <class Head, class Tail, typename R>
class StrictVisitor<Typelist<Head, Tail>, R, Mutable>
    : public StrictVisitor<Head, R, Mutable>,
      public StrictVisitor<Tail, R, Mutable> {
public:
  typedef R ReturnType;
  typedef Head ParamType;
  //	using StrictVisitor<Head, R>::Visit;
  //	using StrictVisitor<Tail, R>::Visit;
};

template <class Head, typename R>
class StrictVisitor<Typelist<Head, Null_type>, R, Mutable>
    : public StrictVisitor<Head, R, Mutable> {
public:
  typedef R ReturnType;
  typedef Head ParamType;
  using StrictVisitor<Head, R, Mutable>::Visit;
};

template <class Head, class Tail, typename R>
class StrictVisitor<Typelist<Head, Tail>, R, Immutable>
    : public StrictVisitor<Head, R, Immutable>,
      public StrictVisitor<Tail, R, Immutable> {
public:
  typedef R ReturnType;
  typedef Head ParamType;
  //	using StrictVisitor<Head, R>::Visit;
  //	using StrictVisitor<Tail, R>::Visit;
};

template <class Head, typename R>
class StrictVisitor<Typelist<Head, Null_type>, R, Immutable>
    : public StrictVisitor<Head, R, Immutable> {
public:
  typedef R ReturnType;
  typedef Head ParamType;
  using StrictVisitor<Head, R, Immutable>::Visit;
};

////////////////////////////////////////////////////////////////////////////////
// class template NonStrictVisitor
// Implements non-strict visitation (you can implement only part of the Visit
//     functions)
//

template <class R> struct DefaultFunctor {
  template <class T> R operator()(T &) { return R(); }
};

template <class T, typename R = void, Visit_type V = Mutable,
          class F = DefaultFunctor<R>>
class BaseVisitorImpl;

template <class Head, class Tail, typename R, Visit_type V, class F>
class BaseVisitorImpl<Typelist<Head, Tail>, R, V, F>
    : public StrictVisitor<Head, R, V>, public BaseVisitorImpl<Tail, R, V, F> {
public:
  typedef typename StrictVisitor<Head, R, V>::ParamType ParamType;
  virtual R Visit(ParamType & h) { return F()(h); }
};

template <class Head, typename R, Visit_type V, class F>
class BaseVisitorImpl<Typelist<Head, Null_type>, R, V, F>
    : public StrictVisitor<Head, R, V> {
public:
  typedef typename StrictVisitor<Head, R, V>::ParamType ParamType;
  virtual R Visit(ParamType & h) { return F()(h); }
};

/// Visitor
template <class R> struct Strict {};

template <typename R, class TList, Visit_type V = Mutable,
          template <class> class FunctorPolicy = DefaultFunctor>
class Visitor : public BaseVisitorImpl<TList, R, V, FunctorPolicy<R>> {
public:
  typedef R ReturnType;

  template <class Visited> ReturnType GenericVisit(Visited & host) {
    StrictVisitor<Visited, ReturnType, V> & subObj = *this;
    return subObj.Visit(host);
  }
};

template <typename R, class TList, Visit_type V>
class Visitor<R, TList, V, Strict> : public StrictVisitor<TList, R, V> {
public:
  typedef R ReturnType;

  template <class Visited> ReturnType GenericVisit(Visited & host) {
    StrictVisitor<Visited, ReturnType, V> & subObj = *this;
    return subObj.Visit(host);
  }
};

} // namespace akantu

#endif /* __AKANTU_VISITOR_HH__ */
