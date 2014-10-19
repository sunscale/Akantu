/**
 * @file   contact_common.hh
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date creation: Fri Jan 04 2013
 * @date last modification: Wed Sep 17 2014
 *
 * @brief  Forward declarations for contact classes
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

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_CONTACT_COMMON_HH__
#define __AKANTU_CONTACT_COMMON_HH__

#include <map>

#include "aka_visitor.hh"
#include <array/expr.hpp>

#include "cppargparse.hh"

//#define DEBUG_CONTACT 1

__BEGIN_AKANTU__

using std::cout;
using std::endl;

typedef array::Array<1, Real> vector_type;
typedef array::Array<2, Real> matrix_type;

using array::transpose;

//! Enumerated type used for the Contact overloaded operator[] that returns real
// values
enum Contact_parameter_type {
  Epsilon,
  Alpha,
  Multiplier_tol,
  Newton_tol,
  Multiplier_max_steps,
  Newton_max_steps
};

//! Enumerated type used for the Contact overloaded operator[] that returns
// boolean values
enum Contact_flag_type {
  Verbose,
  Automatic_penalty_parameter
};

//! Discretization types
enum Discretization_type {
  Node_to_node,
  Node_to_segment,
  Segment_to_segment
};

//! Contact type
enum Contact_type {
  Self_contact,
  No_self_contact
};

struct EmptyType {};
class NullType {};

//! This functor is called when the visitor is not implemented for a particular
// object.
template <class R> struct Discretization_visitor_default {
  template <class T> R operator()(T &t) {
    cout << "*** WARNING *** No implementation for discretization visitor"
         << endl;
  }
};

template <Discretization_type d> class Contact_discretization;

typedef Contact_discretization<Node_to_node> N2N_c;
typedef Contact_discretization<Node_to_segment> N2S_c;
typedef Contact_discretization<Segment_to_segment> S2S_c;

template <AnalysisMethod s, ContactResolutionMethod r,
          ContactImplementationMethod i = _none>
struct SelectResolution {
  constexpr static const AnalysisMethod analysis = s;
  constexpr static const ContactResolutionMethod method = r;
  constexpr static const ContactImplementationMethod implementation = i;
};

template <int Dim, AnalysisMethod s, ContactResolutionMethod r>
class ContactResolution;

template <int Dim, template <int> class Search_policy, class Resolution_policy>
class Contact;

// parsers
extern cppargparse::ArgumentParser contact_argparser;
extern Parser contact_parser;

class ContactParameters {
protected:
  typedef std::map<Contact_parameter_type, Real> options_map;
  typedef std::map<Contact_flag_type, bool> flag_map;

  options_map options_;
  flag_map flags_;

public:
  //! Overloaded operator[] (const) that returns real values
  Real operator[](Contact_parameter_type p) const {
    auto it = options_.find(p);
    assert(it != options_.end());
    return it->second;
  }

  //! Overloaded operator[] (non-const) that returns real values
  Real &operator[](Contact_parameter_type p) { return options_[p]; }

  //! Overloaded operator[] (const) that returns flags
  bool operator[](Contact_flag_type f) const {
    auto it = flags_.find(f);
    return it->second;
  }

  //! Overloaded operator[] (non-const) that returns bool values
  bool &operator[](Contact_flag_type p) { return flags_[p]; }

  virtual void initialize() {}
};

template <typename T> T Heaviside(T v) { return v < 0 ? 0 : 1; }

template <typename T> T Macauley(T v) { return v < 0 ? 0 : v; }

template <typename value_type = Real> array::Array<2, value_type> eye(UInt d) {
  array::Array<2, value_type> I(d);

  for (UInt i = 0; i < d; ++i)
    I(i, i) = 1.;
  return I;
}

__END_AKANTU__

#endif /* __AKANTU_CONTACT_COMMON_HH__ */
