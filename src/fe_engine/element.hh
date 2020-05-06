/**
 * @file   element.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Sep 02 2014
 * @date last modification: Tue Jan 23 2018
 *
 * @brief  Element helper class
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
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_ELEMENT_HH__
#define __AKANTU_ELEMENT_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
/* Element                                                                    */
/* -------------------------------------------------------------------------- */
class Element {
public:
  ElementType type;
  UInt element;
  GhostType ghost_type;

  // ElementKind kind;
  // ElementType type{_not_defined};
  // UInt element{0};
  // GhostType ghost_type{_not_ghost};
  // ElementKind kind{_ek_regular};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  inline ElementKind kind() const;

  inline bool operator==(const Element & elem) const {
    return std::tie(type, element, ghost_type) ==
           std::tie(elem.type, elem.element, elem.ghost_type);
  }

  inline bool operator!=(const Element & elem) const {
    return std::tie(type, element, ghost_type) !=
           std::tie(elem.type, elem.element, elem.ghost_type);
  }

  // inline bool operator==(const Element & elem) const {
  //   return ((element == elem.element) && (type == elem.type) &&
  //           (ghost_type == elem.ghost_type) && (kind == elem.kind));
  // }

  // inline bool operator!=(const Element & elem) const {
  //   return ((element != elem.element) || (type != elem.type) ||
  //           (ghost_type != elem.ghost_type) || (kind != elem.kind));
  // }

  inline bool operator<(const Element & rhs) const;
};

namespace {
  const Element ElementNull{_not_defined, UInt(-1), _casper};
  //      Element{_not_defined, 0, _casper, _ek_not_defined};
} // namespace

/* -------------------------------------------------------------------------- */
inline bool Element::operator<(const Element & rhs) const {
  // bool res =
  //     (rhs == ElementNull) ||
  //     ((this->kind < rhs.kind) ||
  //      ((this->kind == rhs.kind) &&
  //       ((this->ghost_type < rhs.ghost_type) ||
  //        ((this->ghost_type == rhs.ghost_type) &&
  //         ((this->type < rhs.type) ||
  //          ((this->type == rhs.type) && (this->element < rhs.element)))))));
  return ((rhs == ElementNull) ||
          std::tie(ghost_type, type, element) <
              std::tie(rhs.ghost_type, rhs.type, rhs.element));
}

} // namespace akantu

namespace std {
inline string to_string(const akantu::Element & _this) {
  if (_this == akantu::ElementNull) {
    return "ElementNull";
  }

  string str = "Element [" + to_string(_this.type) + ", " +
               to_string(_this.element) + ", " + to_string(_this.ghost_type) +
               "]";
  return str;
}
} // namespace std

namespace akantu {

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream, const Element & _this) {
  stream << std::to_string(_this);
  return stream;
}
} // namespace akantu

#endif /* __AKANTU_ELEMENT_HH__ */
