/**
 * @file   element.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Sep 02 2014
 * @date last modification: Sat Jul 11 2015
 *
 * @brief  Element helper class
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

#ifndef __AKANTU_ELEMENT_HH__
#define __AKANTU_ELEMENT_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
/* Element                                                                    */
/* -------------------------------------------------------------------------- */
class Element;
extern const Element ElementNull;

class Element {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

public:
  explicit Element(ElementType type = _not_defined, UInt element = 0,
          GhostType ghost_type = _not_ghost, ElementKind kind = _ek_regular) :
    type(type), element(element),
    ghost_type(ghost_type), kind(kind) {};

  Element(const Element & element) {
    this->type    = element.type;
    this->element = element.element;
    this->ghost_type = element.ghost_type;
    this->kind = element.kind;
  }

  virtual ~Element() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

public:

  inline bool operator==(const Element & elem) const {
    return ((element == elem.element)
            && (type == elem.type)
            && (ghost_type == elem.ghost_type)
            && (kind == elem.kind));
  }

  inline bool operator!=(const Element & elem) const {
    return ((element != elem.element)
            || (type != elem.type)
            || (ghost_type != elem.ghost_type)
            || (kind != elem.kind));
  }

  bool operator<(const Element& rhs) const {
    bool res = (rhs == ElementNull) || ((this->kind < rhs.kind) ||
                                        ((this->kind == rhs.kind) &&
                                         ((this->ghost_type < rhs.ghost_type) ||
                                          ((this->ghost_type == rhs.ghost_type) &&
                                           ((this->type < rhs.type) ||
                                            ((this->type == rhs.type) &&
                                             (this->element < rhs.element)))))));
    return res;
  }

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
  
public:

  const ElementType & getType(){return type;}
  const UInt & getIndex(){return element;};
  const GhostType & getGhostType(){return ghost_type;}
  const ElementKind & getElementKind(){return kind;}

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

public:
  ElementType type;
  UInt element;
  GhostType ghost_type;
  ElementKind kind;
};

struct CompElementLess {
  bool operator() (const Element& lhs, const Element& rhs) const {
    return lhs < rhs;
  }
};

} // akantu

#endif /* __AKANTU_ELEMENT_HH__ */
