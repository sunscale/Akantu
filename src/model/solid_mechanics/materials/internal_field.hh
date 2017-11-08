/**
 * @file   internal_field.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Dec 08 2015
 *
 * @brief  Material internal properties
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "aka_common.hh"
#include "element_type_map.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_INTERNAL_FIELD_HH__
#define __AKANTU_INTERNAL_FIELD_HH__

namespace akantu {

class Material;
class FEEngine;

/**
 * class for the internal fields of materials
 * to store values for each quadrature
 */
template <typename T> class InternalField : public ElementTypeMapArray<T> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  InternalField(const ID & id, Material & material);
  ~InternalField() override;

  /// This constructor is only here to let cohesive elements compile
  InternalField(const ID & id, Material & material, FEEngine & fem,
                const ElementTypeMapArray<UInt> & element_filter);

  /// More general constructor
  InternalField(const ID & id, Material & material, UInt dim, FEEngine & fem,
                const ElementTypeMapArray<UInt> & element_filter);

  InternalField(const ID & id, const InternalField<T> & other);

private:
  InternalField operator=(__attribute__((unused))
                          const InternalField & other){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// function to reset the FEEngine for the internal field
  virtual void setFEEngine(FEEngine & fe_engine);

  /// function to reset the element kind for the internal
  virtual void setElementKind(ElementKind element_kind);

  /// initialize the field to a given number of component
  virtual void initialize(UInt nb_component);

  /// activate the history of this field
  virtual void initializeHistory();

  /// resize the arrays and set the new element to 0
  virtual void resize();

  /// set the field to a given value v
  virtual void setDefaultValue(const T & v);

  /// reset all the fields to the default value
  virtual void reset();

  /// save the current values in the history
  virtual void saveCurrentValues();

  /// remove the quadrature points corresponding to suppressed elements
  virtual void
  removeIntegrationPoints(const ElementTypeMapArray<UInt> & new_numbering);

  /// print the content
  void printself(std::ostream & stream, int indent = 0) const override;

  /// get the default value
  inline operator T() const;

  virtual FEEngine & getFEEngine() { return *fem; }

  virtual const FEEngine & getFEEngine() const { return *fem; }

  /// AKANTU_GET_MACRO(FEEngine, *fem, FEEngine &);

protected:
  /// initialize the arrays in the ElementTypeMapArray<T>
  void internalInitialize(UInt nb_component);

  /// set the values for new internals
  virtual void setArrayValues(T * begin, T * end);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  typedef typename ElementTypeMapArray<T>::type_iterator type_iterator;
  typedef typename ElementTypeMapArray<UInt>::type_iterator filter_type_iterator;

  /// get the type iterator on all types contained in the internal field
  type_iterator firstType(const GhostType & ghost_type = _not_ghost) const {
    return ElementTypeMapArray<T>::firstType(this->spatial_dimension, ghost_type,
                                             this->element_kind);
  }

  /// get the type iterator on the last type contained in the internal field
  type_iterator lastType(const GhostType & ghost_type = _not_ghost) const {
    return ElementTypeMapArray<T>::lastType(this->spatial_dimension, ghost_type,
                                            this->element_kind);
  }

    /// get the type iterator on all types contained in the internal field
  filter_type_iterator filterFirstType(const GhostType & ghost_type = _not_ghost) const {
    return this->element_filter.firstType(this->spatial_dimension, ghost_type,
                           this->element_kind);
  }

  /// get the type iterator on the last type contained in the internal field
  filter_type_iterator filterLastType(const GhostType & ghost_type = _not_ghost) const {
    return this->element_filter.lastType(this->spatial_dimension, ghost_type,
                          this->element_kind);
  }

  /// get the array for a given type of the element_filter
  const Array<UInt> getFilter(const ElementType & type,
                              const GhostType & ghost_type = _not_ghost) const {
    return this->element_filter(type, ghost_type);
  }

  /// get the Array corresponding to the type en ghost_type specified
  virtual Array<T> & operator()(const ElementType & type,
                                const GhostType & ghost_type = _not_ghost) {
    return ElementTypeMapArray<T>::operator()(type, ghost_type);
  }

  virtual const Array<T> &
  operator()(const ElementType & type,
             const GhostType & ghost_type = _not_ghost) const {
    return ElementTypeMapArray<T>::operator()(type, ghost_type);
  }

  virtual Array<T> & previous(const ElementType & type,
                              const GhostType & ghost_type = _not_ghost) {
    AKANTU_DEBUG_ASSERT(previous_values != nullptr,
                        "The history of the internal "
                            << this->getID() << " has not been activated");
    return this->previous_values->operator()(type, ghost_type);
  }

  virtual const Array<T> &
  previous(const ElementType & type,
           const GhostType & ghost_type = _not_ghost) const {
    AKANTU_DEBUG_ASSERT(previous_values != nullptr,
                        "The history of the internal "
                            << this->getID() << " has not been activated");
    return this->previous_values->operator()(type, ghost_type);
  }

  virtual InternalField<T> & previous() {
    AKANTU_DEBUG_ASSERT(previous_values != nullptr,
                        "The history of the internal "
                            << this->getID() << " has not been activated");
    return *(this->previous_values);
  }

  virtual const InternalField<T> & previous() const {
    AKANTU_DEBUG_ASSERT(previous_values != nullptr,
                        "The history of the internal "
                            << this->getID() << " has not been activated");
    return *(this->previous_values);
  }

  /// check if the history is used or not
  bool hasHistory() const { return (previous_values != nullptr); }

  /// get the kind treated by the internal
  const ElementKind & getElementKind() const { return element_kind; }

  /// return the number of components
  UInt getNbComponent() const { return nb_component; }

  /// return the spatial dimension corresponding to the internal element type
  /// loop filter
  UInt getSpatialDimension() const { return this->spatial_dimension; }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// the material for which this is an internal parameter
  Material & material;

  /// the fem containing the mesh and the element informations
  FEEngine * fem;

  /// Element filter if needed
  const ElementTypeMapArray<UInt> & element_filter;

  /// default value
  T default_value;

  /// spatial dimension of the element to consider
  UInt spatial_dimension;

  /// ElementKind of the element to consider
  ElementKind element_kind;

  /// Number of component of the internal field
  UInt nb_component;

  /// Is the field initialized
  bool is_init;

  /// previous values
  InternalField<T> * previous_values;
};

/// standard output stream operator
template <typename T>
inline std::ostream & operator<<(std::ostream & stream,
                                 const InternalField<T> & _this) {
  _this.printself(stream);
  return stream;
}

} // akantu

#endif /* __AKANTU_INTERNAL_FIELD_HH__ */
