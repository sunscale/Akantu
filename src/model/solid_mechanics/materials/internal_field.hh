/**
 * @file   internal_field.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Nov 13 2013
 * @date last modification: Tue Sep 02 2014
 *
 * @brief  Material internal properties
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
#include "aka_common.hh"
#include "element_type_map.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_INTERNAL_FIELD_HH__
#define __AKANTU_INTERNAL_FIELD_HH__

__BEGIN_AKANTU__

class Material;
class FEEngine;

template<typename T>
class InternalField : public ElementTypeMapArray<T> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  InternalField(const ID & id, Material & material);
  virtual ~InternalField();

protected:
  InternalField(const ID & id, Material & material, FEEngine & fem,
		const ElementTypeMapArray<UInt> & element_filter);

  InternalField(const ID & id, const InternalField<T> & other);

private:
  InternalField operator=(__attribute__((unused)) const InternalField & other) {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
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
  virtual void removeQuadraturePoints(const ElementTypeMapArray<UInt> & new_numbering);

  /// print the content
  virtual void printself(std::ostream & stream, UInt indent = 0) const;

protected:
  /// initialize the arrays in the ElementTypeMapArray<T>
  void internalInitialize(UInt nb_component);

  /// set the values for new internals
  virtual void setArrayValues(T * begin, T * end);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /// get the Array corresponding to the type en ghost_type specified
  virtual Array<T> & operator()(const ElementType & type, const GhostType & ghost_type = _not_ghost) {
    return ElementTypeMapArray<T>::operator()(type, ghost_type);
  }

  virtual const Array<T> & operator()(const ElementType & type, const GhostType & ghost_type = _not_ghost) const {
    return ElementTypeMapArray<T>::operator()(type, ghost_type);
  }

  virtual Array<T> & previous(const ElementType & type, const GhostType & ghost_type = _not_ghost) {
    AKANTU_DEBUG_ASSERT(previous_values != NULL,
			"The history of the internal " << this->getID()
			<< " has not been activated");
    return this->previous_values->operator()(type, ghost_type);
  }

  virtual const Array<T> & previous(const ElementType & type, const GhostType & ghost_type = _not_ghost) const {
    AKANTU_DEBUG_ASSERT(previous_values != NULL,
			"The history of the internal " << this->getID()
			<< " has not been activated");
    return this->previous_values->operator()(type, ghost_type);
  }


  virtual InternalField<T> & previous() {
    AKANTU_DEBUG_ASSERT(previous_values != NULL,
			"The history of the internal " << this->getID()
			<< " has not been activated");
    return *(this->previous_values);
  }

  virtual const InternalField<T> & previous() const {
    AKANTU_DEBUG_ASSERT(previous_values != NULL,
			"The history of the internal " << this->getID()
			<< " has not been activated");
    return *(this->previous_values);
  }

  /// check if the history is used or not
  bool hasHistory() const { return (previous_values != NULL); }

  /// get the kind treated by the internal
  const ElementKind & getElementKind() const {return element_kind;};


  /// return the number of components
  UInt getNbComponent(){return nb_component;}
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// the material for which this is an internal parameter
  Material & material;

  /// the fem containing the mesh and the element informations
  FEEngine & fem;

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
template<typename T>
inline std::ostream & operator <<(std::ostream & stream, const InternalField<T> & _this)
{
  _this.printself(stream);
  return stream;
}

__END_AKANTU__

#endif /* __AKANTU_INTERNAL_FIELD_HH__ */
