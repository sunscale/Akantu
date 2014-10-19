/**
 * @file   cohesive_internal_field.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Nov 13 2013
 * @date last modification: Tue Jul 29 2014
 *
 * @brief  Internal field for cohesive elements
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
#include "internal_field.hh"

#ifndef __AKANTU_COHESIVE_INTERNAL_FIELD_HH__
#define __AKANTU_COHESIVE_INTERNAL_FIELD_HH__

__BEGIN_AKANTU__

template<typename T>
class CohesiveInternalField : public InternalField<T> {
public:
  CohesiveInternalField(const ID & id, Material & material);
  virtual ~CohesiveInternalField();
  void initialize(UInt nb_component);
private:
  CohesiveInternalField operator=(__attribute__((unused)) const CohesiveInternalField & other) {};

};


/* -------------------------------------------------------------------------- */
/* Facet Internal Field                                                       */
/* -------------------------------------------------------------------------- */
template<typename T>
class FacetInternalField : public InternalField<T> {
public:
  FacetInternalField(const ID & id, Material & material);
  virtual ~FacetInternalField();
  void initialize(UInt nb_component);
};

__END_AKANTU__

#endif /* __AKANTU_COHESIVE_INTERNAL_FIELD_HH__ */
