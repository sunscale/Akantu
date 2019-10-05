/**
 * @file   cohesive_internal_field.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  Internal field for cohesive elements
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "internal_field.hh"

#ifndef __AKANTU_COHESIVE_INTERNAL_FIELD_HH__
#define __AKANTU_COHESIVE_INTERNAL_FIELD_HH__

namespace akantu {

/// internal field class for cohesive materials
template <typename T> class CohesiveInternalField : public InternalField<T> {
public:
  CohesiveInternalField(const ID & id, Material & material);
  ~CohesiveInternalField() override;

  /// initialize the field to a given number of component
  void initialize(UInt nb_component) override;

private:
  CohesiveInternalField operator=(__attribute__((unused))
                                  const CohesiveInternalField & other){};
};

/* -------------------------------------------------------------------------- */
/* Facet Internal Field                                                       */
/* -------------------------------------------------------------------------- */
template <typename T> class FacetInternalField : public InternalField<T> {
public:
  FacetInternalField(const ID & id, Material & material);
  ~FacetInternalField() override;

  /// initialize the field to a given number of component
  void initialize(UInt nb_component) override;
};

} // namespace akantu

#endif /* __AKANTU_COHESIVE_INTERNAL_FIELD_HH__ */
