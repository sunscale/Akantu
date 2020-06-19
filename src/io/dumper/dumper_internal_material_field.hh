/**
 * @file   dumper_internal_material_field.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  description of material internal field
 *
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

#ifndef __AKANTU_DUMPER_INTERNAL_MATERIAL_FIELD_HH__
#define __AKANTU_DUMPER_INTERNAL_MATERIAL_FIELD_HH__
/* -------------------------------------------------------------------------- */
#include "dumper_quadrature_point_iterator.hh"
#ifdef AKANTU_IGFEM
#include "dumper_igfem_material_internal_field.hh"
#endif
/* -------------------------------------------------------------------------- */
namespace akantu {
namespace dumpers {
/* -------------------------------------------------------------------------- */

template <typename T, bool filtered = false>
class InternalMaterialField
    : public GenericElementalField<SingleType<T, Vector, filtered>,
                                   quadrature_point_iterator> {

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */

public:
  using types = SingleType<T, Vector, filtered>;
  using parent = GenericElementalField<types, quadrature_point_iterator>;
  using field_type = typename types::field_type;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

  InternalMaterialField(const field_type & field,
                        UInt spatial_dimension = _all_dimensions,
                        GhostType ghost_type = _not_ghost,
                        ElementKind element_kind = _ek_not_defined)
      : parent(field, spatial_dimension, ghost_type, element_kind) {}
};

} // namespace dumpers
} // namespace akantu

#endif /* __AKANTU_DUMPER_INTERNAL_MATERIAL_FIELD_HH__ */
