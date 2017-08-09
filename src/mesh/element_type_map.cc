/**
 * @file   element_type_map.cc
 *
 * @author Nicolas Richart
 *
 * @date creation  Tue Jun 27 2017
 *
 * @brief A Documented file.
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "fe_engine.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

FEEngineElementTypeMapArrayInializer::FEEngineElementTypeMapArrayInializer(
    const FEEngine & fe_engine, UInt nb_component, UInt spatial_dimension,
    const GhostType & ghost_type, const ElementKind & element_kind)
    : MeshElementTypeMapArrayInializer(
          fe_engine.getMesh(),
          nb_component,
          spatial_dimension == UInt(-2)
              ? fe_engine.getMesh().getSpatialDimension()
              : spatial_dimension,
          ghost_type, element_kind, true, false),
      fe_engine(fe_engine) {}

UInt FEEngineElementTypeMapArrayInializer::size(
    const ElementType & type) const {
  return MeshElementTypeMapArrayInializer::size(type) *
         fe_engine.getNbIntegrationPoints(type, this->ghost_type);
}

FEEngineElementTypeMapArrayInializer::ElementTypesIteratorHelper
FEEngineElementTypeMapArrayInializer::elementTypes() const {
  return this->fe_engine.elementTypes(spatial_dimension, ghost_type,
                                      element_kind);
}
} // namespace akantu
