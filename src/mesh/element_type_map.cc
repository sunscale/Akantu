/**
 * @file   element_type_map.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  A Documented file.
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

/* -------------------------------------------------------------------------- */
#include "fe_engine.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

FEEngineElementTypeMapArrayInitializer::FEEngineElementTypeMapArrayInitializer(
    const FEEngine & fe_engine, UInt nb_component, UInt spatial_dimension,
    GhostType ghost_type, ElementKind element_kind)
    : MeshElementTypeMapArrayInitializer(
          fe_engine.getMesh(), nb_component,
          spatial_dimension == UInt(-2)
              ? fe_engine.getMesh().getSpatialDimension()
              : spatial_dimension,
          ghost_type, element_kind, true, false),
      fe_engine(fe_engine) {}

FEEngineElementTypeMapArrayInitializer::FEEngineElementTypeMapArrayInitializer(
    const FEEngine & fe_engine,
    const ElementTypeMapArrayInitializer::CompFunc & nb_component,
    UInt spatial_dimension, GhostType ghost_type,
    ElementKind element_kind)
    : MeshElementTypeMapArrayInitializer(
          fe_engine.getMesh(), nb_component,
          spatial_dimension == UInt(-2)
              ? fe_engine.getMesh().getSpatialDimension()
              : spatial_dimension,
          ghost_type, element_kind, true, false),
      fe_engine(fe_engine) {}

UInt FEEngineElementTypeMapArrayInitializer::size(
    ElementType type) const {
  return MeshElementTypeMapArrayInitializer::size(type) *
         fe_engine.getNbIntegrationPoints(type, this->ghost_type);
}

FEEngineElementTypeMapArrayInitializer::ElementTypesIteratorHelper
FEEngineElementTypeMapArrayInitializer::elementTypes() const {
  return this->fe_engine.elementTypes(spatial_dimension, ghost_type,
                                      element_kind);
}

} // namespace akantu
