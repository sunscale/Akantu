/**
 * @file   fe_engine.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Jul 20 2010
 * @date last modification: Thu Feb 01 2018
 *
 * @brief  Implementation of the FEEngine class
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
#include "aka_memory.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
FEEngine::FEEngine(Mesh & mesh, UInt element_dimension, const ID & id,
                   MemoryID memory_id)
    : Memory(id, memory_id), mesh(mesh),
      normals_on_integration_points("normals_on_quad_points", id, memory_id) {
  AKANTU_DEBUG_IN();
  this->element_dimension = (element_dimension != _all_dimensions)
                                ? element_dimension
                                : mesh.getSpatialDimension();

  this->mesh.registerEventHandler(*this, _ehp_fe_engine);

  init();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FEEngine::init() {}

/* -------------------------------------------------------------------------- */
FEEngine::~FEEngine() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
typename FEEngine::ElementTypesIteratorHelper
FEEngine::elementTypes(UInt dim, GhostType ghost_type, ElementKind kind) const {
  return this->getIntegratorInterface().getJacobians().elementTypes(
      dim, ghost_type, kind);
}

/* -------------------------------------------------------------------------- */
void FEEngine::printself(std::ostream & stream, int indent) const {
  std::string space(indent, AKANTU_INDENT);

  stream << space << "FEEngine [" << std::endl;
  stream << space << " + id                : " << id << std::endl;
  stream << space << " + element dimension : " << element_dimension
         << std::endl;

  stream << space << " + mesh [" << std::endl;
  mesh.printself(stream, indent + 2);
  stream << space << AKANTU_INDENT << "]" << std::endl;

  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
