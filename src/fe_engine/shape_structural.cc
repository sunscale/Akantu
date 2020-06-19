/**
 * @file   shape_structural.cc
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Feb 21 2013
 * @date last modification: Thu Dec 21 2017
 *
 * @brief  ShapeStructural implementation
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
#include "shape_structural.hh"
#include "aka_memory.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <>
ShapeStructural<_ek_structural>::ShapeStructural(Mesh & mesh,
                                                 UInt spatial_dimension,
                                                 const ID & id,
                                                 const MemoryID & memory_id)
    : ShapeFunctions(mesh, spatial_dimension, id, memory_id),
      rotation_matrices("rotation_matrices", id, memory_id) {}

/* -------------------------------------------------------------------------- */
template <> ShapeStructural<_ek_structural>::~ShapeStructural() = default;

/* -------------------------------------------------------------------------- */

} // namespace akantu
