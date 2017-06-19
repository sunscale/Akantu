/**
 * @file   shape_linked.cc
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jul 15 2011
 * @date last modification: Sun Oct 19 2014
 *
 * @brief  ShapeLinked implementation
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
#include "aka_memory.hh"
#include "mesh.hh"
#include "shape_linked.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

#if defined(AKANTU_STRUCTURAL_MECHANICS)
/* -------------------------------------------------------------------------- */
template <>
ShapeLinked<_ek_structural>::ShapeLinked(Mesh & mesh, const ID & id, const MemoryID & memory_id) :
  ShapeFunctions(mesh, id, memory_id)
{

}


/* -------------------------------------------------------------------------- */
template <>
ShapeLinked<_ek_structural>::~ShapeLinked() {
  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {
    GhostType ghost_type = *gt;

    // delete all the shapes id
    ElementTypeMapMultiReal::type_iterator s_type_it  =
      shapes.firstType(_all_dimensions, ghost_type, _ek_structural);
    ElementTypeMapMultiReal::type_iterator s_type_end =
      shapes.lastType (_all_dimensions, ghost_type, _ek_structural);

    for(; s_type_it != s_type_end; ++s_type_it) {
      delete [] shapes(*s_type_it, ghost_type);
    }


    // delete all the shapes derivatives id
    ElementTypeMapMultiReal::type_iterator sd_type_it  =
      shapes_derivatives.firstType(_all_dimensions, ghost_type, _ek_structural);
    ElementTypeMapMultiReal::type_iterator sd_type_end =
      shapes_derivatives.lastType (_all_dimensions, ghost_type, _ek_structural);

    for(; sd_type_it != sd_type_end; ++sd_type_it) {
      delete [] shapes_derivatives(*sd_type_it, ghost_type);
    }
  }
}
#endif

} // akantu
