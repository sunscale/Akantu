/**
 * @file   mesh_geom_factory.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Feb 27 2015
 * @date last modification: Mon Jun 19 2017
 *
 * @brief  Class for constructing the CGAL primitives of a mesh
 *
 * @section LICENSE
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef __AKANTU_MESH_GEOM_FACTORY_HH__
#define __AKANTU_MESH_GEOM_FACTORY_HH__

#include "aka_common.hh"
#include "geom_helper_functions.hh"
#include "mesh.hh"
#include "mesh_geom_abstract.hh"
#include "tree_type_helper.hh"

#include <algorithm>

/* -------------------------------------------------------------------------- */

namespace akantu {

/**
 * @brief Class used to construct AABB tree for intersection computations
 *
 * This class constructs a CGAL AABB tree of one type of element in a mesh
 * for fast intersection computations.
 */
template <UInt dim, ElementType el_type, class Primitive, class Kernel>
class MeshGeomFactory : public MeshGeomAbstract {

public:
  /// Construct from mesh
  explicit MeshGeomFactory(Mesh & mesh);

  /// Desctructor
  virtual ~MeshGeomFactory();

public:
  /// Construct AABB tree for fast intersection computing
  virtual void constructData(GhostType ghost_type = _not_ghost);

  /**
   * @brief Construct a primitive and add it to a list of primitives
   *
   * This function needs to be specialized for every type that is wished to be
   * supported.
   * @param node_coordinates coordinates of the nodes making up the element
   * @param id element number
   * @param list the primitive list (not used inside MeshGeomFactory)
   */
  inline void addPrimitive(
      const Matrix<Real> & node_coordinates, UInt id,
      typename TreeTypeHelper<Primitive, Kernel>::container_type & list);

  inline void addPrimitive(const Matrix<Real> & node_coordinates, UInt id);

  /// Getter for the AABB tree
  const typename TreeTypeHelper<Primitive, Kernel>::tree & getTree() const {
    return *data_tree;
  }

  /// Getter for primitive list
  const typename TreeTypeHelper<Primitive, Kernel>::container_type &
  getPrimitiveList() const {
    return primitive_list;
  }

protected:
  /// AABB data tree
  typename TreeTypeHelper<Primitive, Kernel>::tree * data_tree;

  /// Primitive list
  typename TreeTypeHelper<Primitive, Kernel>::container_type primitive_list;
};

} // namespace akantu

#include "mesh_geom_factory_tmpl.hh"

#endif // __AKANTU_MESH_GEOM_FACTORY_HH__
