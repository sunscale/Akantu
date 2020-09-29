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

#ifndef AKANTU_MESH_GEOM_FACTORY_HH_
#define AKANTU_MESH_GEOM_FACTORY_HH_

/* -------------------------------------------------------------------------- */
#include <algorithm>
/* -------------------------------------------------------------------------- */
#include "geom_helper_functions.hh"
#include "mesh.hh"
#include "mesh_geom_abstract.hh"
#include "tree_type_helper.hh"
/* -------------------------------------------------------------------------- */

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
  ~MeshGeomFactory() override;

  using TreeTypeHelper = TreeTypeHelper<Primitive, Kernel>;
  using TreeType = typename TreeTypeHelper::tree;
  using ContainerType = typename TreeTypeHelper::container_type;

public:
  /// Construct AABB tree for fast intersection computing
  void constructData(GhostType ghost_type = _not_ghost) override;

  /**
   * @brief Construct a primitive and add it to a list of primitives
   *
   * This function needs to be specialized for every type that is wished to be
   * supported.
   * @param node_coordinates coordinates of the nodes making up the element
   * @param id element number
   * @param list the primitive list (not used inside MeshGeomFactory)
   */
  inline void addPrimitive(const Matrix<Real> & /*node_coordinates*/,
                           UInt /*id*/, ContainerType & /*list*/);

  inline void addPrimitive(const Matrix<Real> & node_coordinates, UInt id);

  /// Getter for the AABB tree
  auto getTree() const -> const TreeType & { return *data_tree; }

  /// Getter for primitive list
  auto getPrimitiveList() const -> const ContainerType & {
    return primitive_list;
  }

protected:
  /// AABB data tree
  TreeType * data_tree{nullptr};

  /// Primitive list
  ContainerType primitive_list;
};

} // namespace akantu

#include "mesh_geom_factory_tmpl.hh"

#endif // AKANTU_MESH_GEOM_FACTORY_HH_
