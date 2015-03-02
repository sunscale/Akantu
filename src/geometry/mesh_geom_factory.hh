/**
 * @file   mesh_geom_factory.hh
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Thu Feb 26 2015
 * @date last modification: Mon Mar 2 2015
 *
 * @brief  Class for constructing the CGAL primitives of a mesh
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef __AKANTU_MESH_GEOM_FACTORY_HH__
#define __AKANTU_MESH_GEOM_FACTORY_HH__

#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_geom_abstract.hh"
#include "tree_type_helper.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

template<UInt d, ElementType el_type>
class MeshGeomFactory : public MeshGeomAbstract {

public:
  /// Construct from mesh
  explicit MeshGeomFactory(const Mesh & mesh);

  /// Desctructor
  virtual ~MeshGeomFactory();

public:
  /// Construct AABB tree for fast intersection computing
  virtual void constructData();

  const typename TreeTypeHelper<d, el_type>::tree & getTree() const { return *data_tree; }

protected:
  typename TreeTypeHelper<d, el_type>::tree * data_tree;
  typename TreeTypeHelper<d, el_type>::container_type primitive_list;
};

__END_AKANTU__

#include "mesh_geom_factory_tmpl.hh"

#endif // __AKANTU_MESH_GEOM_FACTORY_HH__