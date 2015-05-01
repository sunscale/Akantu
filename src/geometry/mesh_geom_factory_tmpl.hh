/**
 * @file   mesh_geom_factory_tmpl.hh
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Thu Feb 26 2015
 * @date last modification: Fri Mar 6 2015
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

typedef CGAL::Cartesian<Real> Cartesian;

template<UInt dim, ElementType type, class Primitive, class Kernel>
MeshGeomFactory<dim, type, Primitive, Kernel>::MeshGeomFactory(const Mesh & mesh) :
  MeshGeomAbstract(mesh),
  data_tree(NULL),
  primitive_list()
{}

template<UInt dim, ElementType type, class Primitive, class Kernel>
MeshGeomFactory<dim, type, Primitive, Kernel>::~MeshGeomFactory() {
  delete data_tree;
}

/**
 * This function loops over the elements of `type` in the mesh and creates the
 * AABB tree of geometrical primitves (`data_tree`).
 */
template<UInt dim, ElementType type, class Primitive, class Kernel>
void MeshGeomFactory<dim, type, Primitive, Kernel>::constructData() {
  AKANTU_DEBUG_IN();

  const GhostType ghost_type = _not_ghost;

  primitive_list.clear();

  UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);
  const Array<UInt> & connectivity = mesh.getConnectivity(type, ghost_type);
  const Array<Real> & nodes = mesh.getNodes();

  Array<UInt>::const_vector_iterator begin = connectivity.begin(nb_nodes_per_element);
  Array<UInt>::const_vector_iterator it    = connectivity.begin(nb_nodes_per_element);
  Array<UInt>::const_vector_iterator end   = connectivity.end(nb_nodes_per_element);

  // This loop builds the list of primitives
  for (; it != end ; ++it) {
    const Vector<UInt> & el_connectivity = *it;
    Matrix<Real> node_coordinates(dim, nb_nodes_per_element);

    for (UInt i = 0 ; i < nb_nodes_per_element ; i++)
      for (UInt j = 0 ; j < dim ; j++)
        node_coordinates(j, i) = nodes(el_connectivity(i), j);

    addPrimitive(node_coordinates, it - begin);
  }

  delete data_tree;

  data_tree = new typename TreeTypeHelper<Primitive, Kernel>::tree(primitive_list.begin(), primitive_list.end());

  AKANTU_DEBUG_OUT();
}

