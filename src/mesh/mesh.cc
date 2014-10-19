/**
 * @file   mesh.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Fri Sep 05 2014
 *
 * @brief  class handling meshes
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include <sstream>

#include "aka_config.hh"

/* -------------------------------------------------------------------------- */
#include "mesh.hh"
#include "group_manager_inline_impl.cc"
#include "mesh_io.hh"
#include "element_class.hh"
#include "static_communicator.hh"
#include "element_group.hh"
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "dumper_field.hh"
#  include "dumper_material_internal_field.hh"
#endif
/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

const Element ElementNull(_not_defined, 0);

/* -------------------------------------------------------------------------- */
void Element::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);
  stream << space << "Element [" << type << ", " << element << ", " << ghost_type << "]";
}


/* -------------------------------------------------------------------------- */
Mesh::Mesh(UInt spatial_dimension,
	   const ID & id,
	   const MemoryID & memory_id) :
  Memory(id, memory_id),
  GroupManager(*this, id + ":group_manager", memory_id),
  nodes_global_ids(NULL), nodes_type(0, 1, id + ":nodes_type"),
  created_nodes(true),
  connectivities("connectivities", id),
  normals("normals", id),
  spatial_dimension(spatial_dimension),
  types_offsets(Array<UInt>((UInt) _max_element_type + 1, 1)),
  ghost_types_offsets(Array<UInt>((UInt) _max_element_type + 1, 1)),
  lower_bounds(spatial_dimension,0.),
  upper_bounds(spatial_dimension,0.),
  size(spatial_dimension, 0.),
  local_lower_bounds(spatial_dimension,0.),
  local_upper_bounds(spatial_dimension,0.),
  mesh_data("mesh_data", id, memory_id),
  mesh_facets(NULL) {
  AKANTU_DEBUG_IN();

  this->nodes = &(alloc<Real>(id + ":coordinates", 0, this->spatial_dimension));

  nb_global_nodes = 0;

  init();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Mesh::Mesh(UInt spatial_dimension,
           const ID & nodes_id,
           const ID & id,
           const MemoryID & memory_id) :
  Memory(id, memory_id),
  GroupManager(*this, id + ":group_manager", memory_id),
  nodes_global_ids(NULL), nodes_type(0, 1, id + ":nodes_type"),
  created_nodes(false),
  connectivities("connectivities", id),
  normals("normals", id),
  spatial_dimension(spatial_dimension),
  types_offsets(Array<UInt>((UInt) _max_element_type + 1, 1)),
  ghost_types_offsets(Array<UInt>((UInt) _max_element_type + 1, 1)),
  lower_bounds(spatial_dimension,0.),
  upper_bounds(spatial_dimension,0.),
  size(spatial_dimension, 0.),
  local_lower_bounds(spatial_dimension,0.),
  local_upper_bounds(spatial_dimension,0.),
  mesh_data("mesh_data", id, memory_id),
  mesh_facets(NULL) {
  AKANTU_DEBUG_IN();

  this->nodes = &(getArray<Real>(nodes_id));
  nb_global_nodes = nodes->getSize();

  init();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Mesh::Mesh(UInt spatial_dimension,
           Array<Real> & nodes,
           const ID & id,
           const MemoryID & memory_id) :
  Memory(id, memory_id),
  GroupManager(*this, id + ":group_manager", memory_id),
  nodes_global_ids(NULL), nodes_type(0, 1, id + ":nodes_type"),
  created_nodes(false),
  connectivities("connectivities", id),
  normals("normals", id),
  spatial_dimension(spatial_dimension),
  types_offsets(Array<UInt>(_max_element_type + 1, 1)),
  ghost_types_offsets(Array<UInt>(_max_element_type + 1, 1)),
  lower_bounds(spatial_dimension,0.),
  upper_bounds(spatial_dimension,0.),
  size(spatial_dimension, 0.),
  local_lower_bounds(spatial_dimension,0.),
  local_upper_bounds(spatial_dimension,0.),
  mesh_data("mesh_data", id, memory_id),
  mesh_facets(NULL) {
  AKANTU_DEBUG_IN();

  this->nodes = &(nodes);
  nb_global_nodes = nodes.getSize();

  init();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Mesh & Mesh::initMeshFacets(const ID & id) {
  AKANTU_DEBUG_IN();

  if (!mesh_facets) {
    mesh_facets = new Mesh(spatial_dimension,
                           *(this->nodes),
			   getID()+":"+id,
                           getMemoryID());

    mesh_facets->mesh_parent = this;
    mesh_facets->is_mesh_facets = true;
  }

  AKANTU_DEBUG_OUT();
  return *mesh_facets;
}

/* -------------------------------------------------------------------------- */
void Mesh::defineMeshParent(const Mesh & mesh) {
  AKANTU_DEBUG_IN();

  this->mesh_parent = &mesh;
  this->is_mesh_facets = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Mesh::init() {
  this->is_mesh_facets = false;
  this->mesh_parent = NULL;
  this->is_distributed = false;
  //  computeBoundingBox();
}

/* -------------------------------------------------------------------------- */
Mesh::~Mesh() {
  AKANTU_DEBUG_IN();

  delete mesh_facets;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Mesh::read (const std::string & filename, const MeshIOType & mesh_io_type) {
  MeshIO mesh_io;
  mesh_io.read(filename, *this, mesh_io_type);

  type_iterator it   = this->firstType(spatial_dimension, _not_ghost, _ek_not_defined);
  type_iterator last = this->lastType(spatial_dimension, _not_ghost, _ek_not_defined);
  if(it == last) AKANTU_EXCEPTION("The mesh contained in the file " << filename
				  << " does not seam to be of the good dimension."
				  << " No element of dimension " << spatial_dimension
				  << " where read.");
}

/* -------------------------------------------------------------------------- */
void Mesh::write(const std::string & filename, const MeshIOType & mesh_io_type) {
  MeshIO mesh_io;
  mesh_io.write(filename, *this, mesh_io_type);
}

/* -------------------------------------------------------------------------- */
void Mesh::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "Mesh [" << std::endl;
  stream << space << " + id                : " << getID() << std::endl;
  stream << space << " + spatial dimension : " << this->spatial_dimension << std::endl;
  stream << space << " + nodes [" << std::endl;
  nodes->printself(stream, indent+2);
  stream << space << " + connectivities [" << std::endl;
  connectivities.printself(stream, indent+2);
  stream << space << " ]" << std::endl;

  GroupManager::printself(stream, indent + 1);
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
void Mesh::computeBoundingBox(){
  AKANTU_DEBUG_IN();
  for (UInt k = 0; k < spatial_dimension; ++k) {
    local_lower_bounds(k) =   std::numeric_limits<double>::max();
    local_upper_bounds(k) = - std::numeric_limits<double>::max();
  }

  for (UInt i = 0; i < nodes->getSize(); ++i) {
    //    if(!isPureGhostNode(i))
    for (UInt k = 0; k < spatial_dimension; ++k) {
      local_lower_bounds(k) = std::min(local_lower_bounds[k], (*nodes)(i, k));
      local_upper_bounds(k) = std::max(local_upper_bounds[k], (*nodes)(i, k));
    }
  }

  if (this->is_distributed) {
    StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();

    Real reduce_bounds[2 * spatial_dimension];
    for (UInt k = 0; k < spatial_dimension; ++k) {
      reduce_bounds[2*k    ] =   local_lower_bounds(k);
      reduce_bounds[2*k + 1] = - local_upper_bounds(k);
    }
    
    comm.allReduce(reduce_bounds, 2 * spatial_dimension, _so_min);

    for (UInt k = 0; k < spatial_dimension; ++k) {
      lower_bounds(k) =   reduce_bounds[2*k];
      upper_bounds(k) = - reduce_bounds[2*k + 1];
    }
  }
  else {
    this->lower_bounds = this->local_lower_bounds;
    this->upper_bounds = this->local_upper_bounds;
  }

  size = upper_bounds - lower_bounds;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<typename T>
void Mesh::initElementTypeMapArray(ElementTypeMapArray<T> & vect,
                                  UInt nb_component,
                                  UInt dim,
                                  const bool & flag_nb_node_per_elem_multiply,
                                  ElementKind element_kind,
                                  bool size_to_nb_element) const {
  AKANTU_DEBUG_IN();

  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;
    this->initElementTypeMapArray(vect, nb_component, dim, gt,
                                 flag_nb_node_per_elem_multiply,
                                 element_kind, size_to_nb_element);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<typename T>
void Mesh::initElementTypeMapArray(ElementTypeMapArray<T> & vect,
                                  UInt nb_component,
                                  UInt dim,
                                  GhostType gt,
                                  const bool & flag_nb_node_per_elem_multiply,
                                  ElementKind element_kind,
                                  bool size_to_nb_element) const {
  AKANTU_DEBUG_IN();

  Mesh::type_iterator it  = firstType(dim, gt, element_kind);
  Mesh::type_iterator end = lastType(dim, gt, element_kind);
  for(; it != end; ++it) {
    ElementType type = *it;
    if (flag_nb_node_per_elem_multiply) nb_component *= Mesh::getNbNodesPerElement(*it);
    UInt size = 0;
    if (size_to_nb_element) size = this->getNbElement(type, gt);
    vect.alloc(size, nb_component, type, gt);
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Mesh::initNormals() {
  initElementTypeMapArray(normals, spatial_dimension, spatial_dimension, false, _ek_not_defined);
}

/* -------------------------------------------------------------------------- */
void Mesh::getGlobalConnectivity(ElementTypeMapArray<UInt> & global_connectivity,
				 UInt dimension,
				 GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Mesh::type_iterator it  = firstType(dimension, ghost_type);
  Mesh::type_iterator end = lastType(dimension, ghost_type);

  for(; it != end; ++it) {
    ElementType type = *it;

    Array<UInt> & local_conn = connectivities(type, ghost_type);
    Array<UInt> & g_connectivity = global_connectivity(type, ghost_type);

    if (!nodes_global_ids)
      nodes_global_ids = mesh_parent->nodes_global_ids;

    UInt * local_c = local_conn.storage();
    UInt * global_c = g_connectivity.storage();

    UInt nb_terms = local_conn.getSize() * local_conn.getNbComponent();

    for (UInt i = 0; i < nb_terms; ++i, ++local_c, ++global_c)
      *global_c = (*nodes_global_ids)(*local_c);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

DumperIOHelper & Mesh::getGroupDumper(const std::string & dumper_name, 
				      const std::string & group_name){
  
  if (group_name == "all") return this->getDumper(dumper_name);
  else return element_groups[group_name]->getDumper(dumper_name);

}



/* -------------------------------------------------------------------------- */
#define AKANTU_INSTANTIATE_INIT(type)					\
template void Mesh::initElementTypeMapArray<type>(ElementTypeMapArray<type> & vect, \
                                                 UInt nb_component,	\
                                                 UInt dim,		\
                                                 const bool & flag_nb_elem_multiply, \
                                                 ElementKind element_kind, \
                                                 bool size_to_nb_element) const; \
template void Mesh::initElementTypeMapArray<type>(ElementTypeMapArray<type> & vect, \
                                                 UInt nb_component,	\
                                                 UInt dim,		\
						 GhostType gt,		\
                                                 const bool & flag_nb_elem_multiply, \
                                                 ElementKind element_kind, \
                                                 bool size_to_nb_element) const

AKANTU_INSTANTIATE_INIT(Real);
AKANTU_INSTANTIATE_INIT(UInt);
AKANTU_INSTANTIATE_INIT(Int);
AKANTU_INSTANTIATE_INIT(bool);

/* -------------------------------------------------------------------------- */
template <typename T>
ElementTypeMap<UInt> Mesh::getNbDataPerElem(ElementTypeMapArray<T> & array,
					    const ElementKind & element_kind){

  ElementTypeMap<UInt> nb_data_per_elem;
  
  typename ElementTypeMapArray<T>::type_iterator it = 
    array.firstType(spatial_dimension, _not_ghost,element_kind);
  typename ElementTypeMapArray<T>::type_iterator last_type = 
    array.lastType(spatial_dimension, _not_ghost,element_kind);
  
  for(; it != last_type; ++it) {
    UInt nb_elements = this->getNbElement(*it);
    nb_data_per_elem(*it) = 
      array(*it).getNbComponent() * 
      array(*it).getSize();

    nb_data_per_elem(*it) /= nb_elements;
  }
  return nb_data_per_elem;
}

/* -------------------------------------------------------------------------- */

template 
ElementTypeMap<UInt> Mesh::getNbDataPerElem(ElementTypeMapArray<Real> & array,
					    const ElementKind & element_kind);

template 
ElementTypeMap<UInt> Mesh::getNbDataPerElem(ElementTypeMapArray<UInt> & array,
					    const ElementKind & element_kind);


/* -------------------------------------------------------------------------- */

#ifdef AKANTU_USE_IOHELPER

template <typename T>
dumper::Field * Mesh::createFieldFromAttachedData(const std::string & field_id,
						  const std::string & group_name,
						  const ElementKind & element_kind){


  dumper::Field * field = NULL;
  ElementTypeMapArray<T> * internal = NULL;
  try {
    internal = &(this->getData<T>(field_id));
  }
  catch (...){
    return NULL;
  }

  ElementTypeMap<UInt> nb_data_per_elem = 
    this->getNbDataPerElem(*internal,element_kind);

  field = 
    this->createElementalField<T, dumper::InternalMaterialField>(*internal,
								group_name,
								this->spatial_dimension,
								element_kind,
								nb_data_per_elem);

  return field;
}

template dumper::Field * 
Mesh::createFieldFromAttachedData<Real>(const std::string & field_id,
					const std::string & group_name,
					const ElementKind & element_kind);

template dumper::Field * 
Mesh::createFieldFromAttachedData<UInt>(const std::string & field_id,
					const std::string & group_name,
					const ElementKind & element_kind);


/* -------------------------------------------------------------------------- */

#endif


__END_AKANTU__
