/**
 * @file   mesh.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Fri Jan 22 2016
 *
 * @brief  class handling meshes
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
#include "aka_config.hh"
/* -------------------------------------------------------------------------- */
#include "element_class.hh"
#include "group_manager_inline_impl.cc"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_utils.hh"
/* -------------------------------------------------------------------------- */
#include "element_synchronizer.hh"
#include "facet_synchronizer.hh"
#include "mesh_utils_distribution.hh"
#include "node_synchronizer.hh"
#include "communicator.hh"
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#include "dumper_field.hh"
#include "dumper_internal_material_field.hh"
#endif
/* -------------------------------------------------------------------------- */
#include <sstream>
/* -------------------------------------------------------------------------- */


namespace akantu {

/* -------------------------------------------------------------------------- */
Mesh::Mesh(UInt spatial_dimension, const ID & id, const MemoryID & memory_id,
           Communicator & communicator)
    : Memory(id, memory_id),
      GroupManager(*this, id + ":group_manager", memory_id),
      nodes_type(0, 1, id + ":nodes_type"),
      connectivities("connectivities", id, memory_id),
      ghosts_counters("ghosts_counters", id, memory_id),
      normals("normals", id, memory_id), spatial_dimension(spatial_dimension),
      lower_bounds(spatial_dimension, 0.), upper_bounds(spatial_dimension, 0.),
      size(spatial_dimension, 0.), local_lower_bounds(spatial_dimension, 0.),
      local_upper_bounds(spatial_dimension, 0.),
      mesh_data("mesh_data", id, memory_id), communicator(&communicator) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Mesh::Mesh(UInt spatial_dimension, Communicator & communicator,
           const ID & id, const MemoryID & memory_id)
    : Mesh(spatial_dimension, id, memory_id, communicator) {
  AKANTU_DEBUG_IN();

  this->nodes =
      std::make_shared<Array<Real>>(0, spatial_dimension, id + ":coordinates");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Mesh::Mesh(UInt spatial_dimension, const ID & id, const MemoryID & memory_id)
    : Mesh(spatial_dimension, Communicator::getStaticCommunicator(), id,
           memory_id) {}

/* -------------------------------------------------------------------------- */
Mesh::Mesh(UInt spatial_dimension, std::shared_ptr<Array<Real>> nodes,
           const ID & id, const MemoryID & memory_id)
    : Mesh(spatial_dimension, id, memory_id,
           Communicator::getStaticCommunicator()) {
  this->nodes = nodes;

  this->nb_global_nodes = this->nodes->size();

  this->nodes_to_elements.resize(nodes->size());
  for (auto & node_set : nodes_to_elements) {
    node_set = std::make_unique<std::set<Element>>();
  }

  this->computeBoundingBox();
}

/* -------------------------------------------------------------------------- */
Mesh & Mesh::initMeshFacets(const ID & id) {
  AKANTU_DEBUG_IN();

  if (!mesh_facets) {
    mesh_facets = std::make_unique<Mesh>(spatial_dimension, this->nodes,
                                         getID() + ":" + id, getMemoryID());

    mesh_facets->mesh_parent = this;
    mesh_facets->is_mesh_facets = true;

    MeshUtils::buildAllFacets(*this, *mesh_facets, 0);

    if (mesh.isDistributed()) {
      mesh_facets->is_distributed = true;
      mesh_facets->element_synchronizer = std::make_unique<FacetSynchronizer>(
          *mesh_facets, mesh.getElementSynchronizer());
    }
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
Mesh::~Mesh() = default;

/* -------------------------------------------------------------------------- */
void Mesh::read(const std::string & filename, const MeshIOType & mesh_io_type) {
  MeshIO mesh_io;
  mesh_io.read(filename, *this, mesh_io_type);

  type_iterator it =
      this->firstType(spatial_dimension, _not_ghost, _ek_not_defined);
  type_iterator last =
      this->lastType(spatial_dimension, _not_ghost, _ek_not_defined);
  if (it == last)
    AKANTU_EXCEPTION(
        "The mesh contained in the file "
        << filename << " does not seem to be of the good dimension."
        << " No element of dimension " << spatial_dimension << " where read.");

  this->nb_global_nodes = this->nodes->size();

  this->computeBoundingBox();

  this->nodes_to_elements.resize(nodes->size());
  for (auto & node_set : nodes_to_elements) {
    node_set = std::make_unique<std::set<Element>>();
  }
}

/* -------------------------------------------------------------------------- */
void Mesh::write(const std::string & filename,
                 const MeshIOType & mesh_io_type) {
  MeshIO mesh_io;
  mesh_io.write(filename, *this, mesh_io_type);
}

/* -------------------------------------------------------------------------- */
void Mesh::printself(std::ostream & stream, int indent) const {
  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
    ;

  stream << space << "Mesh [" << std::endl;
  stream << space << " + id                : " << getID() << std::endl;
  stream << space << " + spatial dimension : " << this->spatial_dimension
         << std::endl;
  stream << space << " + nodes [" << std::endl;
  nodes->printself(stream, indent + 2);
  stream << space << " + connectivities [" << std::endl;
  connectivities.printself(stream, indent + 2);
  stream << space << " ]" << std::endl;

  GroupManager::printself(stream, indent + 1);
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
void Mesh::computeBoundingBox() {
  AKANTU_DEBUG_IN();
  for (UInt k = 0; k < spatial_dimension; ++k) {
    local_lower_bounds(k) = std::numeric_limits<double>::max();
    local_upper_bounds(k) = -std::numeric_limits<double>::max();
  }

  for (UInt i = 0; i < nodes->size(); ++i) {
    //    if(!isPureGhostNode(i))
    for (UInt k = 0; k < spatial_dimension; ++k) {
      local_lower_bounds(k) = std::min(local_lower_bounds[k], (*nodes)(i, k));
      local_upper_bounds(k) = std::max(local_upper_bounds[k], (*nodes)(i, k));
    }
  }

  if (this->is_distributed) {
    Matrix<Real> reduce_bounds(spatial_dimension, 2);
    for (UInt k = 0; k < spatial_dimension; ++k) {
      reduce_bounds(k, 0) = local_lower_bounds(k);
      reduce_bounds(k, 1) = -local_upper_bounds(k);
    }

    communicator->allReduce(reduce_bounds, SynchronizerOperation::_min);

    for (UInt k = 0; k < spatial_dimension; ++k) {
      lower_bounds(k) = reduce_bounds(k, 0);
      upper_bounds(k) = -reduce_bounds(k, 1);
    }
  } else {
    this->lower_bounds = this->local_lower_bounds;
    this->upper_bounds = this->local_upper_bounds;
  }

  size = upper_bounds - lower_bounds;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Mesh::initNormals() {
  normals.initialize(*this, _nb_component = spatial_dimension,
                     _spatial_dimension = spatial_dimension,
                     _element_kind = _ek_not_defined);
}

/* -------------------------------------------------------------------------- */
void Mesh::getGlobalConnectivity(
    ElementTypeMapArray<UInt> & global_connectivity) {
  AKANTU_DEBUG_IN();

  for (auto && ghost_type : ghost_types) {
    for (auto type :
         global_connectivity.elementTypes(_ghost_type = ghost_type)) {
      if (not connectivities.exists(type, ghost_type))
        continue;

      auto & local_conn = connectivities(type, ghost_type);
      auto & g_connectivity = global_connectivity(type, ghost_type);

      UInt nb_nodes = local_conn.size() * local_conn.getNbComponent();

      if (not nodes_global_ids && is_mesh_facets) {
        std::transform(
            local_conn.begin_reinterpret(nb_nodes),
            local_conn.end_reinterpret(nb_nodes),
            g_connectivity.begin_reinterpret(nb_nodes),
            [& node_ids = *mesh_parent->nodes_global_ids](UInt l)->UInt {
              return node_ids(l);
            });
      } else {
        std::transform(local_conn.begin_reinterpret(nb_nodes),
                       local_conn.end_reinterpret(nb_nodes),
                       g_connectivity.begin_reinterpret(nb_nodes),
                       [& node_ids = *nodes_global_ids](UInt l)->UInt {
                         return node_ids(l);
                       });
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
DumperIOHelper & Mesh::getGroupDumper(const std::string & dumper_name,
                                      const std::string & group_name) {
  if (group_name == "all")
    return this->getDumper(dumper_name);
  else
    return element_groups[group_name]->getDumper(dumper_name);
}

/* -------------------------------------------------------------------------- */
template <typename T>
ElementTypeMap<UInt> Mesh::getNbDataPerElem(ElementTypeMapArray<T> & arrays,
                                            const ElementKind & element_kind) {
  ElementTypeMap<UInt> nb_data_per_elem;

  for (auto type : elementTypes(spatial_dimension, _not_ghost, element_kind)) {
    UInt nb_elements = this->getNbElement(type);
    auto & array = arrays(type);

    nb_data_per_elem(type) = array.getNbComponent() * array.size();
    nb_data_per_elem(type) /= nb_elements;
  }

  return nb_data_per_elem;
}

/* -------------------------------------------------------------------------- */
template ElementTypeMap<UInt>
Mesh::getNbDataPerElem(ElementTypeMapArray<Real> & array,
                       const ElementKind & element_kind);

template ElementTypeMap<UInt>
Mesh::getNbDataPerElem(ElementTypeMapArray<UInt> & array,
                       const ElementKind & element_kind);

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
template <typename T>
dumper::Field *
Mesh::createFieldFromAttachedData(const std::string & field_id,
                                  const std::string & group_name,
                                  const ElementKind & element_kind) {

  dumper::Field * field = nullptr;
  ElementTypeMapArray<T> * internal = nullptr;
  try {
    internal = &(this->getData<T>(field_id));
  } catch (...) {
    return nullptr;
  }

  ElementTypeMap<UInt> nb_data_per_elem =
      this->getNbDataPerElem(*internal, element_kind);

  field = this->createElementalField<T, dumper::InternalMaterialField>(
      *internal, group_name, this->spatial_dimension, element_kind,
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
#endif

/* -------------------------------------------------------------------------- */
void Mesh::distribute() {
  this->distribute(Communicator::getStaticCommunicator());
}

/* -------------------------------------------------------------------------- */
void Mesh::distribute(Communicator & communicator) {
  AKANTU_DEBUG_ASSERT(is_distributed == false,
                      "This mesh is already distribute");
  this->communicator = &communicator;

  this->element_synchronizer = std::make_unique<ElementSynchronizer>(
      *this, this->getID() + ":element_synchronizer", this->getMemoryID(),
      true);

  this->node_synchronizer = std::make_unique<NodeSynchronizer>(
      *this, this->getID() + ":node_synchronizer", this->getMemoryID(), true);

  Int psize = this->communicator->getNbProc();
#ifdef AKANTU_USE_SCOTCH
  Int prank = this->communicator->whoAmI();
  if (prank == 0) {
    MeshPartitionScotch partition(*this, spatial_dimension);
    partition.partitionate(psize);

    MeshUtilsDistribution::distributeMeshCentralized(*this, 0, partition);
  } else {
    MeshUtilsDistribution::distributeMeshCentralized(*this, 0);
  }
#else
  if (!(psize == 1)) {
    AKANTU_DEBUG_ERROR("Cannot distribute a mesh without a partitioning tool");
  }
#endif

  this->is_distributed = true;
  this->computeBoundingBox();
}

/* -------------------------------------------------------------------------- */
void Mesh::getAssociatedElements(const Array<UInt> & node_list,
                                 Array<Element> & elements) {
  for (const auto & node : node_list)
    for (const auto & element : *nodes_to_elements[node])
      elements.push_back(element);
}

/* -------------------------------------------------------------------------- */
void Mesh::fillNodesToElements() {
  Element e;

  UInt nb_nodes = nodes->size();
  for (UInt n = 0; n < nb_nodes; ++n) {
    if (this->nodes_to_elements[n])
      this->nodes_to_elements[n]->clear();
    else
      this->nodes_to_elements[n] = std::make_unique<std::set<Element>>();
  }

  for (auto ghost_type : ghost_types) {
    e.ghost_type = ghost_type;
    for (const auto & type :
         elementTypes(spatial_dimension, ghost_type, _ek_not_defined)) {
      e.type = type;

      UInt nb_element = this->getNbElement(type, ghost_type);
      Array<UInt>::const_iterator<Vector<UInt>> conn_it =
          connectivities(type, ghost_type)
              .begin(Mesh::getNbNodesPerElement(type));

      for (UInt el = 0; el < nb_element; ++el, ++conn_it) {
        e.element = el;
        const Vector<UInt> & conn = *conn_it;
        for (UInt n = 0; n < conn.size(); ++n)
          nodes_to_elements[conn(n)]->insert(e);
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void Mesh::eraseElements(const Array<Element> & elements) {
  ElementTypeMap<UInt> last_element;

  RemovedElementsEvent event(*this);
  auto & remove_list = event.getList();
  auto & new_numbering = event.getNewNumbering();

  for (auto && el : elements) {
    if (el.ghost_type != _not_ghost) {
      auto & count = ghosts_counters(el);
      --count;
      if (count > 0)
        continue;
    }

    remove_list.push_back(el);
    if (not last_element.exists(el.type, el.ghost_type)) {
      UInt nb_element = mesh.getNbElement(el.type, el.ghost_type);
      last_element(nb_element - 1, el.type, el.ghost_type);
      auto & numbering =
          new_numbering.alloc(nb_element, 1, el.type, el.ghost_type);
      for (auto && pair : enumerate(numbering)) {
        std::get<1>(pair) = std::get<0>(pair);
      }
    }

    UInt & pos = last_element(el.type, el.ghost_type);
    auto & numbering = new_numbering(el.type, el.ghost_type);

    numbering(el.element) = UInt(-1);
    numbering(pos) = el.element;
    --pos;
  }

  this->sendEvent(event);
}

} // namespace akantu
