/**
 * @file   mesh.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  class handling meshes
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
#include "aka_config.hh"
/* -------------------------------------------------------------------------- */
#include "element_class.hh"
#include "group_manager_inline_impl.hh"
#include "mesh.hh"
#include "mesh_global_data_updater.hh"
#include "mesh_io.hh"
#include "mesh_iterators.hh"
#include "mesh_utils.hh"
/* -------------------------------------------------------------------------- */
#include "communicator.hh"
#include "element_synchronizer.hh"
#include "facet_synchronizer.hh"
#include "mesh_utils_distribution.hh"
#include "node_synchronizer.hh"
#include "periodic_node_synchronizer.hh"
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#include "dumper_field.hh"
#include "dumper_internal_material_field.hh"
#endif
/* -------------------------------------------------------------------------- */
#include <limits>
#include <sstream>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
Mesh::Mesh(UInt spatial_dimension, const ID & id, const MemoryID & memory_id,
           Communicator & communicator)
    : Memory(id, memory_id),
      GroupManager(*this, id + ":group_manager", memory_id),
      MeshData("mesh_data", id, memory_id),
      connectivities("connectivities", id, memory_id),
      ghosts_counters("ghosts_counters", id, memory_id),
      normals("normals", id, memory_id), spatial_dimension(spatial_dimension),
      size(spatial_dimension, 0.), bbox(spatial_dimension),
      bbox_local(spatial_dimension), communicator(&communicator) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Mesh::Mesh(UInt spatial_dimension, Communicator & communicator, const ID & id,
           const MemoryID & memory_id)
    : Mesh(spatial_dimension, id, memory_id, communicator) {
  AKANTU_DEBUG_IN();

  this->nodes =
      std::make_shared<Array<Real>>(0, spatial_dimension, id + ":coordinates");
  this->nodes_flags = std::make_shared<Array<NodeFlag>>(0, 1, NodeFlag::_normal,
                                                        id + ":nodes_flags");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Mesh::Mesh(UInt spatial_dimension, const ID & id, const MemoryID & memory_id)
    : Mesh(spatial_dimension, Communicator::getStaticCommunicator(), id,
           memory_id) {}

/* -------------------------------------------------------------------------- */
Mesh::Mesh(UInt spatial_dimension, const std::shared_ptr<Array<Real>> & nodes,
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
void Mesh::getBarycenters(Array<Real> & barycenter, const ElementType & type,
                          const GhostType & ghost_type) const {
  barycenter.resize(getNbElement(type, ghost_type));
  for (auto && data : enumerate(make_view(barycenter, spatial_dimension))) {
    getBarycenter(Element{type, UInt(std::get<0>(data)), ghost_type},
                  std::get<1>(data));
  }
}

class FacetGlobalConnectivityAccessor : public DataAccessor<Element> {
public:
  FacetGlobalConnectivityAccessor(Mesh & mesh)
      : global_connectivity("global_connectivity",
                            "facet_connectivity_synchronizer") {
    global_connectivity.initialize(
        mesh, _spatial_dimension = _all_dimensions, _with_nb_element = true,
        _with_nb_nodes_per_element = true, _element_kind = _ek_regular);
    mesh.getGlobalConnectivity(global_connectivity);
  }

  UInt getNbData(const Array<Element> & elements,
                 const SynchronizationTag & tag) const {
    UInt size = 0;
    if (tag == SynchronizationTag::_smmc_facets_conn) {
      UInt nb_nodes = Mesh::getNbNodesPerElementList(elements);
      size += nb_nodes * sizeof(UInt);
    }
    return size;
  }

  void packData(CommunicationBuffer & buffer, const Array<Element> & elements,
                const SynchronizationTag & tag) const {
    if (tag == SynchronizationTag::_smmc_facets_conn) {
      for (const auto & element : elements) {
        auto & conns = global_connectivity(element.type, element.ghost_type);
        for (auto n : arange(conns.getNbComponent())) {
          buffer << conns(element.element, n);
        }
      }
    }
  }

  void unpackData(CommunicationBuffer & buffer, const Array<Element> & elements,
                  const SynchronizationTag & tag) {
    if (tag == SynchronizationTag::_smmc_facets_conn) {
      for (const auto & element : elements) {
        auto & conns = global_connectivity(element.type, element.ghost_type);
        for (auto n : arange(conns.getNbComponent())) {
          buffer >> conns(element.element, n);
        }
      }
    }
  }

  AKANTU_GET_MACRO(GlobalConnectivity, (global_connectivity), decltype(auto));

protected:
  ElementTypeMapArray<UInt> global_connectivity;
};

/* -------------------------------------------------------------------------- */
Mesh & Mesh::initMeshFacets(const ID & id) {
  AKANTU_DEBUG_IN();

  if (mesh_facets) {
    AKANTU_DEBUG_OUT();
    return *mesh_facets;
  }

  mesh_facets = std::make_unique<Mesh>(spatial_dimension, this->nodes,
                                       getID() + ":" + id, getMemoryID());

  mesh_facets->mesh_parent = this;
  mesh_facets->is_mesh_facets = true;
  mesh_facets->nodes_flags = this->nodes_flags;
  mesh_facets->nodes_global_ids = this->nodes_global_ids;

  MeshUtils::buildAllFacets(*this, *mesh_facets, 0);

  if (mesh.isDistributed()) {
    mesh_facets->is_distributed = true;
    mesh_facets->element_synchronizer = std::make_unique<FacetSynchronizer>(
        *mesh_facets, mesh.getElementSynchronizer());

    FacetGlobalConnectivityAccessor data_accessor(*mesh_facets);
    /// communicate
    mesh_facets->element_synchronizer->synchronizeOnce(
        data_accessor, SynchronizationTag::_smmc_facets_conn);

    /// flip facets
    MeshUtils::flipFacets(*mesh_facets, data_accessor.getGlobalConnectivity(),
                          _ghost);
  }

  /// transfers the the mesh physical names to the mesh facets
  if (not this->hasData("physical_names")) {
    AKANTU_DEBUG_OUT();
    return *mesh_facets;
  }

  auto & mesh_phys_data = this->getData<std::string>("physical_names");
  auto & phys_data = mesh_facets->getData<std::string>("physical_names");
  phys_data.initialize(*mesh_facets, _spatial_dimension = spatial_dimension - 1,
                       _with_nb_element = true);

  ElementTypeMapArray<Real> barycenters(getID(), "temporary_barycenters");
  barycenters.initialize(*mesh_facets, _nb_component = spatial_dimension,
                         _spatial_dimension = spatial_dimension - 1,
                         _with_nb_element = true);

  for (auto && ghost_type : ghost_types) {
    for (auto && type :
         barycenters.elementTypes(spatial_dimension - 1, ghost_type)) {
      mesh_facets->getBarycenters(barycenters(type, ghost_type), type,
                                  ghost_type);
    }
  }

  for_each_element(
      mesh,
      [&](auto && element) {
        Vector<Real> barycenter(spatial_dimension);
        mesh.getBarycenter(element, barycenter);
        auto norm_barycenter = barycenter.norm();
        auto tolerance = Math::getTolerance();
        if (norm_barycenter > tolerance)
          tolerance *= norm_barycenter;

        const auto & element_to_facet = mesh_facets->getElementToSubelement(
            element.type, element.ghost_type);

        Vector<Real> barycenter_facet(spatial_dimension);

        auto range = enumerate(make_view(
            barycenters(element.type, element.ghost_type), spatial_dimension));
#ifndef AKANTU_NDEBUG
        auto min_dist = std::numeric_limits<Real>::max();
#endif
        // this is a spacial search coded the most inefficient way.
        auto facet =
            std::find_if(range.begin(), range.end(), [&](auto && data) {
              auto facet = std::get<0>(data);
              if (element_to_facet(facet)[1] == ElementNull)
                return false;

              auto norm_distance = barycenter.distance(std::get<1>(data));
#ifndef AKANTU_NDEBUG
              min_dist = std::min(min_dist, norm_distance);
#endif
              return (norm_distance < tolerance);
            });

        if (facet == range.end()) {
          AKANTU_DEBUG_INFO("The element "
                            << element
                            << " did not find its associated facet in the "
                               "mesh_facets! Try to decrease math tolerance. "
                               "The closest element was at a distance of "
                            << min_dist);
          return;
        }

        // set physical name
        phys_data(Element{element.type, UInt(std::get<0>(*facet)),
                          element.ghost_type}) = mesh_phys_data(element);
      },
      _spatial_dimension = spatial_dimension - 1);

  mesh_facets->createGroupsFromMeshData<std::string>("physical_names");

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

  AKANTU_DEBUG_ASSERT(not is_distributed,
                      "You cannot read a mesh that is already distributed");

  MeshIO mesh_io;
  mesh_io.read(filename, *this, mesh_io_type);

  auto types =
      this->elementTypes(spatial_dimension, _not_ghost, _ek_not_defined);
  auto it = types.begin();
  auto last = types.end();
  if (it == last) {
    AKANTU_DEBUG_WARNING(
        "The mesh contained in the file "
        << filename << " does not seem to be of the good dimension."
        << " No element of dimension " << spatial_dimension << " were read.");
  }

  this->makeReady();
}

/* -------------------------------------------------------------------------- */
void Mesh::write(const std::string & filename,
                 const MeshIOType & mesh_io_type) {
  MeshIO mesh_io;
  mesh_io.write(filename, *this, mesh_io_type);
}

/* -------------------------------------------------------------------------- */
void Mesh::makeReady() {
  this->nb_global_nodes = this->nodes->size();
  this->computeBoundingBox();
  this->nodes_flags->resize(nodes->size(), NodeFlag::_normal);
  this->nodes_to_elements.resize(nodes->size());
  for (auto & node_set : nodes_to_elements) {
    node_set = std::make_unique<std::set<Element>>();
  }
}

/* -------------------------------------------------------------------------- */
void Mesh::printself(std::ostream & stream, int indent) const {
  std::string space(indent, AKANTU_INDENT);

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

  bbox_local.reset();

  for (auto & pos : make_view(*nodes, spatial_dimension)) {
    //    if(!isPureGhostNode(i))
    bbox_local += pos;
  }

  if (this->is_distributed) {
    bbox = bbox_local.allSum(*communicator);
  } else {
    bbox = bbox_local;
  }

  size = bbox.size();

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
         global_connectivity.elementTypes(_spatial_dimension = _all_dimensions,
         _element_kind = _ek_not_defined, _ghost_type = ghost_type)) {
      if (not connectivities.exists(type, ghost_type))
        continue;

      auto & local_conn = connectivities(type, ghost_type);
      auto & g_connectivity = global_connectivity(type, ghost_type);

      UInt nb_nodes = local_conn.size() * local_conn.getNbComponent();

      std::transform(local_conn.begin_reinterpret(nb_nodes),
                     local_conn.end_reinterpret(nb_nodes),
                     g_connectivity.begin_reinterpret(nb_nodes),
                     [&](UInt l) -> UInt { return this->getNodeGlobalId(l); });
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
ElementTypeMap<UInt> Mesh::getNbDataPerElem(ElementTypeMapArray<T> & arrays) {
  ElementTypeMap<UInt> nb_data_per_elem;

  for (auto type : arrays.elementTypes()) {
    UInt nb_elements = this->getNbElement(type);
    auto & array = arrays(type);

    nb_data_per_elem(type) = array.getNbComponent() * array.size();
    nb_data_per_elem(type) /= nb_elements;
  }

  return nb_data_per_elem;
}

/* -------------------------------------------------------------------------- */
template ElementTypeMap<UInt>
Mesh::getNbDataPerElem(ElementTypeMapArray<Real> & array);

template ElementTypeMap<UInt>
Mesh::getNbDataPerElem(ElementTypeMapArray<UInt> & array);

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
template <typename T>
std::shared_ptr<dumpers::Field>
Mesh::createFieldFromAttachedData(const std::string & field_id,
                                  const std::string & group_name,
                                  const ElementKind & element_kind) {

  std::shared_ptr<dumpers::Field> field;
  ElementTypeMapArray<T> * internal = nullptr;
  try {
    internal = &(this->getData<T>(field_id));
  } catch (...) {
    return nullptr;
  }

  ElementTypeMap<UInt> nb_data_per_elem = this->getNbDataPerElem(*internal);

  field = this->createElementalField<T, dumpers::InternalMaterialField>(
      *internal, group_name, this->spatial_dimension, element_kind,
      nb_data_per_elem);

  return field;
}

template std::shared_ptr<dumpers::Field>
Mesh::createFieldFromAttachedData<Real>(const std::string & field_id,
                                        const std::string & group_name,
                                        const ElementKind & element_kind);

template std::shared_ptr<dumpers::Field>
Mesh::createFieldFromAttachedData<UInt>(const std::string & field_id,
                                        const std::string & group_name,
                                        const ElementKind & element_kind);
#endif

/* -------------------------------------------------------------------------- */
void Mesh::distributeImpl(
    Communicator & communicator,
    std::function<Int(const Element &, const Element &)> edge_weight_function
    [[gnu::unused]],
    std::function<Int(const Element &)> vertex_weight_function
    [[gnu::unused]]) {
  AKANTU_DEBUG_ASSERT(is_distributed == false,
                      "This mesh is already distribute");
  this->communicator = &communicator;

  this->element_synchronizer = std::make_unique<ElementSynchronizer>(
      *this, this->getID() + ":element_synchronizer", this->getMemoryID(),
      true);

  this->node_synchronizer = std::make_unique<NodeSynchronizer>(
      *this, this->getID() + ":node_synchronizer", this->getMemoryID(), true);

  Int psize = this->communicator->getNbProc();

  if (psize > 1) {
#ifdef AKANTU_USE_SCOTCH
    Int prank = this->communicator->whoAmI();
    if (prank == 0) {
      MeshPartitionScotch partition(*this, spatial_dimension);
      partition.partitionate(psize, edge_weight_function,
                             vertex_weight_function);

      MeshUtilsDistribution::distributeMeshCentralized(*this, 0, partition);
    } else {
      MeshUtilsDistribution::distributeMeshCentralized(*this, 0);
    }
#else
    if (psize > 1) {
      AKANTU_ERROR("Cannot distribute a mesh without a partitioning tool");
    }
#endif
  }

  // if (psize > 1)
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
std::tuple<UInt, UInt>
Mesh::updateGlobalData(NewNodesEvent & nodes_event,
                       NewElementsEvent & elements_event) {
  if (global_data_updater)
    return this->global_data_updater->updateData(nodes_event, elements_event);
  else {
    return std::make_tuple(nodes_event.getList().size(),
                           elements_event.getList().size());
  }
}

/* -------------------------------------------------------------------------- */
void Mesh::registerGlobalDataUpdater(
    std::unique_ptr<MeshGlobalDataUpdater> && global_data_updater) {
  this->global_data_updater = std::move(global_data_updater);
}

/* -------------------------------------------------------------------------- */
void Mesh::eraseElements(const Array<Element> & elements) {
  ElementTypeMap<UInt> last_element;

  RemovedElementsEvent event(*this, "new_numbering", AKANTU_CURRENT_FUNCTION);
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
