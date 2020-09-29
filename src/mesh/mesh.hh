/**
 * @file   mesh.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Dana Christen <dana.christen@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Mon Feb 19 2018
 *
 * @brief  the class representing the meshes
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
#ifndef AKANTU_MESH_HH_
#define AKANTU_MESH_HH_

/* -------------------------------------------------------------------------- */
#include "aka_array.hh"
#include "aka_bbox.hh"
#include "aka_event_handler_manager.hh"
#include "aka_memory.hh"
#include "communicator.hh"
#include "dumpable.hh"
#include "element.hh"
#include "element_class.hh"
#include "element_type_map.hh"
#include "group_manager.hh"
#include "mesh_data.hh"
#include "mesh_events.hh"
/* -------------------------------------------------------------------------- */
#include <functional>
#include <set>
#include <unordered_map>
/* -------------------------------------------------------------------------- */

namespace akantu {
class ElementSynchronizer;
class NodeSynchronizer;
class PeriodicNodeSynchronizer;
class MeshGlobalDataUpdater;
} // namespace akantu

namespace akantu {

namespace {
  DECLARE_NAMED_ARGUMENT(communicator);
  DECLARE_NAMED_ARGUMENT(edge_weight_function);
  DECLARE_NAMED_ARGUMENT(vertex_weight_function);
} // namespace

/* -------------------------------------------------------------------------- */
/* Mesh                                                                       */
/* -------------------------------------------------------------------------- */

/**
 * @class  Mesh mesh.hh
 *
 * This class contaisn the coordinates of the nodes in the Mesh.nodes
 * akant::Array, and the connectivity. The connectivity are stored in by element
 * types.
 *
 * In order to loop on all element you have to loop on all types like this :
 * @code{.cpp}
 for(auto & type : mesh.elementTypes()) {
   UInt nb_element  = mesh.getNbElement(type);
   const Array<UInt> & conn = mesh.getConnectivity(type);
   for(UInt e = 0; e < nb_element; ++e) {
     ...
   }
 }

 or

 for_each_element(mesh, [](Element & element) {
    std::cout << element << std::endl
  });
 @endcode
*/
class Mesh : protected Memory,
             public EventHandlerManager<MeshEventHandler>,
             public GroupManager,
             public MeshData,
             public Dumpable {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
private:
  /// default constructor used for chaining, the last parameter is just to
  /// differentiate constructors
  Mesh(UInt spatial_dimension, const ID & id, const MemoryID & memory_id,
       Communicator & communicator);

public:
  /// constructor that create nodes coordinates array
  Mesh(UInt spatial_dimension, const ID & id = "mesh",
       const MemoryID & memory_id = 0);

  /// mesh not distributed and not using the default communicator
  Mesh(UInt spatial_dimension, Communicator & communicator,
       const ID & id = "mesh", const MemoryID & memory_id = 0);

  /**
   * constructor that use an existing nodes coordinates
   * array, by getting the vector of coordinates
   */
  Mesh(UInt spatial_dimension, const std::shared_ptr<Array<Real>> & nodes,
       const ID & id = "mesh", const MemoryID & memory_id = 0);

  ~Mesh() override;

  /// read the mesh from a file
  void read(const std::string & filename,
            const MeshIOType & mesh_io_type = _miot_auto);
  /// write the mesh to a file
  void write(const std::string & filename,
             const MeshIOType & mesh_io_type = _miot_auto);

protected:
  void makeReady();

private:
  /// initialize the connectivity to NULL and other stuff
  void init();

  /// function that computes the bounding box (fills xmin, xmax)
  void computeBoundingBox();

  /* ------------------------------------------------------------------------ */
  /* Distributed memory methods and accessors                                 */
  /* ------------------------------------------------------------------------ */
public:
protected:
  /// patitionate the mesh among the processors involved in their computation
  virtual void distributeImpl(
      Communicator & communicator,
      const std::function<Int(const Element &, const Element &)> &
          edge_weight_function,
      const std::function<Int(const Element &)> & vertex_weight_function);

public:
  /// with the arguments to pass to the partitionner
  template <typename... pack>
  std::enable_if_t<are_named_argument<pack...>::value>
  distribute(pack &&... _pack) {
    distributeImpl(
        OPTIONAL_NAMED_ARG(communicator, Communicator::getStaticCommunicator()),
        OPTIONAL_NAMED_ARG(edge_weight_function,
                           [](auto &&, auto &&) { return 1; }),
        OPTIONAL_NAMED_ARG(vertex_weight_function, [](auto &&) { return 1; }));
  }

  /// defines is the mesh is distributed or not
  inline bool isDistributed() const { return this->is_distributed; }

  /* ------------------------------------------------------------------------ */
  /* Periodicity methods and accessors                                        */
  /* ------------------------------------------------------------------------ */
public:
  /// set the periodicity in a given direction
  void makePeriodic(const SpatialDirection & direction);
  void makePeriodic(const SpatialDirection & direction, const ID & list_1,
                    const ID & list_2);

protected:
  void makePeriodic(const SpatialDirection & direction,
                    const Array<UInt> & list_left,
                    const Array<UInt> & list_right);

  /// Removes the face that the mesh is periodic
  void wipePeriodicInfo();

  inline void addPeriodicSlave(UInt slave, UInt master);

  template <typename T>
  void synchronizePeriodicSlaveDataWithMaster(Array<T> & data);

  // update the periodic synchronizer (creates it if it does not exists)
  void updatePeriodicSynchronizer();

public:
  /// defines if the mesh is periodic or not
  inline bool isPeriodic() const { return this->is_periodic; }

  inline bool isPeriodic(const SpatialDirection & /*direction*/) const {
    return this->is_periodic;
  }

  class PeriodicSlaves;

  /// get the master node for a given slave nodes, except if node not a slave
  inline UInt getPeriodicMaster(UInt slave) const;

  /// get an iterable list of slaves for a given master node
  inline decltype(auto) getPeriodicSlaves(UInt master) const;

  /* ------------------------------------------------------------------------ */
  /* General Methods                                                          */
  /* ------------------------------------------------------------------------ */
public:
  /// function to print the containt of the class
  void printself(std::ostream & stream, int indent = 0) const override;

  /// extract coordinates of nodes from an element
  template <typename T>
  inline void
  extractNodalValuesFromElement(const Array<T> & nodal_values, T * local_coord,
                                const UInt * connectivity, UInt n_nodes,
                                UInt nb_degree_of_freedom) const;

  // /// extract coordinates of nodes from a reversed element
  // inline void extractNodalCoordinatesFromPBCElement(Real * local_coords,
  //                                                   UInt * connectivity,
  //                                                   UInt n_nodes);

  /// add a Array of connectivity for the given ElementType and GhostType .
  inline void addConnectivityType(ElementType type,
                                  GhostType ghost_type = _not_ghost);

  /* ------------------------------------------------------------------------ */
  template <class Event> inline void sendEvent(Event & event) {
    //    if(event.getList().size() != 0)
    EventHandlerManager<MeshEventHandler>::sendEvent<Event>(event);
  }

  /// prepare the  event to remove the elements listed
  void eraseElements(const Array<Element> & elements);

  /* ------------------------------------------------------------------------ */
  template <typename T>
  inline void removeNodesFromArray(Array<T> & vect,
                                   const Array<UInt> & new_numbering);

  /// initialize normals
  void initNormals();

  /// init facets' mesh
  Mesh & initMeshFacets(const ID & id = "mesh_facets");

  /// define parent mesh
  void defineMeshParent(const Mesh & mesh);

  /// get global connectivity array
  void getGlobalConnectivity(ElementTypeMapArray<UInt> & global_connectivity);

public:
  void getAssociatedElements(const Array<UInt> & node_list,
                             Array<Element> & elements);

private:
  /// fills the nodes_to_elements structure
  void fillNodesToElements();

  /// update the global ids, nodes type, ...
  std::tuple<UInt, UInt> updateGlobalData(NewNodesEvent & nodes_event,
                                          NewElementsEvent & elements_event);

  void registerGlobalDataUpdater(
      std::unique_ptr<MeshGlobalDataUpdater> && global_data_updater);
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the id of the mesh
  AKANTU_GET_MACRO(ID, Memory::id, const ID &);

  /// get the id of the mesh
  AKANTU_GET_MACRO(MemoryID, Memory::memory_id, const MemoryID &);

  /// get the spatial dimension of the mesh = number of component of the
  /// coordinates
  AKANTU_GET_MACRO(SpatialDimension, spatial_dimension, UInt);

  /// get the nodes Array aka coordinates
  AKANTU_GET_MACRO(Nodes, *nodes, const Array<Real> &);
  AKANTU_GET_MACRO_NOT_CONST(Nodes, *nodes, Array<Real> &);

  /// get the normals for the elements
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(Normals, normals, Real);

  /// get the number of nodes
  AKANTU_GET_MACRO(NbNodes, nodes->size(), UInt);

  /// get the Array of global ids of the nodes (only used in parallel)
  AKANTU_GET_MACRO(GlobalNodesIds, *nodes_global_ids, const Array<UInt> &);
  // AKANTU_GET_MACRO_NOT_CONST(GlobalNodesIds, *nodes_global_ids, Array<UInt>
  // &);

  /// get the global id of a node
  inline UInt getNodeGlobalId(UInt local_id) const;

  /// get the global id of a node
  inline UInt getNodeLocalId(UInt global_id) const;

  /// get the global number of nodes
  inline UInt getNbGlobalNodes() const;

  /// get the nodes type Array
  AKANTU_GET_MACRO(NodesFlags, *nodes_flags, const Array<NodeFlag> &);

protected:
  AKANTU_GET_MACRO_NOT_CONST(NodesFlags, *nodes_flags, Array<NodeFlag> &);

public:
  inline NodeFlag getNodeFlag(UInt local_id) const;
  inline Int getNodePrank(UInt local_id) const;

  /// say if a node is a pure ghost node
  inline bool isPureGhostNode(UInt n) const;

  /// say if a node is pur local or master node
  inline bool isLocalOrMasterNode(UInt n) const;

  inline bool isLocalNode(UInt n) const;
  inline bool isMasterNode(UInt n) const;
  inline bool isSlaveNode(UInt n) const;

  inline bool isPeriodicSlave(UInt n) const;
  inline bool isPeriodicMaster(UInt n) const;

  const Vector<Real> & getLowerBounds() const { return bbox.getLowerBounds(); }
  const Vector<Real> & getUpperBounds() const { return bbox.getUpperBounds(); }
  AKANTU_GET_MACRO(BBox, bbox, const BBox &);

  const Vector<Real> & getLocalLowerBounds() const {
    return bbox_local.getLowerBounds();
  }
  const Vector<Real> & getLocalUpperBounds() const {
    return bbox_local.getUpperBounds();
  }
  AKANTU_GET_MACRO(LocalBBox, bbox_local, const BBox &);

  /// get the connectivity Array for a given type
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Connectivity, connectivities, UInt);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(Connectivity, connectivities, UInt);
  AKANTU_GET_MACRO(Connectivities, connectivities,
                   const ElementTypeMapArray<UInt> &);

  /// get the number of element of a type in the mesh
  inline UInt getNbElement(ElementType type,
                           GhostType ghost_type = _not_ghost) const;

  /// get the number of element for a given ghost_type and a given dimension
  inline UInt getNbElement(UInt spatial_dimension = _all_dimensions,
                           GhostType ghost_type = _not_ghost,
                           ElementKind kind = _ek_not_defined) const;

  /// compute the barycenter of a given element
  inline void getBarycenter(const Element & element,
                            Vector<Real> & barycenter) const;

  void getBarycenters(Array<Real> & barycenter, ElementType type,
                      GhostType ghost_type) const;

  /// get the element connected to a subelement (element of lower dimension)
  const auto & getElementToSubelement() const;

  /// get the element connected to a subelement
  const auto & getElementToSubelement(ElementType el_type,
                                      GhostType ghost_type = _not_ghost) const;

  /// get the element connected to a subelement
  auto & getElementToSubelement(ElementType el_type,
                                GhostType ghost_type = _not_ghost);

  /// get the elements connected to a subelement
  const auto & getElementToSubelement(const Element & element) const;

  /// get the subelement (element of lower dimension) connected to a element
  const auto & getSubelementToElement() const;

  /// get the subelement connected to an element
  const auto & getSubelementToElement(ElementType el_type,
                                      GhostType ghost_type = _not_ghost) const;
  /// get the subelement connected to an element
  auto & getSubelementToElement(ElementType el_type,
                                GhostType ghost_type = _not_ghost);

  /// get the subelement (element of lower dimension) connected to a element
  VectorProxy<Element> getSubelementToElement(const Element & element) const;

  /// get connectivity of a given element
  inline VectorProxy<UInt> getConnectivity(const Element & element) const;
  inline Vector<UInt>
  getConnectivityWithPeriodicity(const Element & element) const;

protected:
  inline auto & getElementToSubelement(const Element & element);
  inline VectorProxy<Element> getSubelementToElement(const Element & element);
  inline VectorProxy<UInt> getConnectivity(const Element & element);

public:
  /// get a name field associated to the mesh
  template <typename T>
  inline const Array<T> & getData(const ID & data_name, ElementType el_type,
                                  GhostType ghost_type = _not_ghost) const;

  /// get a name field associated to the mesh
  template <typename T>
  inline Array<T> & getData(const ID & data_name, ElementType el_type,
                            GhostType ghost_type = _not_ghost);

  /// get a name field associated to the mesh
  template <typename T>
  inline const ElementTypeMapArray<T> & getData(const ID & data_name) const;

  /// get a name field associated to the mesh
  template <typename T>
  inline ElementTypeMapArray<T> & getData(const ID & data_name);

  template <typename T>
  ElementTypeMap<UInt> getNbDataPerElem(ElementTypeMapArray<T> & array);

  template <typename T>
  std::shared_ptr<dumpers::Field>
  createFieldFromAttachedData(const std::string & field_id,
                              const std::string & group_name,
                              ElementKind element_kind);

  /// templated getter returning the pointer to data in MeshData (modifiable)
  template <typename T>
  inline Array<T> &
  getDataPointer(const std::string & data_name, ElementType el_type,
                 GhostType ghost_type = _not_ghost, UInt nb_component = 1,
                 bool size_to_nb_element = true,
                 bool resize_with_parent = false);

  template <typename T>
  inline Array<T> & getDataPointer(const ID & data_name, ElementType el_type,
                                   GhostType ghost_type, UInt nb_component,
                                   bool size_to_nb_element,
                                   bool resize_with_parent, const T & defaul_);

  /// Facets mesh accessor
  inline const Mesh & getMeshFacets() const;
  inline Mesh & getMeshFacets();

  /// Parent mesh accessor
  inline const Mesh & getMeshParent() const;

  inline bool isMeshFacets() const { return this->is_mesh_facets; }

  /// return the dumper from a group and and a dumper name
  DumperIOHelper & getGroupDumper(const std::string & dumper_name,
                                  const std::string & group_name);

  /* ------------------------------------------------------------------------ */
  /* Wrappers on ElementClass functions                                       */
  /* ------------------------------------------------------------------------ */
public:
  /// get the number of nodes per element for a given element type
  static inline UInt getNbNodesPerElement(ElementType type);

  /// get the number of nodes per element for a given element type considered as
  /// a first order element
  static inline ElementType getP1ElementType(ElementType type);

  /// get the kind of the element type
  static inline ElementKind getKind(ElementType type);

  /// get spatial dimension of a type of element
  static inline UInt getSpatialDimension(ElementType type);

  /// get number of facets of a given element type
  static inline UInt getNbFacetsPerElement(ElementType type);

  /// get number of facets of a given element type
  static inline UInt getNbFacetsPerElement(ElementType type, UInt t);

  /// get local connectivity of a facet for a given facet type
  static inline auto getFacetLocalConnectivity(ElementType type, UInt t = 0);

  /// get connectivity of facets for a given element
  inline auto getFacetConnectivity(const Element & element, UInt t = 0) const;

  /// get the number of type of the surface element associated to a given
  /// element type
  static inline UInt getNbFacetTypes(ElementType type, UInt t = 0);

  /// get the type of the surface element associated to a given element
  static inline constexpr auto getFacetType(ElementType type, UInt t = 0);

  /// get all the type of the surface element associated to a given element
  static inline constexpr auto getAllFacetTypes(ElementType type);

  /// get the number of nodes in the given element list
  static inline UInt getNbNodesPerElementList(const Array<Element> & elements);

  /* ------------------------------------------------------------------------ */
  /* Element type Iterator                                                    */
  /* ------------------------------------------------------------------------ */

  using type_iterator [[deprecated]] =
      ElementTypeMapArray<UInt, ElementType>::type_iterator;
  using ElementTypesIteratorHelper =
      ElementTypeMapArray<UInt, ElementType>::ElementTypesIteratorHelper;

  template <typename... pack>
  ElementTypesIteratorHelper elementTypes(pack &&... _pack) const;

  [[deprecated("Use elementTypes instead")]] inline decltype(auto)
  firstType(UInt dim = _all_dimensions, GhostType ghost_type = _not_ghost,
            ElementKind kind = _ek_regular) const {
    return connectivities.elementTypes(dim, ghost_type, kind).begin();
  }

  [[deprecated("Use elementTypes instead")]] inline decltype(auto)
  lastType(UInt dim = _all_dimensions, GhostType ghost_type = _not_ghost,
           ElementKind kind = _ek_regular) const {
    return connectivities.elementTypes(dim, ghost_type, kind).end();
  }

  AKANTU_GET_MACRO(ElementSynchronizer, *element_synchronizer,
                   const ElementSynchronizer &);
  AKANTU_GET_MACRO_NOT_CONST(ElementSynchronizer, *element_synchronizer,
                             ElementSynchronizer &);
  AKANTU_GET_MACRO(NodeSynchronizer, *node_synchronizer,
                   const NodeSynchronizer &);
  AKANTU_GET_MACRO_NOT_CONST(NodeSynchronizer, *node_synchronizer,
                             NodeSynchronizer &);
  AKANTU_GET_MACRO(PeriodicNodeSynchronizer, *periodic_node_synchronizer,
                   const PeriodicNodeSynchronizer &);
  AKANTU_GET_MACRO_NOT_CONST(PeriodicNodeSynchronizer,
                             *periodic_node_synchronizer,
                             PeriodicNodeSynchronizer &);

  // AKANTU_GET_MACRO_NOT_CONST(Communicator, *communicator, StaticCommunicator
  // &);
  AKANTU_GET_MACRO(Communicator, *communicator, const auto &);
  AKANTU_GET_MACRO_NOT_CONST(Communicator, *communicator, auto &);
  AKANTU_GET_MACRO(PeriodicMasterSlaves, periodic_master_slave, const auto &);

  /* ------------------------------------------------------------------------ */
  /* Private methods for friends                                              */
  /* ------------------------------------------------------------------------ */
private:
  friend class MeshAccessor;
  friend class MeshUtils;

  AKANTU_GET_MACRO(NodesPointer, *nodes, Array<Real> &);

  /// get a pointer to the nodes_global_ids Array<UInt> and create it if
  /// necessary
  inline Array<UInt> & getNodesGlobalIdsPointer();

  /// get a pointer to the nodes_type Array<Int> and create it if necessary
  inline Array<NodeFlag> & getNodesFlagsPointer();

  /// get a pointer to the connectivity Array for the given type and create it
  /// if necessary
  inline Array<UInt> &
  getConnectivityPointer(ElementType type, GhostType ghost_type = _not_ghost);

  /// get the ghost element counter
  inline Array<UInt> & getGhostsCounters(ElementType type,
                                         GhostType ghost_type = _ghost) {
    AKANTU_DEBUG_ASSERT(ghost_type != _not_ghost,
                        "No ghost counter for _not_ghost elements");
    return ghosts_counters(type, ghost_type);
  }

  /// get a pointer to the element_to_subelement Array for the given type and
  /// create it if necessary
  inline Array<std::vector<Element>> &
  getElementToSubelementPointer(ElementType type,
                                GhostType ghost_type = _not_ghost);

  /// get a pointer to the subelement_to_element Array for the given type and
  /// create it if necessary
  inline Array<Element> &
  getSubelementToElementPointer(ElementType type,
                                GhostType ghost_type = _not_ghost);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// array of the nodes coordinates
  std::shared_ptr<Array<Real>> nodes;

  /// global node ids
  std::shared_ptr<Array<UInt>> nodes_global_ids;

  /// node flags (shared/periodic/...)
  std::shared_ptr<Array<NodeFlag>> nodes_flags;

  /// processor handling the node when not local or master
  std::unordered_map<UInt, Int> nodes_prank;

  /// global number of nodes;
  UInt nb_global_nodes{0};

  /// all class of elements present in this mesh (for heterogenous meshes)
  ElementTypeMapArray<UInt> connectivities;

  /// count the references on ghost elements
  ElementTypeMapArray<UInt> ghosts_counters;

  /// map to normals for all class of elements present in this mesh
  ElementTypeMapArray<Real> normals;

  /// the spatial dimension of this mesh
  UInt spatial_dimension{0};

  /// size covered by the mesh on each direction
  Vector<Real> size;
  /// global bounding box
  BBox bbox;

  /// local bounding box
  BBox bbox_local;

  /// Extra data loaded from the mesh file
  // MeshData mesh_data;

  /// facets' mesh
  std::unique_ptr<Mesh> mesh_facets;

  /// parent mesh (this is set for mesh_facets meshes)
  const Mesh * mesh_parent{nullptr};

  /// defines if current mesh is mesh_facets or not
  bool is_mesh_facets{false};

  /// defines if the mesh is centralized or distributed
  bool is_distributed{false};

  /// defines if the mesh is periodic
  bool is_periodic{false};

  /// Communicator on which mesh is distributed
  Communicator * communicator;

  /// Element synchronizer
  std::unique_ptr<ElementSynchronizer> element_synchronizer;

  /// Node synchronizer
  std::unique_ptr<NodeSynchronizer> node_synchronizer;

  /// Node synchronizer for periodic nodes
  std::unique_ptr<PeriodicNodeSynchronizer> periodic_node_synchronizer;

  using NodesToElements = std::vector<std::unique_ptr<std::set<Element>>>;

  /// class to update global data using external knowledge
  std::unique_ptr<MeshGlobalDataUpdater> global_data_updater;

  /// This info is stored to simplify the dynamic changes
  NodesToElements nodes_to_elements;

  /// periodicity local info
  std::unordered_map<UInt, UInt> periodic_slave_master;
  std::unordered_multimap<UInt, UInt> periodic_master_slave;
};

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream, const Mesh & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* Inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "element_type_map_tmpl.hh"
#include "mesh_inline_impl.hh"

#endif /* AKANTU_MESH_HH_ */
