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
 * @date last modification: Thu Jan 14 2016
 *
 * @brief  the class representing the meshes
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
#ifndef __AKANTU_MESH_HH__
#define __AKANTU_MESH_HH__

/* -------------------------------------------------------------------------- */
#include "aka_array.hh"
#include "aka_event_handler_manager.hh"
#include "aka_memory.hh"
#include "dumpable.hh"
#include "element.hh"
#include "element_class.hh"
#include "element_type_map.hh"
#include "group_manager.hh"
#include "mesh_data.hh"
#include "mesh_events.hh"
/* -------------------------------------------------------------------------- */
#include <set>
/* -------------------------------------------------------------------------- */

namespace akantu {
class Communicator;
class ElementSynchronizer;
class NodeSynchronizer;
} // namespace akantu

namespace akantu {

/* -------------------------------------------------------------------------- */
/* Mesh                                                                       */
/* -------------------------------------------------------------------------- */

/**
 * @class  Mesh this  contain the  coordinates of  the nodes  in  the Mesh.nodes
 * Array,  and the  connectivity. The  connectivity  are stored  in by  element
 * types.
 *
 * In order to loop on all element you have to loop on all types like this :
 * @code{.cpp}
 Mesh::type_iterator it = mesh.firstType(dim, ghost_type);
 Mesh::type_iterator end = mesh.lastType(dim, ghost_type);
 for(; it != end; ++it) {
   UInt nb_element  = mesh.getNbElement(*it);
   const Array<UInt> & conn = mesh.getConnectivity(*it);
   for(UInt e = 0; e < nb_element; ++e) {
     ...
   }
 }
 @endcode
*/
class Mesh : protected Memory,
             public EventHandlerManager<MeshEventHandler>,
             public GroupManager,
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

  /// constructor that use an existing nodes coordinates array, by knowing its
  /// ID
  // Mesh(UInt spatial_dimension, const ID & nodes_id, const ID & id,
  //      const MemoryID & memory_id = 0);

  /**
   * constructor that use an existing nodes coordinates
   * array, by getting the vector of coordinates
   */
  Mesh(UInt spatial_dimension, const std::shared_ptr<Array<Real>> & nodes,
       const ID & id = "mesh", const MemoryID & memory_id = 0);

  ~Mesh() override;

  /// @typedef ConnectivityTypeList list of the types present in a Mesh
  using ConnectivityTypeList = std::set<ElementType>;

  /// read the mesh from a file
  void read(const std::string & filename,
            const MeshIOType & mesh_io_type = _miot_auto);
  /// write the mesh to a file
  void write(const std::string & filename,
             const MeshIOType & mesh_io_type = _miot_auto);

private:
  /// initialize the connectivity to NULL and other stuff
  void init();

  /// function that computes the bounding box (fills xmin, xmax)
  void computeBoundingBox();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// patitionate the mesh among the processors involved in their computation
  virtual void distribute(Communicator & communicator);
  virtual void distribute();

  /// function to print the containt of the class
  void printself(std::ostream & stream, int indent = 0) const override;

  /// extract coordinates of nodes from an element
  template <typename T>
  inline void extractNodalValuesFromElement(const Array<T> & nodal_values,
                                            T * elemental_values,
                                            UInt * connectivity, UInt n_nodes,
                                            UInt nb_degree_of_freedom) const;

  // /// extract coordinates of nodes from a reversed element
  // inline void extractNodalCoordinatesFromPBCElement(Real * local_coords,
  //                                                   UInt * connectivity,
  //                                                   UInt n_nodes);

  /// convert a element to a linearized element
  inline UInt elementToLinearized(const Element & elem) const;

  /// convert a linearized element to an element
  inline Element linearizedToElement(UInt linearized_element) const;

  /// update the types offsets array for the conversions
  // inline void updateTypesOffsets(const GhostType & ghost_type);

  /// add a Array of connectivity for the type <type>.
  inline void addConnectivityType(const ElementType & type,
                                  const GhostType & ghost_type = _not_ghost);

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
  AKANTU_GET_MACRO_NOT_CONST(GlobalNodesIds, *nodes_global_ids, Array<UInt> &);

  /// get the global id of a node
  inline UInt getNodeGlobalId(UInt local_id) const;

  /// get the global id of a node
  inline UInt getNodeLocalId(UInt global_id) const;

  /// get the global number of nodes
  inline UInt getNbGlobalNodes() const;

  /// get the nodes type Array
  AKANTU_GET_MACRO(NodesType, nodes_type, const Array<NodeType> &);

protected:
  AKANTU_GET_MACRO_NOT_CONST(NodesType, nodes_type, Array<NodeType> &);

public:
  inline NodeType getNodeType(UInt local_id) const;

  /// say if a node is a pure ghost node
  inline bool isPureGhostNode(UInt n) const;

  /// say if a node is pur local or master node
  inline bool isLocalOrMasterNode(UInt n) const;

  inline bool isLocalNode(UInt n) const;
  inline bool isMasterNode(UInt n) const;
  inline bool isSlaveNode(UInt n) const;

  AKANTU_GET_MACRO(LowerBounds, lower_bounds, const Vector<Real> &);
  AKANTU_GET_MACRO(UpperBounds, upper_bounds, const Vector<Real> &);
  AKANTU_GET_MACRO(LocalLowerBounds, local_lower_bounds, const Vector<Real> &);
  AKANTU_GET_MACRO(LocalUpperBounds, local_upper_bounds, const Vector<Real> &);

  /// get the connectivity Array for a given type
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Connectivity, connectivities, UInt);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(Connectivity, connectivities, UInt);
  AKANTU_GET_MACRO(Connectivities, connectivities,
                   const ElementTypeMapArray<UInt> &);

  /// get the number of element of a type in the mesh
  inline UInt getNbElement(const ElementType & type,
                           const GhostType & ghost_type = _not_ghost) const;

  /// get the number of element for a given ghost_type and a given dimension
  inline UInt getNbElement(const UInt spatial_dimension = _all_dimensions,
                           const GhostType & ghost_type = _not_ghost,
                           const ElementKind & kind = _ek_not_defined) const;

  // /// get the connectivity list either for the elements or the ghost elements
  // inline const ConnectivityTypeList &
  // getConnectivityTypeList(const GhostType & ghost_type = _not_ghost) const;

  /// compute the barycenter of a given element
private:
  inline void getBarycenter(UInt element, const ElementType & type,
                            Real * barycenter,
                            GhostType ghost_type = _not_ghost) const;

public:
  inline void getBarycenter(const Element & element,
                            Vector<Real> & barycenter) const;

  /// get the element connected to a subelement
  const auto & getElementToSubelement() const;

  /// get the element connected to a subelement
  const auto &
  getElementToSubelement(const ElementType & el_type,
                         const GhostType & ghost_type = _not_ghost) const;

  /// get the element connected to a subelement
  auto & getElementToSubelement(const ElementType & el_type,
                                const GhostType & ghost_type = _not_ghost);

  /// get the subelement connected to an element
  const auto &
  getSubelementToElement(const ElementType & el_type,
                         const GhostType & ghost_type = _not_ghost) const;
  /// get the subelement connected to an element
  auto & getSubelementToElement(const ElementType & el_type,
                                const GhostType & ghost_type = _not_ghost);

  /// register a new ElementalTypeMap in the MeshData
  template <typename T>
  inline ElementTypeMapArray<T> & registerData(const std::string & data_name);

  template <typename T>
  inline bool hasData(const ID & data_name, const ElementType & el_type,
                      const GhostType & ghost_type = _not_ghost) const;

  /// get a name field associated to the mesh
  template <typename T>
  inline const Array<T> &
  getData(const ID & data_name, const ElementType & el_type,
          const GhostType & ghost_type = _not_ghost) const;

  /// get a name field associated to the mesh
  template <typename T>
  inline Array<T> & getData(const ID & data_name, const ElementType & el_type,
                            const GhostType & ghost_type = _not_ghost);

  /// get a name field associated to the mesh
  template <typename T>
  inline const ElementTypeMapArray<T> & getData(const ID & data_name) const;

  /// get a name field associated to the mesh
  template <typename T>
  inline ElementTypeMapArray<T> & getData(const ID & data_name);

  template <typename T>
  ElementTypeMap<UInt> getNbDataPerElem(ElementTypeMapArray<T> & array,
                                        const ElementKind & element_kind);

  template <typename T>
  dumper::Field * createFieldFromAttachedData(const std::string & field_id,
                                              const std::string & group_name,
                                              const ElementKind & element_kind);

  /// templated getter returning the pointer to data in MeshData (modifiable)
  template <typename T>
  inline Array<T> &
  getDataPointer(const std::string & data_name, const ElementType & el_type,
                 const GhostType & ghost_type = _not_ghost,
                 UInt nb_component = 1, bool size_to_nb_element = true,
                 bool resize_with_parent = false);

  /// Facets mesh accessor
  inline const Mesh & getMeshFacets() const;
  inline Mesh & getMeshFacets();

  /// Parent mesh accessor
  inline const Mesh & getMeshParent() const;

  inline bool isMeshFacets() const { return this->is_mesh_facets; }

  /// defines is the mesh is distributed or not
  inline bool isDistributed() const { return this->is_distributed; }

#ifndef SWIG
  /// return the dumper from a group and and a dumper name
  DumperIOHelper & getGroupDumper(const std::string & dumper_name,
                                  const std::string & group_name);
#endif
  /* ------------------------------------------------------------------------ */
  /* Wrappers on ElementClass functions                                       */
  /* ------------------------------------------------------------------------ */
public:
  /// get the number of nodes per element for a given element type
  static inline UInt getNbNodesPerElement(const ElementType & type);

  /// get the number of nodes per element for a given element type considered as
  /// a first order element
  static inline ElementType getP1ElementType(const ElementType & type);

  /// get the kind of the element type
  static inline ElementKind getKind(const ElementType & type);

  /// get spatial dimension of a type of element
  static inline UInt getSpatialDimension(const ElementType & type);

  /// get number of facets of a given element type
  static inline UInt getNbFacetsPerElement(const ElementType & type);

  /// get number of facets of a given element type
  static inline UInt getNbFacetsPerElement(const ElementType & type, UInt t);

  /// get local connectivity of a facet for a given facet type
  static inline auto getFacetLocalConnectivity(const ElementType & type,
                                               UInt t = 0);

  /// get connectivity of facets for a given element
  inline auto getFacetConnectivity(const Element & element, UInt t = 0) const;

  /// get the number of type of the surface element associated to a given
  /// element type
  static inline UInt getNbFacetTypes(const ElementType & type, UInt t = 0);

  /// get the type of the surface element associated to a given element
  static inline constexpr auto getFacetType(const ElementType & type,
                                            UInt t = 0);

  /// get all the type of the surface element associated to a given element
  static inline constexpr auto getAllFacetTypes(const ElementType & type);

  /// get the number of nodes in the given element list
  static inline UInt getNbNodesPerElementList(const Array<Element> & elements);

  /* ------------------------------------------------------------------------ */
  /* Element type Iterator                                                    */
  /* ------------------------------------------------------------------------ */
  using type_iterator = ElementTypeMapArray<UInt, ElementType>::type_iterator;
  using ElementTypesIteratorHelper =
      ElementTypeMapArray<UInt, ElementType>::ElementTypesIteratorHelper;

  template <typename... pack>
  ElementTypesIteratorHelper elementTypes(pack &&... _pack) const;

  inline type_iterator firstType(UInt dim = _all_dimensions,
                                 GhostType ghost_type = _not_ghost,
                                 ElementKind kind = _ek_regular) const {
    return connectivities.firstType(dim, ghost_type, kind);
  }

  inline type_iterator lastType(UInt dim = _all_dimensions,
                                GhostType ghost_type = _not_ghost,
                                ElementKind kind = _ek_regular) const {
    return connectivities.lastType(dim, ghost_type, kind);
  }

  AKANTU_GET_MACRO(ElementSynchronizer, *element_synchronizer,
                   const ElementSynchronizer &);
  AKANTU_GET_MACRO_NOT_CONST(ElementSynchronizer, *element_synchronizer,
                             ElementSynchronizer &);
  AKANTU_GET_MACRO(NodeSynchronizer, *node_synchronizer,
                   const NodeSynchronizer &);
  AKANTU_GET_MACRO_NOT_CONST(NodeSynchronizer, *node_synchronizer,
                             NodeSynchronizer &);

// AKANTU_GET_MACRO_NOT_CONST(Communicator, *communicator, StaticCommunicator
// &);
#ifndef SWIG
  AKANTU_GET_MACRO(Communicator, *communicator, const auto &);
  AKANTU_GET_MACRO_NOT_CONST(Communicator, *communicator, auto &);
#endif
  /* ------------------------------------------------------------------------ */
  /* Private methods for friends                                              */
  /* ------------------------------------------------------------------------ */
private:
  friend class MeshAccessor;

  AKANTU_GET_MACRO(NodesPointer, *nodes, Array<Real> &);

  /// get a pointer to the nodes_global_ids Array<UInt> and create it if
  /// necessary
  inline Array<UInt> & getNodesGlobalIdsPointer();

  /// get a pointer to the nodes_type Array<Int> and create it if necessary
  inline Array<NodeType> & getNodesTypePointer();

  /// get a pointer to the connectivity Array for the given type and create it
  /// if necessary
  inline Array<UInt> &
  getConnectivityPointer(const ElementType & type,
                         const GhostType & ghost_type = _not_ghost);

  /// get the ghost element counter
  inline Array<UInt> &
  getGhostsCounters(const ElementType & type,
                    const GhostType & ghost_type = _ghost) {
    AKANTU_DEBUG_ASSERT(ghost_type != _not_ghost,
                        "No ghost counter for _not_ghost elements");
    return ghosts_counters(type, ghost_type);
  }

  /// get a pointer to the element_to_subelement Array for the given type and
  /// create it if necessary
  inline Array<std::vector<Element>> &
  getElementToSubelementPointer(const ElementType & type,
                                const GhostType & ghost_type = _not_ghost);

  /// get a pointer to the subelement_to_element Array for the given type and
  /// create it if necessary
  inline Array<Element> &
  getSubelementToElementPointer(const ElementType & type,
                                const GhostType & ghost_type = _not_ghost);

  AKANTU_GET_MACRO_NOT_CONST(MeshData, mesh_data, MeshData &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// array of the nodes coordinates
  std::shared_ptr<Array<Real>> nodes;

  /// global node ids
  std::unique_ptr<Array<UInt>> nodes_global_ids;

  /// node type,  -3 pure ghost, -2  master for the  node, -1 normal node,  i in
  /// [0-N] slave node and master is proc i
  Array<NodeType> nodes_type;

  /// global number of nodes;
  UInt nb_global_nodes{0};

  /// all class of elements present in this mesh (for heterogenous meshes)
  ElementTypeMapArray<UInt> connectivities;

  /// count the references on ghost elements
  ElementTypeMapArray<UInt> ghosts_counters;

  /// map to normals for all class of elements present in this mesh
  ElementTypeMapArray<Real> normals;

  /// list of all existing types in the mesh
  // ConnectivityTypeList type_set;

  /// the spatial dimension of this mesh
  UInt spatial_dimension{0};

  /// types offsets
  // Array<UInt> types_offsets;

  // /// list of all existing types in the mesh
  // ConnectivityTypeList ghost_type_set;
  // /// ghost types offsets
  // Array<UInt> ghost_types_offsets;

  /// min of coordinates
  Vector<Real> lower_bounds;
  /// max of coordinates
  Vector<Real> upper_bounds;
  /// size covered by the mesh on each direction
  Vector<Real> size;

  /// local min of coordinates
  Vector<Real> local_lower_bounds;
  /// local max of coordinates
  Vector<Real> local_upper_bounds;

  /// Extra data loaded from the mesh file
  MeshData mesh_data;

  /// facets' mesh
  std::unique_ptr<Mesh> mesh_facets;

  /// parent mesh (this is set for mesh_facets meshes)
  const Mesh * mesh_parent{nullptr};

  /// defines if current mesh is mesh_facets or not
  bool is_mesh_facets{false};

  /// defines if the mesh is centralized or distributed
  bool is_distributed{false};

  /// Communicator on which mesh is distributed
  Communicator * communicator;

  /// Element synchronizer
  std::unique_ptr<ElementSynchronizer> element_synchronizer;

  /// Node synchronizer
  std::unique_ptr<NodeSynchronizer> node_synchronizer;

  using NodesToElements = std::vector<std::unique_ptr<std::set<Element>>>;

  /// This info is stored to simplify the dynamic changes
  NodesToElements nodes_to_elements;
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
#include "mesh_inline_impl.cc"

#endif /* __AKANTU_MESH_HH__ */
