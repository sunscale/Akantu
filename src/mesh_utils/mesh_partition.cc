/**
 * @file   mesh_partition.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Aug 17 2010
 * @date last modification: Mon Jul 21 2014
 *
 * @brief  implementation of common part of all partitioner
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "mesh_partition.hh"
#include "mesh_utils.hh"
#include "aka_types.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
MeshPartition::MeshPartition(const Mesh & mesh, UInt spatial_dimension,
			     const ID & id,
			     const MemoryID & memory_id) :
  Memory(id, memory_id),
  mesh(mesh), spatial_dimension(spatial_dimension),
  partitions             ("partition"             , id, memory_id),
  ghost_partitions       ("ghost_partition"       , id, memory_id),
  ghost_partitions_offset("ghost_partition_offset", id, memory_id),
  saved_connectivity("saved_connectivity", id, memory_id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
MeshPartition::~MeshPartition() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * conversion in c++ of the GENDUALMETIS (mesh.c) function wrote by George in
 * Metis (University of Minnesota)
 */
void MeshPartition::buildDualGraph(Array<Int> & dxadj, Array<Int> & dadjncy,
				   Array<Int> & edge_loads,
				   const EdgeLoadFunctor & edge_load_func) {
  AKANTU_DEBUG_IN();

  // tweak mesh;
  const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
  UInt nb_types = type_list.size();
  UInt nb_good_types = 0;

  UInt nb_nodes_per_element_p1[nb_types];

  UInt magic_number[nb_types];

  //  UInt * conn_val[nb_types];
  UInt nb_element[nb_types];

  Array<UInt> * conn[nb_types];

  Array<Element> lin_to_element;

  Element elem;
  elem.ghost_type = _not_ghost;

  const_cast<Mesh &>(mesh).updateTypesOffsets(_not_ghost);
  const_cast<Mesh &>(mesh).updateTypesOffsets(_ghost);

  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    ElementType type = *it;
    if(Mesh::getSpatialDimension(type) != mesh.getSpatialDimension()) continue;
    elem.type = type;

    ElementType type_p1 = Mesh::getP1ElementType(type);

    nb_nodes_per_element_p1[nb_good_types] = Mesh::getNbNodesPerElement(type_p1);
    nb_element[nb_good_types]              = mesh.getConnectivity(type, _not_ghost).getSize();
    magic_number[nb_good_types]            =
      Mesh::getNbNodesPerElement(Mesh::getFacetType(type_p1));

    conn[nb_good_types] = &const_cast<Array<UInt> &>(mesh.getConnectivity(type, _not_ghost));

    for (UInt i = 0; i < nb_element[nb_good_types]; ++i) {
      elem.element = i;
      lin_to_element.push_back(elem);
    }

    nb_good_types++;
  }

  CSR<UInt> node_to_elem;

  MeshUtils::buildNode2Elements(mesh, node_to_elem);

  UInt nb_total_element = 0;
  UInt nb_total_node_element = 0;
  for (UInt t = 0; t < nb_good_types; ++t) {
    nb_total_element += nb_element[t];
    nb_total_node_element += nb_element[t]*nb_nodes_per_element_p1[t];
  }

  dxadj.resize(nb_total_element + 1);

  /// initialize the dxadj array
  for (UInt t = 0, linerized_el = 0; t < nb_good_types; ++t)
    for (UInt el = 0; el < nb_element[t]; ++el, ++linerized_el)
      dxadj(linerized_el) = nb_nodes_per_element_p1[t];

  /// convert the dxadj_val array in a csr one
  for (UInt i = 1; i < nb_total_element; ++i) dxadj(i) += dxadj(i-1);
  for (UInt i = nb_total_element; i > 0; --i) dxadj(i)  = dxadj(i-1);
  dxadj(0) = 0;

  dadjncy.resize(2*dxadj(nb_total_element));
  edge_loads.resize(2*dxadj(nb_total_element));

  /// weight map to determine adjacency
  unordered_map<UInt, UInt>::type weight_map;

  for (UInt t = 0, linerized_el = 0; t < nb_good_types; ++t) {
    for (UInt el = 0; el < nb_element[t]; ++el, ++linerized_el) {
      /// fill the weight map
      for (UInt n = 0; n < nb_nodes_per_element_p1[t]; ++n) {
	UInt node = (*conn[t])(el, n);
	CSR<UInt>::iterator k;
	for (k = node_to_elem.rbegin(node); k != node_to_elem.rend(node); --k) {
	  UInt current_el = *k;
	  if(current_el <= linerized_el) break;

	  unordered_map<UInt, UInt>::type::iterator it_w;
	  it_w = weight_map.find(current_el);

	  if(it_w == weight_map.end()) {
	    weight_map[current_el] = 1;
	  } else {
	    it_w->second++;
	  }
	}
      }
      /// each element with a weight of the size of a facet are adjacent
      unordered_map<UInt, UInt>::type::iterator it_w;
      for(it_w = weight_map.begin(); it_w != weight_map.end(); ++it_w) {
	if(it_w->second == magic_number[t]) {
	  UInt adjacent_el = it_w->first;

#if defined(AKANTU_COHESIVE_ELEMENT)
	  /// Patch in order to prevent neighboring cohesive elements
	  /// from detecting each other
	  ElementKind linearized_el_kind = mesh.linearizedToElement(linerized_el).kind;
	  ElementKind adjacent_el_kind = mesh.linearizedToElement(adjacent_el).kind;

	  if (linearized_el_kind == adjacent_el_kind &&
	      linearized_el_kind == _ek_cohesive) continue;
#endif

	  UInt index_adj = dxadj(adjacent_el )++;
	  UInt index_lin = dxadj(linerized_el)++;

	  dadjncy(index_lin) = adjacent_el;
	  dadjncy(index_adj) = linerized_el;
	}
      }

      weight_map.clear();
    }
  }

  Int k_start = 0;
  for (UInt t = 0, linerized_el = 0, j = 0; t < nb_good_types; ++t)
    for (UInt el = 0; el < nb_element[t]; ++el, ++linerized_el) {
      for (Int k = k_start; k < dxadj(linerized_el); ++k, ++j)
	dadjncy(j) = dadjncy(k);
      dxadj(linerized_el) = j;
      k_start += nb_nodes_per_element_p1[t];
    }

  for (UInt i = nb_total_element; i > 0; --i) dxadj(i) = dxadj(i - 1);
  dxadj(0) = 0;

  UInt adj = 0;
  for (UInt i = 0; i < nb_total_element; ++i) {
    UInt nb_adj = dxadj(i + 1) - dxadj(i);
    for (UInt j = 0; j < nb_adj; ++j, ++adj) {
      Int el_adj_id = dadjncy(dxadj(i) + j);
      Element el     = lin_to_element(i);
      Element el_adj = lin_to_element(el_adj_id);

      Int load = edge_load_func(el, el_adj);
      edge_loads(adj) = load;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshPartition::fillPartitionInformation(const Mesh & mesh,
					     const Int * linearized_partitions) {
  AKANTU_DEBUG_IN();

  CSR<UInt> node_to_elem;

  MeshUtils::buildNode2Elements(mesh, node_to_elem);

  Mesh::type_iterator it  = mesh.firstType(spatial_dimension,
					   _not_ghost,
					   _ek_not_defined);
  Mesh::type_iterator end = mesh.lastType(spatial_dimension,
					  _not_ghost,
					  _ek_not_defined);

  UInt linearized_el = 0;
  for(; it != end; ++it) {
    ElementType type = *it;

    UInt nb_element = mesh.getNbElement(type);
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

    partitions             .alloc(nb_element,     1, type, _not_ghost);
    CSR<UInt> & ghost_part_csr = ghost_partitions_csr(type, _not_ghost);
    ghost_part_csr.resizeRows(nb_element);

    ghost_partitions_offset.alloc(nb_element + 1, 1, type, _ghost);
    ghost_partitions       .alloc(0,              1, type, _ghost);

    const Array<UInt> & connectivity = mesh.getConnectivity(type, _not_ghost);

    for (UInt el = 0; el < nb_element; ++el, ++linearized_el) {
      UInt part = linearized_partitions[linearized_el];

      partitions(type, _not_ghost)(el) = part;
      std::list<UInt> list_adj_part;
      for (UInt n = 0; n < nb_nodes_per_element; ++n) {
	UInt node = connectivity.storage()[el * nb_nodes_per_element + n];
	CSR<UInt>::iterator ne;
	for (ne = node_to_elem.begin(node); ne != node_to_elem.end(node); ++ne) {
	  UInt adj_el = *ne;
	  UInt adj_part = linearized_partitions[adj_el];
	  if(part != adj_part) {
	    list_adj_part.push_back(adj_part);
	  }
	}
      }

      list_adj_part.sort();
      list_adj_part.unique();

      for(std::list<UInt>::iterator adj_it = list_adj_part.begin();
	  adj_it != list_adj_part.end();
	  ++adj_it) {
	ghost_part_csr.getRows().push_back(*adj_it);
	ghost_part_csr.rowOffset(el)++;

	ghost_partitions(type, _ghost).push_back(*adj_it);
	ghost_partitions_offset(type, _ghost)(el)++;
      }
    }

    ghost_part_csr.countToCSR();

    /// convert the ghost_partitions_offset array in an offset array
    Array<UInt> & ghost_partitions_offset_ptr = ghost_partitions_offset(type, _ghost);
    for (UInt i = 1; i < nb_element; ++i)
      ghost_partitions_offset_ptr(i) += ghost_partitions_offset_ptr(i-1);
    for (UInt i = nb_element; i > 0; --i)
      ghost_partitions_offset_ptr(i)  = ghost_partitions_offset_ptr(i-1);
    ghost_partitions_offset_ptr(0) = 0;
  }

  // All Facets
  for(Int sp = spatial_dimension - 1; sp >= 0; --sp) {
    Mesh::type_iterator fit  = mesh.firstType(sp,
					      _not_ghost,
					      _ek_not_defined);
    Mesh::type_iterator fend = mesh.lastType(sp,
					     _not_ghost,
					     _ek_not_defined);

    for(; fit != fend; ++fit) {
      ElementType type = *fit;

      UInt nb_element = mesh.getNbElement(type);

      partitions             .alloc(nb_element,     1, type, _not_ghost);
      AKANTU_DEBUG_INFO("Allocating partitions for " << type);
      CSR<UInt> & ghost_part_csr = ghost_partitions_csr(type, _not_ghost);
      ghost_part_csr.resizeRows(nb_element);

      ghost_partitions_offset.alloc(nb_element + 1, 1, type, _ghost);
      ghost_partitions       .alloc(0,              1, type, _ghost);
      AKANTU_DEBUG_INFO("Allocating ghost_partitions for " << type);
      const Array< std::vector<Element> > & elem_to_subelem = mesh.getElementToSubelement(type, _not_ghost);
      for(UInt i(0); i < mesh.getNbElement(type, _not_ghost); ++i) { // Facet loop

	const std::vector<Element> & adjacent_elems = elem_to_subelem(i);
        if(!adjacent_elems.empty()) {
          Element min_elem;
          UInt min_part(std::numeric_limits<UInt>::max());
          std::set<UInt> adjacent_parts;

          for(UInt j(0); j < adjacent_elems.size(); ++j) {
            UInt adjacent_elem_id = adjacent_elems[j].element;
            UInt adjacent_elem_part = partitions(adjacent_elems[j].type, adjacent_elems[j].ghost_type)(adjacent_elem_id);
            if(adjacent_elem_part < min_part) {
              min_part = adjacent_elem_part;
              min_elem = adjacent_elems[j];
            }
            adjacent_parts.insert(adjacent_elem_part);
          }
          partitions(type, _not_ghost)(i) = min_part;

          CSR<UInt>::iterator git = ghost_partitions_csr(min_elem.type, _not_ghost).begin(min_elem.element);
          CSR<UInt>::iterator gend = ghost_partitions_csr(min_elem.type, _not_ghost).end(min_elem.element);
          for(; git != gend; ++git) {
            adjacent_parts.insert(*git);
          }
          adjacent_parts.erase(min_part);
          std::set<UInt>::const_iterator pit = adjacent_parts.begin();
          std::set<UInt>::const_iterator pend = adjacent_parts.end();
          for(; pit != pend; ++pit) {
            ghost_part_csr.getRows().push_back(*pit);
            ghost_part_csr.rowOffset(i)++;
            ghost_partitions(type, _ghost).push_back(*pit);
          }

          ghost_partitions_offset(type, _ghost)(i+1) = ghost_partitions_offset(type, _ghost)(i+1)
            + adjacent_elems.size();
        } else {
          partitions(type, _not_ghost)(i) = 0;
        }
      }
      ghost_part_csr.countToCSR();
    }
  }
  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void MeshPartition::tweakConnectivity(const Array<UInt> & pairs) {
  AKANTU_DEBUG_IN();

  if(pairs.getSize() == 0) return;

  Mesh::type_iterator it  = mesh.firstType(spatial_dimension,
					   _not_ghost,
					   _ek_not_defined);
  Mesh::type_iterator end = mesh.lastType(spatial_dimension,
					  _not_ghost,
					  _ek_not_defined);

  for(; it != end; ++it) {
    ElementType type = *it;

    Array<UInt> & conn = const_cast<Array<UInt> &>(mesh.getConnectivity(type, _not_ghost));
    UInt nb_nodes_per_element = conn.getNbComponent();
    UInt nb_element           = conn.getSize();

    Array<UInt> & saved_conn = saved_connectivity.alloc(nb_element, nb_nodes_per_element, type, _not_ghost);
    saved_conn.copy(conn);

    for (UInt i = 0; i < pairs.getSize(); ++i) {
      for (UInt el = 0; el < nb_element; ++el) {
	for (UInt n = 0; n < nb_nodes_per_element; ++n) {
	  if(pairs(i, 1) == conn(el, n))
	    conn(el, n) = pairs(i, 0);
	}
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshPartition::restoreConnectivity() {
  AKANTU_DEBUG_IN();

  ElementTypeMapArray<UInt>::type_iterator it   = saved_connectivity.firstType(spatial_dimension,
								       _not_ghost,
								       _ek_not_defined);
  ElementTypeMapArray<UInt>::type_iterator end = saved_connectivity.lastType(spatial_dimension,
								     _not_ghost,
								     _ek_not_defined);
  for(; it != end; ++it) {
    ElementType type = *it;

    Array<UInt> & conn = const_cast<Array<UInt> &>(mesh.getConnectivity(type, _not_ghost));
    Array<UInt> & saved_conn = saved_connectivity(type, _not_ghost);
    conn.copy(saved_conn);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__
