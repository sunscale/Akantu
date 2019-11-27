/**
 * @file   master_element_info_per_processor.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Mar 16 2016
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Helper class to distribute a mesh
 *
 * @section LICENSE
 *
 * Copyright (©) 2016-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_iterators.hh"
#include "communicator.hh"
#include "element_group.hh"
#include "element_info_per_processor.hh"
#include "element_synchronizer.hh"
#include "mesh_iterators.hh"
#include "mesh_utils.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
#include <iostream>
#include <map>
#include <tuple>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
MasterElementInfoPerProc::MasterElementInfoPerProc(
    ElementSynchronizer & synchronizer, UInt message_cnt, UInt root,
    ElementType type, const MeshPartition & partition)
    : ElementInfoPerProc(synchronizer, message_cnt, root, type),
      partition(partition), all_nb_local_element(nb_proc, 0),
      all_nb_ghost_element(nb_proc, 0), all_nb_element_to_send(nb_proc, 0) {
  Vector<UInt> size(5);
  size(0) = (UInt)type;

  if (type != _not_defined) {
    nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    nb_element = mesh.getNbElement(type);
    const auto & partition_num =
        this->partition.getPartition(this->type, _not_ghost);
    const auto & ghost_partition =
        this->partition.getGhostPartitionCSR()(this->type, _not_ghost);

    for (UInt el = 0; el < nb_element; ++el) {
      this->all_nb_local_element[partition_num(el)]++;
      for (auto part = ghost_partition.begin(el);
           part != ghost_partition.end(el); ++part) {
        this->all_nb_ghost_element[*part]++;
      }
      this->all_nb_element_to_send[partition_num(el)] +=
          ghost_partition.getNbCols(el) + 1;
    }

    /// tag info
    auto && tag_names = this->mesh.getTagNames(type);
    this->nb_tags = tag_names.size();

    size(4) = nb_tags;

    for (UInt p = 0; p < nb_proc; ++p) {
      if (p != root) {
        size(1) = this->all_nb_local_element[p];
        size(2) = this->all_nb_ghost_element[p];
        size(3) = this->all_nb_element_to_send[p];
        AKANTU_DEBUG_INFO(
            "Sending connectivities informations to proc "
            << p << " TAG("
            << Tag::genTag(this->rank, this->message_count, Tag::_SIZES)
            << ")");
        comm.send(size, p,
                  Tag::genTag(this->rank, this->message_count, Tag::_SIZES));
      } else {
        this->nb_local_element = this->all_nb_local_element[p];
        this->nb_ghost_element = this->all_nb_ghost_element[p];
      }
    }
  } else {
    for (UInt p = 0; p < this->nb_proc; ++p) {
      if (p != this->root) {
        AKANTU_DEBUG_INFO(
            "Sending empty connectivities informations to proc "
            << p << " TAG("
            << Tag::genTag(this->rank, this->message_count, Tag::_SIZES)
            << ")");
        comm.send(size, p,
                  Tag::genTag(this->rank, this->message_count, Tag::_SIZES));
      }
    }
  }
}

/* ------------------------------------------------------------------------ */
void MasterElementInfoPerProc::synchronizeConnectivities() {
  const auto & partition_num =
      this->partition.getPartition(this->type, _not_ghost);
  const auto & ghost_partition =
      this->partition.getGhostPartitionCSR()(this->type, _not_ghost);

  std::vector<Array<UInt>> buffers(this->nb_proc);

  const auto & connectivities =
      this->mesh.getConnectivity(this->type, _not_ghost);

  /// copying the local connectivity
  for (auto && part_conn :
       zip(partition_num,
           make_view(connectivities, this->nb_nodes_per_element))) {
    auto && part = std::get<0>(part_conn);
    auto && conn = std::get<1>(part_conn);
    for (UInt i = 0; i < conn.size(); ++i) {
      buffers[part].push_back(conn[i]);
    }
  }

  /// copying the connectivity of ghost element
  for (auto && tuple :
       enumerate(make_view(connectivities, this->nb_nodes_per_element))) {
    auto && el = std::get<0>(tuple);
    auto && conn = std::get<1>(tuple);
    for (auto part = ghost_partition.begin(el); part != ghost_partition.end(el);
         ++part) {
      UInt proc = *part;
      for (UInt i = 0; i < conn.size(); ++i) {
        buffers[proc].push_back(conn[i]);
      }
    }
  }

#ifndef AKANTU_NDEBUG
  for (auto p : arange(this->nb_proc)) {
    UInt size = this->nb_nodes_per_element *
                (this->all_nb_local_element[p] + this->all_nb_ghost_element[p]);
    AKANTU_DEBUG_ASSERT(
        buffers[p].size() == size,
        "The connectivity data packed in the buffer are not correct");
  }
#endif

  /// send all connectivity and ghost information to all processors
  std::vector<CommunicationRequest> requests;
  for (auto p : arange(this->nb_proc)) {
    if (p == this->root)
      continue;
    auto && tag =
        Tag::genTag(this->rank, this->message_count, Tag::_CONNECTIVITY);
    AKANTU_DEBUG_INFO("Sending connectivities to proc " << p << " TAG(" << tag
                                                        << ")");
    requests.push_back(comm.asyncSend(buffers[p], p, tag));
  }

  Array<UInt> & old_nodes = this->getNodesGlobalIds();

  /// create the renumbered connectivity
  AKANTU_DEBUG_INFO("Renumbering local connectivities");
  MeshUtils::renumberMeshNodes(mesh, buffers[root], all_nb_local_element[root],
                               all_nb_ghost_element[root], type, old_nodes);

  comm.waitAll(requests);
  comm.freeCommunicationRequest(requests);
}

/* ------------------------------------------------------------------------ */
void MasterElementInfoPerProc::synchronizePartitions() {
  const auto & partition_num =
      this->partition.getPartition(this->type, _not_ghost);
  const auto & ghost_partition =
      this->partition.getGhostPartitionCSR()(this->type, _not_ghost);

  std::vector<Array<UInt>> buffers(this->partition.getNbPartition());

  /// splitting the partition information to send them to processors
  Vector<UInt> count_by_proc(nb_proc, 0);
  for (UInt el = 0; el < nb_element; ++el) {
    UInt proc = partition_num(el);
    buffers[proc].push_back(ghost_partition.getNbCols(el));

    UInt i(0);
    for (auto part = ghost_partition.begin(el); part != ghost_partition.end(el);
         ++part, ++i) {
      buffers[proc].push_back(*part);
    }
  }

  for (UInt el = 0; el < nb_element; ++el) {
    UInt i(0);
    for (auto part = ghost_partition.begin(el); part != ghost_partition.end(el);
         ++part, ++i) {
      buffers[*part].push_back(partition_num(el));
    }
  }

#ifndef AKANTU_NDEBUG
  for (UInt p = 0; p < this->nb_proc; ++p) {
    AKANTU_DEBUG_ASSERT(buffers[p].size() == (this->all_nb_ghost_element[p] +
                                              this->all_nb_element_to_send[p]),
                        "Data stored in the buffer are most probably wrong");
  }
#endif

  std::vector<CommunicationRequest> requests;
  /// last data to compute the communication scheme
  for (UInt p = 0; p < this->nb_proc; ++p) {
    if (p == this->root)
      continue;

    auto && tag =
        Tag::genTag(this->rank, this->message_count, Tag::_PARTITIONS);
    AKANTU_DEBUG_INFO("Sending partition informations to proc " << p << " TAG("
                                                                << tag << ")");
    requests.push_back(comm.asyncSend(buffers[p], p, tag));
  }

  if (Mesh::getSpatialDimension(this->type) ==
      this->mesh.getSpatialDimension()) {
    AKANTU_DEBUG_INFO("Creating communications scheme");
    this->fillCommunicationScheme(buffers[this->rank]);
  }

  comm.waitAll(requests);
  comm.freeCommunicationRequest(requests);
}

/* -------------------------------------------------------------------------- */
void MasterElementInfoPerProc::synchronizeTags() {
  AKANTU_DEBUG_IN();

  if (this->nb_tags == 0) {
    AKANTU_DEBUG_OUT();
    return;
  }

  /// tag info
  auto tag_names = mesh.getTagNames(type);

  // Make sure the tags are sorted (or at least not in random order),
  // because they come from a map !!
  std::sort(tag_names.begin(), tag_names.end());

  // Sending information about the tags in mesh_data: name, data type and
  // number of components of the underlying array associated to the current
  // type
  DynamicCommunicationBuffer mesh_data_sizes_buffer;
  for (auto && tag_name : tag_names) {
    mesh_data_sizes_buffer << tag_name;
    mesh_data_sizes_buffer << mesh.getTypeCode(tag_name);
    mesh_data_sizes_buffer << mesh.getNbComponent(tag_name, type);
  }

  AKANTU_DEBUG_INFO(
      "Broadcasting the size of the information about the mesh data tags: ("
      << mesh_data_sizes_buffer.size() << ").");
  AKANTU_DEBUG_INFO(
      "Broadcasting the information about the mesh data tags, addr "
      << (void *)mesh_data_sizes_buffer.storage());

  comm.broadcast(mesh_data_sizes_buffer, root);

  if (mesh_data_sizes_buffer.size() == 0)
    return;

  // Sending the actual data to each processor
  std::vector<DynamicCommunicationBuffer> buffers(nb_proc);
  // Loop over each tag for the current type
  for (auto && tag_name : tag_names) {
    // Type code of the current tag (i.e. the tag named *names_it)
    this->fillTagBuffer(buffers, tag_name);
  }

  std::vector<CommunicationRequest> requests;
  for (UInt p = 0; p < nb_proc; ++p) {
    if (p == root)
      continue;

    auto && tag = Tag::genTag(this->rank, this->message_count, Tag::_MESH_DATA);
    AKANTU_DEBUG_INFO("Sending " << buffers[p].size()
                                 << " bytes of mesh data to proc " << p
                                 << " TAG(" << tag << ")");

    requests.push_back(comm.asyncSend(buffers[p], p, tag));
  }

  // Loop over each tag for the current type
  for (auto && tag_name : tag_names) {
    // Reinitializing the mesh data on the master
    this->fillMeshData(buffers[root], tag_name, mesh.getTypeCode(tag_name),
                       mesh.getNbComponent(tag_name, type));
  }

  comm.waitAll(requests);
  comm.freeCommunicationRequest(requests);
  requests.clear();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <typename T>
void MasterElementInfoPerProc::fillTagBufferTemplated(
    std::vector<DynamicCommunicationBuffer> & buffers,
    const std::string & tag_name) {
  const auto & data = mesh.getElementalDataArray<T>(tag_name, type);
  const auto & partition_num =
      this->partition.getPartition(this->type, _not_ghost);
  const auto & ghost_partition =
      this->partition.getGhostPartitionCSR()(this->type, _not_ghost);

  // Not possible to use the iterator because it potentially triggers the
  // creation of complex
  // type templates (such as akantu::Vector< std::vector<Element> > which don't
  // implement the right interface
  // (e.g. operator<< in that case).
  // typename Array<T>::template const_iterator< Vector<T> > data_it  =
  // data.begin(data.getNbComponent());
  // typename Array<T>::template const_iterator< Vector<T> > data_end =
  // data.end(data.getNbComponent());

  const T * data_it = data.storage();
  const T * data_end = data.storage() + data.size() * data.getNbComponent();
  const UInt * part = partition_num.storage();

  /// copying the data, element by element
  for (; data_it != data_end; ++part) {
    for (UInt j(0); j < data.getNbComponent(); ++j, ++data_it) {
      buffers[*part] << *data_it;
    }
  }

  data_it = data.storage();
  /// copying the data for the ghost element
  for (UInt el(0); data_it != data_end;
       data_it += data.getNbComponent(), ++el) {
    auto it = ghost_partition.begin(el);
    auto end = ghost_partition.end(el);
    for (; it != end; ++it) {
      UInt proc = *it;
      for (UInt j(0); j < data.getNbComponent(); ++j) {
        buffers[proc] << data_it[j];
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void MasterElementInfoPerProc::fillTagBuffer(
    std::vector<DynamicCommunicationBuffer> & buffers,
    const std::string & tag_name) {
#define AKANTU_DISTRIBUTED_SYNHRONIZER_TAG_DATA(r, extra_param, elem)          \
  case MeshDataTypeCode::BOOST_PP_TUPLE_ELEM(2, 0, elem): {                    \
    this->fillTagBufferTemplated<BOOST_PP_TUPLE_ELEM(2, 1, elem)>(buffers,     \
                                                                  tag_name);   \
    break;                                                                     \
  }

  MeshDataTypeCode data_type_code = mesh.getTypeCode(tag_name);
  switch (data_type_code) {
    BOOST_PP_SEQ_FOR_EACH(AKANTU_DISTRIBUTED_SYNHRONIZER_TAG_DATA, ,
                          AKANTU_MESH_DATA_TYPES)
  default:
    AKANTU_ERROR("Could not obtain the type of tag" << tag_name << "!");
    break;
  }
#undef AKANTU_DISTRIBUTED_SYNHRONIZER_TAG_DATA
}

/* -------------------------------------------------------------------------- */
void MasterElementInfoPerProc::synchronizeGroups() {
  AKANTU_DEBUG_IN();

  std::vector<DynamicCommunicationBuffer> buffers(nb_proc);

  using ElementToGroup = std::vector<std::vector<std::string>>;
  ElementToGroup element_to_group(nb_element);

  for (auto & eg : mesh.iterateElementGroups()) {
    const auto & name = eg.getName();

    for (const auto & element : eg.getElements(type, _not_ghost)) {
      element_to_group[element].push_back(name);
    }

    auto eit = eg.begin(type, _not_ghost);
    if (eit != eg.end(type, _not_ghost))
      const_cast<Array<UInt> &>(eg.getElements(type)).empty();
  }

  const auto & partition_num =
      this->partition.getPartition(this->type, _not_ghost);
  const auto & ghost_partition =
      this->partition.getGhostPartitionCSR()(this->type, _not_ghost);

  /// copying the data, element by element
  for (auto && pair : zip(partition_num, element_to_group)) {
    buffers[std::get<0>(pair)] << std::get<1>(pair);
  }

  /// copying the data for the ghost element
  for (auto && pair : enumerate(element_to_group)) {
    auto && el = std::get<0>(pair);
    auto it = ghost_partition.begin(el);
    auto end = ghost_partition.end(el);
    for (; it != end; ++it) {
      UInt proc = *it;
      buffers[proc] << std::get<1>(pair);
    }
  }

  std::vector<CommunicationRequest> requests;
  for (UInt p = 0; p < this->nb_proc; ++p) {
    if (p == this->rank)
      continue;

    auto && tag = Tag::genTag(this->rank, p, Tag::_ELEMENT_GROUP);
    AKANTU_DEBUG_INFO("Sending element groups to proc " << p << " TAG(" << tag
                                                        << ")");
    requests.push_back(comm.asyncSend(buffers[p], p, tag));
  }

  this->fillElementGroupsFromBuffer(buffers[this->rank]);

  comm.waitAll(requests);
  comm.freeCommunicationRequest(requests);
  requests.clear();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
