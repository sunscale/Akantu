/**
 * @file   master_element_info_per_processor.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Mar 11 14:57:13 2016
 *
 * @brief
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_zip.hh"
#include "element_group.hh"
#include "element_info_per_processor.hh"
#include "element_synchronizer.hh"
#include "mesh_utils.hh"
#include "static_communicator.hh"
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
    const Array<UInt> & partition_num =
        this->partition.getPartition(this->type, _not_ghost);
    const CSR<UInt> & ghost_partition =
        this->partition.getGhostPartitionCSR()(this->type, _not_ghost);

    for (UInt el = 0; el < nb_element; ++el) {
      this->all_nb_local_element[partition_num(el)]++;
      for (CSR<UInt>::const_iterator part = ghost_partition.begin(el);
           part != ghost_partition.end(el); ++part) {
        this->all_nb_ghost_element[*part]++;
      }
      this->all_nb_element_to_send[partition_num(el)] +=
          ghost_partition.getNbCols(el) + 1;
    }

    /// tag info
    std::vector<std::string> tag_names;
    this->getMeshData().getTagNames(tag_names, type);
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
  const Array<UInt> & partition_num =
      this->partition.getPartition(this->type, _not_ghost);
  const CSR<UInt> & ghost_partition =
      this->partition.getGhostPartitionCSR()(this->type, _not_ghost);

  std::vector<Array<UInt>> buffers(this->nb_proc);

  auto conn_it = this->mesh.getConnectivity(this->type, _not_ghost)
                     .begin(this->nb_nodes_per_element);
  auto conn_end = this->mesh.getConnectivity(this->type, _not_ghost)
                      .end(this->nb_nodes_per_element);

  /// copying the local connectivity
  auto part_it = partition_num.begin();
  for (; conn_it != conn_end; ++conn_it, ++part_it) {
    const auto & conn = *conn_it;
    for (UInt i = 0; i < conn.size(); ++i) {
      buffers[*part_it].push_back(conn[i]);
    }
  }

  /// copying the connectivity of ghost element
  conn_it = this->mesh.getConnectivity(this->type, _not_ghost)
                .begin(this->nb_nodes_per_element);
  for (UInt el = 0; conn_it != conn_end; ++el, ++conn_it) {
    for (auto part = ghost_partition.begin(el); part != ghost_partition.end(el);
         ++part) {
      UInt proc = *part;
      const Vector<UInt> & conn = *conn_it;
      for (UInt i = 0; i < conn.size(); ++i) {
        buffers[proc].push_back(conn[i]);
      }
    }
  }

#ifndef AKANTU_NDEBUG
  for (UInt p = 0; p < this->nb_proc; ++p) {
    UInt size = this->nb_nodes_per_element *
                (this->all_nb_local_element[p] + this->all_nb_ghost_element[p]);
    AKANTU_DEBUG_ASSERT(
        buffers[p].getSize() == size,
        "The connectivity data packed in the buffer are not correct");
  }
#endif

  /// send all connectivity and ghost information to all processors
  std::vector<CommunicationRequest> requests;
  for (UInt p = 0; p < this->nb_proc; ++p) {
    if (p != this->root) {
      AKANTU_DEBUG_INFO(
          "Sending connectivities to proc "
          << p << " TAG("
          << Tag::genTag(this->rank, this->message_count, Tag::_CONNECTIVITY)
          << ")");
      requests.push_back(comm.asyncSend(
          buffers[p], p,
          Tag::genTag(this->rank, this->message_count, Tag::_CONNECTIVITY)));
    }
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
  const Array<UInt> & partition_num =
      this->partition.getPartition(this->type, _not_ghost);
  const CSR<UInt> & ghost_partition =
      this->partition.getGhostPartitionCSR()(this->type, _not_ghost);

  std::vector<Array<UInt>> buffers(this->partition.getNbPartition());

  /// splitting the partition information to send them to processors
  Vector<UInt> count_by_proc(nb_proc, 0);
  for (UInt el = 0; el < nb_element; ++el) {
    UInt proc = partition_num(el);
    buffers[proc].push_back(ghost_partition.getNbCols(el));

    UInt i(0);
    for (CSR<UInt>::const_iterator part = ghost_partition.begin(el);
         part != ghost_partition.end(el); ++part, ++i) {
      buffers[proc].push_back(*part);
    }
  }

  for (UInt el = 0; el < nb_element; ++el) {
    UInt i(0);
    for (CSR<UInt>::const_iterator part = ghost_partition.begin(el);
         part != ghost_partition.end(el); ++part, ++i) {
      buffers[*part].push_back(partition_num(el));
    }
  }

#ifndef AKANTU_NDEBUG
  for (UInt p = 0; p < this->nb_proc; ++p) {
    AKANTU_DEBUG_ASSERT(
        buffers[p].getSize() ==
            (this->all_nb_ghost_element[p] + this->all_nb_element_to_send[p]),
        "Data stored in the buffer are most probably wrong");
  }
#endif

  std::vector<CommunicationRequest> requests;
  /// last data to compute the communication scheme
  for (UInt p = 0; p < this->nb_proc; ++p) {
    if (p != this->root) {
      AKANTU_DEBUG_INFO(
          "Sending partition informations to proc "
          << p << " TAG("
          << Tag::genTag(this->rank, this->message_count, Tag::_PARTITIONS)
          << ")");
      requests.push_back(comm.asyncSend(
          buffers[p], p,
          Tag::genTag(this->rank, this->message_count, Tag::_PARTITIONS)));
    }
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

  UInt mesh_data_sizes_buffer_length;
  MeshData & mesh_data = this->getMeshData();

  /// tag info
  std::vector<std::string> tag_names;
  mesh_data.getTagNames(tag_names, type);
  // Make sure the tags are sorted (or at least not in random order),
  // because they come from a map !!
  std::sort(tag_names.begin(), tag_names.end());

  // Sending information about the tags in mesh_data: name, data type and
  // number of components of the underlying array associated to the current
  // type
  DynamicCommunicationBuffer mesh_data_sizes_buffer;
  std::vector<std::string>::const_iterator names_it = tag_names.begin();
  std::vector<std::string>::const_iterator names_end = tag_names.end();
  for (; names_it != names_end; ++names_it) {
    mesh_data_sizes_buffer << *names_it;
    mesh_data_sizes_buffer << mesh_data.getTypeCode(*names_it);
    mesh_data_sizes_buffer << mesh_data.getNbComponent(*names_it, type);
  }

  mesh_data_sizes_buffer_length = mesh_data_sizes_buffer.getSize();
  AKANTU_DEBUG_INFO(
      "Broadcasting the size of the information about the mesh data tags: ("
      << mesh_data_sizes_buffer_length << ").");
  comm.broadcast(mesh_data_sizes_buffer_length, root);
  AKANTU_DEBUG_INFO(
      "Broadcasting the information about the mesh data tags, addr "
      << (void *)mesh_data_sizes_buffer.storage());

  if (mesh_data_sizes_buffer_length != 0)
    comm.broadcast(mesh_data_sizes_buffer, root);

  if (mesh_data_sizes_buffer_length != 0) {
    // Sending the actual data to each processor
    DynamicCommunicationBuffer * buffers =
        new DynamicCommunicationBuffer[nb_proc];
    std::vector<std::string>::const_iterator names_it = tag_names.begin();
    std::vector<std::string>::const_iterator names_end = tag_names.end();

    // Loop over each tag for the current type
    for (; names_it != names_end; ++names_it) {
      // Type code of the current tag (i.e. the tag named *names_it)
      this->fillTagBuffer(buffers, *names_it);
    }

    std::vector<CommunicationRequest> requests;
    for (UInt p = 0; p < nb_proc; ++p) {
      if (p != root) {
        AKANTU_DEBUG_INFO("Sending "
                          << buffers[p].getSize()
                          << " bytes of mesh data to proc " << p << " TAG("
                          << Tag::genTag(this->rank, this->message_count,
                                         Tag::_MESH_DATA)
                          << ")");

        requests.push_back(comm.asyncSend(
            buffers[p], p,
            Tag::genTag(this->rank, this->message_count, Tag::_MESH_DATA)));
      }
    }

    names_it = tag_names.begin();
    // Loop over each tag for the current type
    for (; names_it != names_end; ++names_it) {
      // Reinitializing the mesh data on the master
      this->fillMeshData(buffers[root], *names_it,
                         mesh_data.getTypeCode(*names_it),
                         mesh_data.getNbComponent(*names_it, type));
    }

    comm.waitAll(requests);
    comm.freeCommunicationRequest(requests);
    requests.clear();
    delete[] buffers;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <typename T>
void MasterElementInfoPerProc::fillTagBufferTemplated(
    DynamicCommunicationBuffer * buffers, const std::string & tag_name) {
  MeshData & mesh_data = this->getMeshData();

  const Array<T> & data = mesh_data.getElementalDataArray<T>(tag_name, type);
  const Array<UInt> & partition_num =
      this->partition.getPartition(this->type, _not_ghost);
  const CSR<UInt> & ghost_partition =
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
  const T * data_end = data.storage() + data.getSize() * data.getNbComponent();
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
    CSR<UInt>::const_iterator it = ghost_partition.begin(el);
    CSR<UInt>::const_iterator end = ghost_partition.end(el);
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
    DynamicCommunicationBuffer * buffers, const std::string & tag_name) {
  MeshData & mesh_data = this->getMeshData();

#define AKANTU_DISTRIBUTED_SYNHRONIZER_TAG_DATA(r, extra_param, elem)          \
  case BOOST_PP_TUPLE_ELEM(2, 0, elem): {                                      \
    this->fillTagBufferTemplated<BOOST_PP_TUPLE_ELEM(2, 1, elem)>(buffers,     \
                                                                  tag_name);   \
    break;                                                                     \
  }

  MeshDataTypeCode data_type_code = mesh_data.getTypeCode(tag_name);
  switch (data_type_code) {
    BOOST_PP_SEQ_FOR_EACH(AKANTU_DISTRIBUTED_SYNHRONIZER_TAG_DATA, ,
                          AKANTU_MESH_DATA_TYPES)
  default:
    AKANTU_DEBUG_ERROR("Could not obtain the type of tag" << tag_name << "!");
    break;
  }
#undef AKANTU_DISTRIBUTED_SYNHRONIZER_TAG_DATA
}

/* -------------------------------------------------------------------------- */
void MasterElementInfoPerProc::synchronizeGroups() {
  AKANTU_DEBUG_IN();

  DynamicCommunicationBuffer * buffers =
      new DynamicCommunicationBuffer[nb_proc];

  using ElementToGroup = std::vector<std::vector<std::string>>;
  ElementToGroup element_to_group;
  element_to_group.resize(nb_element);

  auto egi = mesh.element_group_begin();
  auto ege = mesh.element_group_end();
  for (; egi != ege; ++egi) {
    ElementGroup & eg = *(egi->second);

    std::string name = egi->first;

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
  ElementToGroup::const_iterator data_it = element_to_group.begin();
  ElementToGroup::const_iterator data_end = element_to_group.end();
  for (auto pair : zip(partition_num, element_to_group)) {
    buffers[std::get<0>(pair)] << std::get<1>(pair);
  }

  data_it = element_to_group.begin();
  /// copying the data for the ghost element
  for (UInt el(0); data_it != data_end; ++data_it, ++el) {
    CSR<UInt>::const_iterator it = ghost_partition.begin(el);
    CSR<UInt>::const_iterator end = ghost_partition.end(el);
    for (; it != end; ++it) {
      UInt proc = *it;
      buffers[proc] << *data_it;
    }
  }

  std::vector<CommunicationRequest> requests;
  for (UInt p = 0; p < this->nb_proc; ++p) {
    if (p == this->rank)
      continue;
    AKANTU_DEBUG_INFO("Sending element groups to proc "
                      << p << " TAG("
                      << Tag::genTag(this->rank, p, Tag::_ELEMENT_GROUP)
                      << ")");
    requests.push_back(comm.asyncSend(
        buffers[p], p, Tag::genTag(this->rank, p, Tag::_ELEMENT_GROUP)));
  }

  this->fillElementGroupsFromBuffer(buffers[this->rank]);

  comm.waitAll(requests);
  comm.freeCommunicationRequest(requests);
  requests.clear();
  delete[] buffers;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
