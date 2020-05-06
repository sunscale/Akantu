/**
 * @file   slave_element_info_per_processor.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Mar 16 2016
 * @date last modification: Tue Nov 07 2017
 *
 * @brief  Helper class to distribute a mesh
 *
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
#include "communicator.hh"
#include "element_info_per_processor.hh"
#include "element_synchronizer.hh"
#include "mesh_utils.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
#include <iostream>
#include <map>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
SlaveElementInfoPerProc::SlaveElementInfoPerProc(
    ElementSynchronizer & synchronizer, UInt message_cnt, UInt root)
    : ElementInfoPerProc(synchronizer, message_cnt, root, _not_defined) {

  Vector<UInt> size(5);
  comm.receive(size, this->root,
               Tag::genTag(this->root, this->message_count, Tag::_SIZES));

  this->type = (ElementType)size[0];
  this->nb_local_element = size[1];
  this->nb_ghost_element = size[2];
  this->nb_element_to_receive = size[3];
  this->nb_tags = size[4];

  if (this->type != _not_defined)
    this->nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
}

/* -------------------------------------------------------------------------- */

bool SlaveElementInfoPerProc::needSynchronize() {
  return this->type != _not_defined;
}

/* -------------------------------------------------------------------------- */
void SlaveElementInfoPerProc::synchronizeConnectivities() {
  Array<UInt> local_connectivity(
      (this->nb_local_element + this->nb_ghost_element) *
      this->nb_nodes_per_element);

  AKANTU_DEBUG_INFO("Receiving connectivities from proc " << root);
  comm.receive(
      local_connectivity, this->root,
      Tag::genTag(this->root, this->message_count, Tag::_CONNECTIVITY));

  auto & old_nodes = this->getNodesGlobalIds();
  AKANTU_DEBUG_INFO("Renumbering local connectivities");
  MeshUtils::renumberMeshNodes(this->mesh, local_connectivity,
                               this->nb_local_element, this->nb_ghost_element,
                               this->type, old_nodes);
}

/* -------------------------------------------------------------------------- */
void SlaveElementInfoPerProc::synchronizePartitions() {
  Array<UInt> local_partitions(this->nb_element_to_receive +
                               this->nb_ghost_element * 2);
  AKANTU_DEBUG_INFO("Receiving partition informations from proc " << root);
  this->comm.receive(local_partitions, this->root,
                     Tag::genTag(root, this->message_count, Tag::_PARTITIONS));

  if (Mesh::getSpatialDimension(this->type) ==
      this->mesh.getSpatialDimension()) {
    AKANTU_DEBUG_INFO("Creating communications scheme");
    this->fillCommunicationScheme(local_partitions);
  }
}

/* -------------------------------------------------------------------------- */
void SlaveElementInfoPerProc::synchronizeTags() {
  AKANTU_DEBUG_IN();

  if (this->nb_tags == 0) {
    AKANTU_DEBUG_OUT();
    return;
  }

  /* --------<<<<-TAGS------------------------------------------------- */
  DynamicCommunicationBuffer mesh_data_sizes_buffer;
  comm.broadcast(mesh_data_sizes_buffer, root);
  AKANTU_DEBUG_INFO("Size of the information about the mesh data: "
                    << mesh_data_sizes_buffer.size());

  if (mesh_data_sizes_buffer.size() == 0)
    return;

  AKANTU_DEBUG_INFO("Receiving the information about the mesh data tags, addr "
                    << (void *)mesh_data_sizes_buffer.storage());

  std::vector<std::string> tag_names;
  std::vector<MeshDataTypeCode> tag_type_codes;
  std::vector<UInt> tag_nb_component;
  tag_names.resize(nb_tags);
  tag_type_codes.resize(nb_tags);
  tag_nb_component.resize(nb_tags);
  CommunicationBuffer mesh_data_buffer;
  UInt type_code_int;
  for (UInt i(0); i < nb_tags; ++i) {
    mesh_data_sizes_buffer >> tag_names[i];
    mesh_data_sizes_buffer >> type_code_int;
    tag_type_codes[i] = static_cast<MeshDataTypeCode>(type_code_int);
    mesh_data_sizes_buffer >> tag_nb_component[i];
  }

  std::vector<std::string>::const_iterator names_it = tag_names.begin();
  std::vector<std::string>::const_iterator names_end = tag_names.end();

  CommunicationStatus mesh_data_comm_status;
  AKANTU_DEBUG_INFO("Checking size of data to receive for mesh data TAG("
                    << Tag::genTag(root, this->message_count, Tag::_MESH_DATA)
                    << ")");
  comm.probe<char>(root,
                   Tag::genTag(root, this->message_count, Tag::_MESH_DATA),
                   mesh_data_comm_status);
  UInt mesh_data_buffer_size(mesh_data_comm_status.size());
  AKANTU_DEBUG_INFO("Receiving "
                    << mesh_data_buffer_size << " bytes of mesh data TAG("
                    << Tag::genTag(root, this->message_count, Tag::_MESH_DATA)
                    << ")");
  mesh_data_buffer.resize(mesh_data_buffer_size);
  comm.receive(mesh_data_buffer, root,
               Tag::genTag(root, this->message_count, Tag::_MESH_DATA));

  // Loop over each tag for the current type
  UInt k(0);
  for (; names_it != names_end; ++names_it, ++k) {
    this->fillMeshData(mesh_data_buffer, *names_it, tag_type_codes[k],
                       tag_nb_component[k]);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SlaveElementInfoPerProc::synchronizeGroups() {
  AKANTU_DEBUG_IN();

  const Communicator & comm = mesh.getCommunicator();
  UInt my_rank = comm.whoAmI();

  AKANTU_DEBUG_INFO("Receiving element groups from proc "
                    << root << " TAG("
                    << Tag::genTag(root, my_rank, Tag::_ELEMENT_GROUP) << ")");

  CommunicationStatus status;
  comm.probe<char>(root, Tag::genTag(root, my_rank, Tag::_ELEMENT_GROUP),
                   status);

  CommunicationBuffer buffer(status.size());
  comm.receive(buffer, root, Tag::genTag(root, my_rank, Tag::_ELEMENT_GROUP));

  this->fillElementGroupsFromBuffer(buffer);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
