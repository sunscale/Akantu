/**
 * @file   distributed_synchronizer.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Dana Christen <dana.christen@epfl.ch>
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Sep 01 2010
 * @date last modification: Fri Jan 22 2016
 *
 * @brief  implementation of a  communicator using a static_communicator for
 * real
 * send/receive
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
#include "aka_common.hh"
#include "distributed_synchronizer.hh"
#include "element_info_per_processor.hh"
#include "node_info_per_processor.hh"
#include "static_communicator.hh"
#include "mesh_utils.hh"
/* -------------------------------------------------------------------------- */
#include <map>
#include <iostream>
#include <algorithm>
/* -------------------------------------------------------------------------- */

#if defined(AKANTU_DEBUG_TOOLS)
#include "aka_debug_tools.hh"
#endif

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
DistributedSynchronizer::DistributedSynchronizer(Mesh & mesh,
						 SynchronizerID id,
						 MemoryID memory_id,
						 const bool register_to_event_manager) :
  Synchronizer(id, memory_id),
  mesh(mesh),
  prank_to_element("prank_to_element", id, memory_id)
{

  AKANTU_DEBUG_IN();

  nb_proc = static_communicator->getNbProc();
  rank = static_communicator->whoAmI();

  send_element = new Array<Element>[nb_proc];
  recv_element = new Array<Element>[nb_proc];

  for (UInt p = 0; p < nb_proc; ++p) {
    std::stringstream sstr;
    sstr << p;
    send_element[p].setID(id + ":send_elements_" + sstr.str());
    recv_element[p].setID(id + ":recv_elements_" + sstr.str());
  }

  if (register_to_event_manager)
    mesh.registerEventHandler(*this);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
DistributedSynchronizer::~DistributedSynchronizer() {
  AKANTU_DEBUG_IN();

  for (UInt p = 0; p < nb_proc; ++p) {
    send_element[p].clear();
    recv_element[p].clear();
  }

  delete[] send_element;
  delete[] recv_element;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/// TODO check what it would imply to rewrite this creation as a distributed
/// process
DistributedSynchronizer *
DistributedSynchronizer::createDistributedSynchronizerMesh(
    Mesh & mesh, const MeshPartition * partition, UInt root, SynchronizerID id,
    MemoryID memory_id) {
  AKANTU_DEBUG_IN();

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();

  UInt nb_proc = comm.getNbProc();
  UInt my_rank = comm.whoAmI();

  DistributedSynchronizer & synchronizer =
      *(new DistributedSynchronizer(mesh, id, memory_id));

  if (nb_proc == 1)
    return &synchronizer;


  MeshAccessor mesh_accessor(mesh);
  mesh_accessor.getNodesGlobalIds().resize(0);

  mesh.synchronizeGroupNames();

  /* ------------------------------------------------------------------------ */
  /*  Local (rank == root)                                                    */
  /* ------------------------------------------------------------------------ */
  if (my_rank == root) {
    AKANTU_DEBUG_ASSERT(
        partition->getNbPartition() == nb_proc,
        "The number of partition does not match the number of processors: "
            << partition->getNbPartition() << " != " << nb_proc);

    /**
     * connectivity and communications scheme construction
     */
    Mesh::type_iterator it =
        mesh.firstType(_all_dimensions, _not_ghost, _ek_not_defined);
    Mesh::type_iterator end =
        mesh.lastType(_all_dimensions, _not_ghost, _ek_not_defined);
    UInt count = 0;
    /* --- MAIN LOOP ON TYPES --- */
    for (; it != end; ++it) {
      ElementType type = *it;

      /// \todo change this ugly way to avoid a problem if an element
      /// type is present in the mesh but not in the partitions
      try {
        partition->getPartition(type, _not_ghost);
      } catch (...) {
        continue;
      }

      MasterElementInfoPerProc proc_infos(synchronizer, comm, count, root, mesh, type,
                                          *partition);
      proc_infos.synchronizeConnectivities();
      proc_infos.synchronizePartitions();
      proc_infos.synchronizeTags();
      proc_infos.synchronizeGroups();
      ++count;
    }

    { /// Ending the synchronization of elements by sending a stop message
      MasterElementInfoPerProc proc_infos(synchronizer, comm, count, root, mesh,
                                          _not_defined, *partition);
      ++count;
    }

    /**
     * Nodes synchronization
     */
    MasterNodeInfoPerProc node_proc_infos(synchronizer, comm, count, root, mesh);

    node_proc_infos.synchronizeNodes();
    node_proc_infos.synchronizeTypes();
    node_proc_infos.synchronizeGroups();

    /* ---------------------------------------------------------------------- */
    /*  Distant (rank != root)                                                */
    /* ---------------------------------------------------------------------- */
  } else {
    /**
     * connectivity and communications scheme construction on distant processors
     */
    UInt count = 0;
    bool need_synchronize = true;
    do {
      /* --------<<<<-SIZE--------------------------------------------------- */
      SlaveElementInfoPerProc proc_infos(synchronizer, comm, count, root, mesh);
      need_synchronize = proc_infos.needSynchronize();

      if (need_synchronize) {
        proc_infos.synchronizeConnectivities();
        proc_infos.synchronizePartitions();
        proc_infos.synchronizeTags();
        proc_infos.synchronizeGroups();
      }
      ++count;
    } while (need_synchronize);

    /**
     * Nodes synchronization
     */

    SlaveNodeInfoPerProc node_proc_infos(synchronizer, comm, count, root, mesh);

    node_proc_infos.synchronizeNodes();
    node_proc_infos.synchronizeTypes();
    node_proc_infos.synchronizeGroups();
  }

  MeshUtils::fillElementToSubElementsData(mesh);

  mesh_accessor.setDistributed();

  AKANTU_DEBUG_OUT();
  return &synchronizer;
}

/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::asynchronousSynchronize(
    DataAccessor & data_accessor, SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  if (communications.find(tag) == communications.end())
    computeBufferSize(data_accessor, tag);

  Communication & communication = communications[tag];

  AKANTU_DEBUG_ASSERT(
      communication.send_requests.size() == 0,
      "There must be some pending sending communications. Tag is " << tag);

  std::map<SynchronizationTag, UInt>::iterator t_it = tag_counter.find(tag);
  UInt counter = 0;
  if (t_it == tag_counter.end()) {
    tag_counter[tag] = 0;
  } else {
    counter = ++(t_it->second);
  }

  for (UInt p = 0; p < nb_proc; ++p) {
    UInt ssize = communication.size_to_send[p];
    if (p == rank || ssize == 0)
      continue;

    CommunicationBuffer & buffer = communication.send_buffer[p];
    buffer.resize(ssize);

    Tag comm_tag = Tag::genTag(rank, counter, tag);
    buffer << int(comm_tag);

#ifndef AKANTU_NDEBUG
    UInt nb_elements = send_element[p].getSize();
    AKANTU_DEBUG_INFO("Packing data for proc " << p << " (" << ssize << "/"
                                               << nb_elements
                                               << " data to send/elements)");

    /// pack barycenters in debug mode
    Array<Element>::const_iterator<Element> bit = send_element[p].begin();
    Array<Element>::const_iterator<Element> bend = send_element[p].end();
    for (; bit != bend; ++bit) {
      const Element & element = *bit;
      Vector<Real> barycenter(mesh.getSpatialDimension());
      mesh.getBarycenter(element.element, element.type, barycenter.storage(),
                         element.ghost_type);
      buffer << barycenter;
    }
#endif

    data_accessor.packElementData(buffer, send_element[p], tag);

    AKANTU_DEBUG_ASSERT(buffer.getPackedSize() == ssize,
                        "a problem have been introduced with "
                            << "false sent sizes declaration "
                            << buffer.getPackedSize() << " != " << ssize);
    AKANTU_DEBUG_INFO("Posting send to proc "
                      << p << " (tag: " << tag << " - " << ssize
                      << " data to send)"
                      << " [ " << comm_tag << ":" << std::hex
                      << int(this->genTagFromID(tag)) << " ]");
    communication.send_requests.push_back(static_communicator->asyncSend(
        buffer, p, this->genTagFromID(tag)));
  }

  AKANTU_DEBUG_ASSERT(communication.recv_requests.size() == 0,
                      "There must be some pending receive communications");

  for (UInt p = 0; p < nb_proc; ++p) {
    UInt rsize = communication.size_to_receive[p];
    if (p == rank || rsize == 0)
      continue;
    CommunicationBuffer & buffer = communication.recv_buffer[p];
    buffer.resize(rsize);

    Tag comm_tag = Tag::genTag(rank, counter, tag);
    buffer << int(comm_tag);

    AKANTU_DEBUG_INFO("Posting receive from proc "
                      << p << " (tag: " << tag << " - " << rsize
                      << " data to receive) "
                      << " [ " << comm_tag << ":" << std::hex
                      << int(this->genTagFromID(tag)) << " ]");
    communication.recv_requests.push_back(static_communicator->asyncReceive(
        buffer, p, this->genTagFromID(tag)));
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::waitEndSynchronize(DataAccessor & data_accessor,
                                                 SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(communications.find(tag) != communications.end(),
                      "No communication with the tag \"" << tag
                                                         << "\" started");

  Communication & communication = communications[tag];

  std::vector<CommunicationRequest *> req_not_finished;
  std::vector<CommunicationRequest *> * req_not_finished_tmp =
      &req_not_finished;
  std::vector<CommunicationRequest *> * recv_requests_tmp =
      &(communication.recv_requests);

  //  static_communicator->waitAll(recv_requests);
  while (!recv_requests_tmp->empty()) {
    for (std::vector<CommunicationRequest *>::iterator req_it =
             recv_requests_tmp->begin();
         req_it != recv_requests_tmp->end(); ++req_it) {
      CommunicationRequest * req = *req_it;

      if (static_communicator->testRequest(req)) {
        UInt proc = req->getSource();
        AKANTU_DEBUG_INFO("Unpacking data coming from proc " << proc);
        CommunicationBuffer & buffer = communication.recv_buffer[proc];

        int _tag;
        buffer >> _tag;
        Tag comm_tag(_tag);

#ifndef AKANTU_NDEBUG
        Array<Element>::const_iterator<Element> bit =
            recv_element[proc].begin();
        Array<Element>::const_iterator<Element> bend = recv_element[proc].end();

        UInt spatial_dimension = mesh.getSpatialDimension();

        for (; bit != bend; ++bit) {
          const Element & element = *bit;

          Vector<Real> barycenter_loc(spatial_dimension);
          mesh.getBarycenter(element.element, element.type,
                             barycenter_loc.storage(), element.ghost_type);
          Vector<Real> barycenter(spatial_dimension);
          buffer >> barycenter;
          for (UInt i = 0; i < spatial_dimension; ++i) {
            if (! Math::are_float_equal(barycenter_loc(i), barycenter(i)))
	      AKANTU_DEBUG_ERROR("Unpacking an unknown value for the element: "
				 << element << "(barycenter[" << i << "] = " << barycenter_loc(i)
				 << " and buffer[" << i << "] = " << barycenter(i) << ") ["
				 << std::abs(barycenter(i) - barycenter_loc(i))
				 << "] - tag: " << tag << " comm_tag[ " << comm_tag << " ]");
          }
        }
#endif

        data_accessor.unpackElementData(buffer, recv_element[proc], tag);
        buffer.resize(0);

        AKANTU_DEBUG_ASSERT(buffer.getLeftToUnpack() == 0,
                            "all data have not been unpacked: "
                                << buffer.getLeftToUnpack() << " bytes left"
                                << " [ " << comm_tag << " ]");
        static_communicator->freeCommunicationRequest(req);
      } else {
        req_not_finished_tmp->push_back(req);
      }
    }

    std::vector<CommunicationRequest *> * swap = req_not_finished_tmp;
    req_not_finished_tmp = recv_requests_tmp;
    recv_requests_tmp = swap;

    req_not_finished_tmp->clear();
  }

  AKANTU_DEBUG_INFO("Waiting that every send requests are received");
  static_communicator->waitAll(communication.send_requests);
  for (std::vector<CommunicationRequest *>::iterator req_it =
           communication.send_requests.begin();
       req_it != communication.send_requests.end(); ++req_it) {
    CommunicationRequest & req = *(*req_it);

    if (static_communicator->testRequest(&req)) {
      UInt proc = req.getDestination();
      CommunicationBuffer & buffer = communication.send_buffer[proc];
      buffer.resize(0);
      static_communicator->freeCommunicationRequest(&req);
    }
  }
  communication.send_requests.clear();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::computeBufferSize(DataAccessor & data_accessor,
                                                SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  communications[tag].resize(nb_proc);

  for (UInt p = 0; p < nb_proc; ++p) {
    UInt ssend = 0;
    UInt sreceive = 0;
    if (p != rank) {
      if (send_element[p].getSize() != 0) {
        ssend += sizeof(int); // sizeof(int) is for the communication tag
#ifndef AKANTU_NDEBUG
        ssend += send_element[p].getSize() * mesh.getSpatialDimension() *
                 sizeof(Real);
#endif
        ssend += data_accessor.getNbDataForElements(send_element[p], tag);
        AKANTU_DEBUG_INFO("I have " << ssend << "(" << ssend / 1024. << "kB - "
                                    << send_element[p].getSize()
                                    << " element(s)) data to send to " << p
                                    << " for tag " << tag);
      }

      if (recv_element[p].getSize() != 0) {
        sreceive += sizeof(int); // sizeof(int) is for the communication tag
#ifndef AKANTU_NDEBUG
        sreceive += recv_element[p].getSize() * mesh.getSpatialDimension() *
                    sizeof(Real);
#endif
        sreceive += data_accessor.getNbDataForElements(recv_element[p], tag);
        AKANTU_DEBUG_INFO("I have " << sreceive << "(" << sreceive / 1024.
                                    << "kB - " << recv_element[p].getSize()
                                    << " element(s)) data to receive for tag "
                                    << tag);
      }
    }

    communications[tag].size_to_send[p] = ssend;
    communications[tag].size_to_receive[p] = sreceive;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::computeAllBufferSizes(
    DataAccessor & data_accessor) {
  std::map<SynchronizationTag, Communication>::iterator it =
      this->communications.begin();
  std::map<SynchronizationTag, Communication>::iterator end =
      this->communications.end();

  for (; it != end; ++it) {
    SynchronizationTag tag = it->first;
    this->computeBufferSize(data_accessor, tag);
  }
}

/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::printself(std::ostream & stream,
                                        int indent) const {
  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
    ;

  Int prank = StaticCommunicator::getStaticCommunicator().whoAmI();
  Int psize = StaticCommunicator::getStaticCommunicator().getNbProc();
  stream << "[" << prank << "/" << psize << "]" << space
         << "DistributedSynchronizer [" << std::endl;
  for (UInt p = 0; p < nb_proc; ++p) {
    if (p == UInt(prank))
      continue;
    stream << "[" << prank << "/" << psize << "]" << space
           << " + Communication to proc " << p << " [" << std::endl;
    if (AKANTU_DEBUG_TEST(dblDump)) {
      stream << "[" << prank << "/" << psize << "]" << space
             << "    - Element to send to proc " << p << " [" << std::endl;

      Array<Element>::iterator<Element> it_el = send_element[p].begin();
      Array<Element>::iterator<Element> end_el = send_element[p].end();
      for (; it_el != end_el; ++it_el)
        stream << "[" << prank << "/" << psize << "]" << space << "       "
               << *it_el << std::endl;
      stream << "[" << prank << "/" << psize << "]" << space << "   ]"
             << std::endl;

      stream << "[" << prank << "/" << psize << "]" << space
             << "    - Element to recv from proc " << p << " [" << std::endl;

      it_el = recv_element[p].begin();
      end_el = recv_element[p].end();
      for (; it_el != end_el; ++it_el)
        stream << "[" << prank << "/" << psize << "]" << space << "       "
               << *it_el << std::endl;

      stream << "[" << prank << "/" << psize << "]" << space << "   ]"
             << std::endl;
    }

    std::map<SynchronizationTag, Communication>::const_iterator it =
        communications.begin();
    std::map<SynchronizationTag, Communication>::const_iterator end =
        communications.end();
    for (; it != end; ++it) {
      const SynchronizationTag & tag = it->first;
      const Communication & communication = it->second;
      UInt ssend = communication.size_to_send[p];
      UInt sreceive = communication.size_to_receive[p];
      stream << "[" << prank << "/" << psize << "]" << space << "     - Tag "
             << tag << " -> " << ssend << "byte(s) -- <- " << sreceive
             << "byte(s)" << std::endl;
    }
  }
}

/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::substituteElements(
    const std::map<Element, Element> & old_to_new_elements) {
  // substitute old elements with new ones
  std::map<Element, Element>::const_iterator found_element_it;
  std::map<Element, Element>::const_iterator found_element_end =
      old_to_new_elements.end();

  for (UInt p = 0; p < nb_proc; ++p) {
    if (p == rank)
      continue;

    Array<Element> & recv = recv_element[p];
    for (UInt el = 0; el < recv.getSize(); ++el) {
      found_element_it = old_to_new_elements.find(recv(el));
      if (found_element_it != found_element_end)
        recv(el) = found_element_it->second;
    }

    Array<Element> & send = send_element[p];
    for (UInt el = 0; el < send.getSize(); ++el) {
      found_element_it = old_to_new_elements.find(send(el));
      if (found_element_it != found_element_end)
        send(el) = found_element_it->second;
    }
  }
}

/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::onElementsChanged(
    const Array<Element> & old_elements_list,
    const Array<Element> & new_elements_list,
    __attribute__((unused)) const ElementTypeMapArray<UInt> & new_numbering,
    __attribute__((unused)) const ChangedElementsEvent & event) {
  // create a map to link old elements to new ones
  std::map<Element, Element> old_to_new_elements;

  for (UInt el = 0; el < old_elements_list.getSize(); ++el) {
    AKANTU_DEBUG_ASSERT(old_to_new_elements.find(old_elements_list(el)) ==
                            old_to_new_elements.end(),
                        "The same element cannot appear twice in the list");

    old_to_new_elements[old_elements_list(el)] = new_elements_list(el);
  }

  substituteElements(old_to_new_elements);
}

/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::onElementsRemoved(
    const Array<Element> & element_to_remove,
    const ElementTypeMapArray<UInt> & new_numbering,
    __attribute__((unused)) const RemovedElementsEvent & event) {
  AKANTU_DEBUG_IN();
  this->removeElements(element_to_remove);
  this->renumberElements(new_numbering);
  AKANTU_DEBUG_OUT();
}

// void DistributedSynchronizer::checkCommunicationScheme() {
//   for (UInt p = 0; p < psize; ++p) {
//     if (p == prank) continue;
//     for(UInt e(0), e < recv_element.getSize())
//   }
// }

/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::buildPrankToElement() {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();

  mesh.initElementTypeMapArray(prank_to_element, 1, spatial_dimension, false,
                               _ek_not_defined, true);

  Mesh::type_iterator it =
      mesh.firstType(spatial_dimension, _not_ghost, _ek_not_defined);

  Mesh::type_iterator end =
      mesh.lastType(spatial_dimension, _not_ghost, _ek_not_defined);

  /// assign prank to all not ghost elements
  for (; it != end; ++it) {
    UInt nb_element = mesh.getNbElement(*it);
    Array<UInt> & prank_to_el = prank_to_element(*it);
    for (UInt el = 0; el < nb_element; ++el) {
      prank_to_el(el) = rank;
    }
  }

  /// assign prank to all ghost elements
  for (UInt p = 0; p < nb_proc; ++p) {
    UInt nb_ghost_element = recv_element[p].getSize();

    for (UInt el = 0; el < nb_ghost_element; ++el) {
      UInt element = recv_element[p](el).element;
      ElementType type = recv_element[p](el).type;
      GhostType ghost_type = recv_element[p](el).ghost_type;

      Array<UInt> & prank_to_el = prank_to_element(type, ghost_type);
      prank_to_el(element) = p;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::filterElementsByKind(
    DistributedSynchronizer * new_synchronizer, ElementKind kind) {
  AKANTU_DEBUG_IN();

  Array<Element> * newsy_send_element = new_synchronizer->send_element;
  Array<Element> * newsy_recv_element = new_synchronizer->recv_element;

  Array<Element> * new_send_element = new Array<Element>[nb_proc];
  Array<Element> * new_recv_element = new Array<Element>[nb_proc];

  for (UInt p = 0; p < nb_proc; ++p) {

    /// send element copying part
    new_send_element[p].resize(0);

    for (UInt el = 0; el < send_element[p].getSize(); ++el) {
      Element & element = send_element[p](el);

      if (element.kind == kind)
        newsy_send_element[p].push_back(element);
      else
        new_send_element[p].push_back(element);
    }

    /// recv element copying part
    new_recv_element[p].resize(0);

    for (UInt el = 0; el < recv_element[p].getSize(); ++el) {
      Element & element = recv_element[p](el);

      if (element.kind == kind)
        newsy_recv_element[p].push_back(element);
      else
        new_recv_element[p].push_back(element);
    }
  }

  /// deleting and reassigning old pointers
  delete[] send_element;
  delete[] recv_element;

  send_element = new_send_element;
  recv_element = new_recv_element;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::reset() {
  AKANTU_DEBUG_IN();

  for (UInt p = 0; p < nb_proc; ++p) {
    send_element[p].resize(0);
    recv_element[p].resize(0);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::removeElements(
    const Array<Element> & element_to_remove) {
  AKANTU_DEBUG_IN();

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();

  std::vector<CommunicationRequest *> isend_requests;
  Array<UInt> * list_of_el = new Array<UInt>[nb_proc];
  // Handling ghost elements
  for (UInt p = 0; p < nb_proc; ++p) {
    if (p == rank)
      continue;

    Array<Element> & recv = recv_element[p];
    if (recv.getSize() == 0)
      continue;

    Array<Element>::iterator<Element> recv_begin = recv.begin();
    Array<Element>::iterator<Element> recv_end = recv.end();

    Array<Element>::const_iterator<Element> er_it = element_to_remove.begin();
    Array<Element>::const_iterator<Element> er_end = element_to_remove.end();

    Array<UInt> & list = list_of_el[p];
    for (UInt i = 0; recv_begin != recv_end; ++i, ++recv_begin) {
      const Element & el = *recv_begin;
      Array<Element>::const_iterator<Element> pos =
          std::find(er_it, er_end, el);
      if (pos == er_end) {
        list.push_back(i);
      }
    }

    if (list.getSize() == recv.getSize())
      list.push_back(UInt(0));
    else
      list.push_back(UInt(-1));

    AKANTU_DEBUG_INFO("Sending a message of size "
                      << list.getSize() << " to proc " << p << " TAG("
                      << this->genTagFromID(0) << ")");
    isend_requests.push_back(comm.asyncSend(list, p, this->genTagFromID(0)));

    list.erase(list.getSize() - 1);

    if (list.getSize() == recv.getSize())
      continue;

    Array<Element> new_recv;
    for (UInt nr = 0; nr < list.getSize(); ++nr) {
      Element & el = recv(list(nr));
      new_recv.push_back(el);
    }

    AKANTU_DEBUG_INFO("I had " << recv.getSize()
                               << " elements to recv from proc " << p << " and "
                               << list.getSize() << " elements to keep. I have "
                               << new_recv.getSize() << " elements left.");
    recv.copy(new_recv);
  }

  for (UInt p = 0; p < nb_proc; ++p) {
    if (p == rank)
      continue;
    Array<Element> & send = send_element[p];

    if (send.getSize() == 0)
      continue;

    CommunicationStatus status;
    AKANTU_DEBUG_INFO("Getting number of elements of proc "
                      << p << " not needed anymore TAG("
                      << this->genTagFromID(0) << ")");
    comm.probe<UInt>(p, this->genTagFromID(0), status);
    Array<UInt> list(status.getSize());

    AKANTU_DEBUG_INFO("Receiving list of elements ("
                      << status.getSize() - 1
                      << " elements) no longer needed by proc " << p << " TAG("
                      << this->genTagFromID(0) << ")");
    comm.receive(list, p, this->genTagFromID(0));

    if (list.getSize() == 1 && list(0) == 0)
      continue;

    list.erase(list.getSize() - 1);

    if (list.getSize() == send.getSize())
      continue;

    Array<Element> new_send;
    for (UInt ns = 0; ns < list.getSize(); ++ns) {
      Element & el = send(list(ns));
      new_send.push_back(el);
    }

    AKANTU_DEBUG_INFO("I had " << send.getSize() << " elements to send to proc "
                               << p << " and " << list.getSize()
                               << " elements to keep. I have "
                               << new_send.getSize() << " elements left.");
    send.copy(new_send);
  }

  comm.waitAll(isend_requests);
  comm.freeCommunicationRequest(isend_requests);

  delete[] list_of_el;
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DistributedSynchronizer::renumberElements(
    const ElementTypeMapArray<UInt> & new_numbering) {
  // Handling ghost elements
  for (UInt p = 0; p < nb_proc; ++p) {
    if (p == rank)
      continue;

    Array<Element> & recv = recv_element[p];

    for (UInt i = 0; i < recv.getSize(); ++i) {
      Element & el = recv(i);
      el.element = new_numbering(el.type, el.ghost_type)(el.element);
    }

    Array<Element> & send = send_element[p];

    for (UInt i = 0; i < send.getSize(); ++i) {
      Element & el = send(i);
      el.element = new_numbering(el.type, el.ghost_type)(el.element);
    }
  }
}

__END_AKANTU__
