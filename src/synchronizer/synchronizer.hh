/**
 * @file   synchronizer.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Thu Dec 10 2015
 *
 * @brief  Common interface for synchronizers
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

#ifndef __AKANTU_SYNCHRONIZER_HH__
#define __AKANTU_SYNCHRONIZER_HH__

/* -------------------------------------------------------------------------- */
#include "aka_memory.hh"
#include "data_accessor.hh"
#include "real_static_communicator.hh"
#include "static_communicator.hh"

/* -------------------------------------------------------------------------- */
#include <map>
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/**
 * tag = |__________20_________|___8____|_4_|
 *       |          proc       | num mes| ct|
 */
class Tag {
public:
  Tag() : tag(0) {}
  Tag(int val) : tag(val) {}

  operator int() { return int(tag); } // remove the sign bit

  template <typename CommTag>
  static inline Tag genTag(int proc, UInt msg_count, CommTag tag) {
    int max_tag = StaticCommunicator::getStaticCommunicator().getMaxTag();
    int _tag = ((((proc & 0xFFFFF) << 12) + ((msg_count & 0xFF) << 4) +
                 ((Int)tag & 0xF)));
    Tag t(max_tag == 0 ? _tag : (_tag % max_tag));
    return t;
  }

  virtual void printself(std::ostream & stream,
                         __attribute__((unused)) int indent = 0) const {
    stream << (tag >> 12) << ":" << (tag >> 4 & 0xFF) << ":" << (tag & 0xF);
  }


  enum CommTags {
    _SIZES = 0,
    _CONNECTIVITY = 1,
    _DATA = 2,
    _PARTITIONS = 3,
    _NB_NODES = 4,
    _NODES = 5,
    _COORDINATES = 6,
    _NODES_TYPE = 7,
    _MESH_DATA = 8,
    _ELEMENT_GROUP = 9,
    _NODE_GROUP = 10,
  };
private:
  int tag;
};
/* -------------------------------------------------------------------------- */

class Synchronizer : protected Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  Synchronizer(SynchronizerID id = "synchronizer", MemoryID memory_id = 0,
	       StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator());

  virtual ~Synchronizer(){};

  virtual void printself(__attribute__((unused)) std::ostream & stream,
                         __attribute__((unused)) int indent = 0) const {};
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// synchronize ghosts
  void synchronize(DataAccessor & data_accessor, SynchronizationTag tag);

  /// asynchronous synchronization of ghosts
  virtual void asynchronousSynchronize(DataAccessor & data_accessor,
                                       SynchronizationTag tag) = 0;

  /// wait end of asynchronous synchronization of ghosts
  virtual void waitEndSynchronize(DataAccessor & data_accessor,
                                  SynchronizationTag tag) = 0;

  /// compute buffer size for a given tag and data accessor
  virtual void computeBufferSize(DataAccessor & data_accessor,
                                 SynchronizationTag tag) = 0;

  /// generate the tag from the ID
  template <typename CommTag> inline Tag genTagFromID(CommTag tag) {
    int max_tag = StaticCommunicator::getStaticCommunicator().getMaxTag();
    int _tag =
        std::abs((int(hash<std::string>(this->getID())) << 4) + (tag & 0xF));
    Tag t(max_tag == 0 ? _tag : (_tag % max_tag));
    return t;
  }

protected:
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  class Communication {
  public:
    void resize(UInt size) {
      send_buffer.resize(size);
      recv_buffer.resize(size);
      size_to_send.resize(size);
      size_to_receive.resize(size);
    }

  public:
    /// size of data to send to each processor
    std::vector<UInt> size_to_send;
    /// size of data to recv to each processor
    std::vector<UInt> size_to_receive;
    std::vector<CommunicationBuffer> send_buffer;
    std::vector<CommunicationBuffer> recv_buffer;

    std::vector<CommunicationRequest *> send_requests;
    std::vector<CommunicationRequest *> recv_requests;
  };

  /// id of the synchronizer
  SynchronizerID id;

  /// message counter per tag
  std::map<SynchronizationTag, UInt> tag_counter;

  /// the static memory instance
  StaticCommunicator * static_communicator;
};

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream,
                                 const Synchronizer & _this) {
  _this.printself(stream);
  return stream;
}

inline std::ostream & operator<<(std::ostream & stream,
                                 const Tag & _this) {
  _this.printself(stream);
  return stream;
}

__END_AKANTU__

#endif /* __AKANTU_SYNCHRONIZER_HH__ */
