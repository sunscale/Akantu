/**
 * @file   communication_tag.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Sep 07 2016
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  Description of the communication tags
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
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_COMMUNICATION_TAG_HH__
#define __AKANTU_COMMUNICATION_TAG_HH__

namespace akantu {
/**
 * tag = |__________20_________|___8____|_4_|
 *       |          proc       | num mes| ct|
 */
class Tag {
public:
  Tag() = default;
  Tag(int val) : tag(val) {}
  Tag(int val, int hash) : tag(val), hash(hash) {}

  operator int() const {
    return int(max_tag == 0 ? tag : (uint32_t(tag) % max_tag));
  }

  /// generates a tag
  template <typename CommTag>
  static inline Tag genTag(int proc, UInt msg_count, CommTag tag) {
    int _tag = ((((proc & 0xFFFFF) << 12) + ((msg_count & 0xFF) << 4) +
                 ((int)tag & 0xF)));
    Tag t(_tag);
    return t;
  }

  /// generates a tag and hashes it
  template <typename CommTag>
  static inline Tag genTag(int proc, UInt msg_count, CommTag tag, int hash) {
    Tag t = genTag(proc, msg_count, tag);
    t.tag = t.tag ^ hash;
    t.hash = hash;
    return t;
  }

  virtual void printself(std::ostream & stream, int) const {
    int t = tag;
    stream << "TAG(";
    if (hash != 0)
      t = t ^ hash;
    stream << (t >> 12) << ":" << (t >> 4 & 0xFF) << ":" << (t & 0xF) << " -> "
           << std::hex << "0x" << int(*this);
    if (hash != 0)
      stream << " {hash: 0x" << hash << "}";
    stream << " [0x" << this->max_tag << "]";
    stream << ")" << std::dec;
  }

  enum CommTags : int {
    _SIZES = 1,
    _CONNECTIVITY = 2,
    _DATA = 3,
    _PARTITIONS = 4,
    _NB_NODES = 5,
    _NODES = 6,
    _COORDINATES = 7,
    _NODES_TYPE = 8,
    _MESH_DATA = 9,
    _ELEMENT_GROUP = 10,
    _NODE_GROUP = 11,
    _MODIFY_SCHEME = 12,
    _GATHER_INITIALIZATION = 1,
    _GATHER = 2,
    _SCATTER = 3,
    _SYNCHRONIZE = 15,
    _REDUCE,
    _PERIODIC_SLAVES,
    _PERIODIC_NODES,
  };

private:
  static void setMaxTag(int _max_tag) { max_tag = _max_tag; }
  friend void initialize(const std::string &, int &, char **&);

private:
  int tag{0};
  int hash{0};
  static int max_tag;
};

/* -------------------------------------------------------------------------- */
inline std::ostream & operator<<(std::ostream & stream, const Tag & _this) {
  _this.printself(stream, 0);
  return stream;
}

/* -------------------------------------------------------------------------- */
} // namespace akantu

#endif /* __AKANTU_COMMUNICATION_TAG_HH__ */
