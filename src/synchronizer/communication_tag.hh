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

#ifndef AKANTU_COMMUNICATION_TAG_HH_
#define AKANTU_COMMUNICATION_TAG_HH_

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

  virtual void printself(std::ostream & stream, int /*unused*/) const {
    int t = tag;
    stream << "TAG(";
    if (hash != 0) {
      t = t ^ hash;
    }
    stream << (t >> 12) << ":" << (t >> 4 & 0xFF) << ":" << (t & 0xF) << " -> "
           << std::hex << "0x" << int(*this);
    if (hash != 0) {
      stream << " {hash: 0x" << hash << "}";
    }
    stream << " [0x" << max_tag << "]";
    stream << ")" << std::dec;
  }

  enum CommTags : int {
    _sizes = 1,
    _connectivity = 2,
    _data = 3,
    _partitions = 4,
    _nb_nodes = 5,
    _nodes = 6,
    _coordinates = 7,
    _nodes_type = 8,
    _mesh_data = 9,
    _element_group = 10,
    _node_group = 11,
    _modify_scheme = 12,
    _gather_initialization = 1,
    _gather = 2,
    _scatter = 3,
    _synchronize = 15,
    _reduce,
    _periodic_slaves,
    _periodic_nodes,
  };

private:
  static void setMaxTag(int _max_tag) { max_tag = _max_tag; }
  friend void initialize(const std::string & /*input_file*/, int & /*argc*/,
                         char **& /*argv*/);

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

#endif /* AKANTU_COMMUNICATION_TAG_HH_ */
