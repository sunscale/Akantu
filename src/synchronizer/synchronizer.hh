/**
 * @file   synchronizer.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Sep 01 2010
 * @date last modification: Tue Apr 30 2013
 *
 * @brief  interface for communicator and pbc synchronizers
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

#ifndef __AKANTU_SYNCHRONIZER_HH__
#define __AKANTU_SYNCHRONIZER_HH__

/* -------------------------------------------------------------------------- */
#include "aka_memory.hh"
#include "data_accessor.hh"
/* -------------------------------------------------------------------------- */
#include <map>
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class Synchronizer : protected Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Synchronizer(SynchronizerID id = "synchronizer", MemoryID memory_id = 0);

  virtual ~Synchronizer() { };

  virtual void printself(__attribute__((unused)) std::ostream & stream,
			 __attribute__((unused)) int indent = 0) const {};
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// synchronize ghosts
  void synchronize(DataAccessor & data_accessor,SynchronizationTag tag);

  /// asynchronous synchronization of ghosts
  virtual void asynchronousSynchronize(DataAccessor & data_accessor,SynchronizationTag tag) = 0;

  /// wait end of asynchronous synchronization of ghosts
  virtual void waitEndSynchronize(DataAccessor & data_accessor,SynchronizationTag tag) = 0;

  /// compute buffer size for a given tag and data accessor
  virtual void computeBufferSize(DataAccessor & data_accessor, SynchronizationTag tag)=0;

  /**
   * tag = |__________20_________|___8____|_4_|
   *       |          proc       | num mes| ct|
   */
  class Tag {
  public:
    operator int() { return int(tag); }

    template<typename CommTag>
    static inline Tag genTag(int proc, UInt msg_count, CommTag tag) {
      Tag t;
      t.tag = (((proc & 0xFFFFF) << 12) + ((msg_count & 0xFF) << 4) + ((Int)tag & 0xF));
      return t;
    }

    virtual void printself(std::ostream & stream,
			   __attribute__((unused)) int indent = 0) const {
      stream << (tag >> 12) << ":" << (tag >> 4 & 0xFF) << ":" << (tag & 0xF);
    }
  private:
    UInt tag;
  };

protected:

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// id of the synchronizer
  SynchronizerID id;

  /// message counter per tag
  std::map<SynchronizationTag, UInt> tag_counter;
};


/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const Synchronizer & _this)
{
  _this.printself(stream);
  return stream;
}

inline std::ostream & operator <<(std::ostream & stream, const Synchronizer::Tag & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__

#endif /* __AKANTU_SYNCHRONIZER_HH__ */
