/**
 * @file   dumper_lammps.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date creation: Thu Nov 25 2010
 * @date last modification: Tue Jun 04 2013
 *
 * @brief  header for lammps dumper
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * IOHelper is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * IOHelper is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with IOHelper. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#ifndef IOHELPER_DUMPER_LAMMPS_H_
#define IOHELPER_DUMPER_LAMMPS_H_
/* -------------------------------------------------------------------------- */
#include "dumper.hh"
#include <fstream>
/* -------------------------------------------------------------------------- */


namespace iohelper {

enum LammpsAtomStyle {atomic, bond}; //please extend ad libidum

template<LammpsAtomStyle style>
class DumperLammps: public Dumper, public Visitor {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  DumperLammps(Real * bounds = nullptr, const std::string & prefix = "./");

  ~DumperLammps() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  //! dump to file
  void dump(const std::string & current_name = std::string(),
            UInt count = UInt(-1)) override;
  void dumpHead(Real * bounds = nullptr);
  template<typename T>
  void visitField(T & visited);

  void dumpFinalize();
  //! set mode for file creation : TEXT, BASE64, COMPRESSED
  void setMode(int mode) override { Dumper::setMode(mode); }
  void dumpAdd(int grain_id = 1);
  void setEmbeddedValue(__attribute__((unused)) const std::string & name,
			__attribute__((unused)) int value){}

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                 */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  //position of where the number of atoms is printed;
  std::streampos nb_atom_position;
  //current number of atoms printed to the file
  unsigned long int curr_nb_atom;
  std::fstream lammps_dump_file;
  Real * bounds;

  //! flag to produce zipped files
  bool flag_compressed;
  //! current values
  int grain_id;
};

/* -------------------------------------------------------------------------- */
template<>
template <typename T>
void DumperLammps<bond>::visitField(T & visited) {

  typename T::iterator it = visited.begin();
  typename T::iterator end = visited.end();

  UInt dim = visited.getDim();

  for (; it != end ; ++it) {
    this->lammps_dump_file << this->curr_nb_atom + 1 << " "
                           << this->grain_id + 2 << " 1 ";
    for (UInt i = 0 ;  i <  dim ; ++i) {
      this->lammps_dump_file << (*it)[i] << " ";
    }
    this->lammps_dump_file << std::endl;
    ++this->curr_nb_atom;
  }
}

/* -------------------------------------------------------------------------- */
template<>
template <typename T>
void DumperLammps<atomic>::visitField(T & visited) {

  typename T::iterator it = visited.begin();
  typename T::iterator end = visited.end();

  UInt dim = visited.getDim();

  for (; it != end ; ++it) {
    this->lammps_dump_file << this->curr_nb_atom + 1 << " 1 ";
    for (UInt i = 0 ;  i <  dim ; ++i) {
      this->lammps_dump_file << (*it)[i] << " ";
    }
    this->lammps_dump_file << std::endl;
    ++this->curr_nb_atom;
  }
}

/* -------------------------------------------------------------------------- */
}

/* -------------------------------------------------------------------------- */
#include "field_inline_impl.hh"
/* -------------------------------------------------------------------------- */



#endif /* IOHELPER_DUMPER_LAMMPS_H_ */
