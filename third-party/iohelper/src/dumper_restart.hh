/**
 * @file   dumper_restart.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Thu Mar 11 2010
 * @date last modification: Tue Jun 04 2013
 *
 * @brief  header for dumper restart
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
#ifndef IOHELPER_DUMPER_RESTART_H_
#define IOHELPER_DUMPER_RESTART_H_
/* -------------------------------------------------------------------------- */
#include <map>
#include <string>
#include "dumper.hh"
				    //#include "field.hh"
#include "file_manager.hh"
/* -------------------------------------------------------------------------- */

namespace iohelper {


/** Class DumperRestart
 * Implementation of a dumper to restart
 */

class DumperRestart : public Dumper, public Visitor {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

public:

  DumperRestart(std::string prefix = "./"):Dumper(prefix){};
  ~DumperRestart(){};


  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

  void dump(const std::string & name);
  void setEmbeddedValue(__attribute__((unused)) const std::string & name,
			__attribute__((unused)) int value) {};

  //! visitor system
  template <typename T> void visit(T & visited);

private:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

private:

};

/* -------------------------------------------------------------------------- */


template <typename T>
void DumperRestart::visit(T & visited){
  File file;

  int nb = visited.size();
  /* node coordinates */
  // std::stringstream filename;
  // filename << this->getFullName(visited.getName());
  // filename << ".restart";

  // file.open(filename.str(),std::fstream::out | std::fstream::binary);
  // file << nb << "\t" << visited.getDim() << std::endl;

  std::string fname = this->getAbsoluteFilePath(visited.getName(),
						"",
						".restart");
  file.open(fname, std::fstream::out | std::fstream::binary);
  file << nb << "\t" << visited.getDim() << std::endl;

  typename T::iterator it = visited.begin();
  typename T::iterator end = visited.end();

  for (; it != end ; ++it) {
    IOHelperVector<T> data = *it;
    T * ptr = data.getPtr();
    UInt nb = data.size();
    file.write((char*)ptr,nb*visited.getDim()*sizeof(T));
  }
  file.close();

}

/* -------------------------------------------------------------------------- */


}


#endif /* IOHELPER_DUMPER_RESTART_H_ */
