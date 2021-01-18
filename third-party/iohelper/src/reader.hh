/**
 * @file   reader.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Thu Mar 11 2010
 * @date last modification: Thu Nov 01 2012
 *
 * @brief  reader description
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * IOHelper is free  software: you can redistribute it and/or  modify it under
 * the terms  of the  GNU Lesser  General Public  License as  published by  the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * IOHelper is  distributed in the  hope that it  will be useful, but  WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for
 * more details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with IOHelper. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef IOHELPER_READER_H_
#define IOHELPER_READER_H_
/* -------------------------------------------------------------------------- */
#include "iohelper_common.hh"
#include <map>
#include <string>
/* -------------------------------------------------------------------------- */

namespace iohelper {

/** Class Reader
 * Interface of a reader
 */

class Reader {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

public:
  Reader();
  virtual ~Reader();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

  //! dump to file
  virtual void Read() = 0;
  //! do the allocation
  void Init();
  //! give vector with coordinates
  void SetPoints(const std::string & name);
  //! give vector with coordinates
  void SetConnectivity(int elem_type);
  //! give vector to per node data
  void AddNodeDataField(const std::string & name);
  //! give vector to per element data
  void AddElemDataField(const char * name);
  //! set prefix directory
  void SetPrefix(const std::string & dir);
  // ! set parallel context
  void SetParallelContext(int me, int wld_size);
  //! set mode
  virtual void SetMode(int mode) { flag_compressed = mode & COMPRESSED; }

  //! give vector with coordinates
  double * GetPoints();
  //! give vector to connectivity
  int * GetConnectivity();
  //! give vector to per node data
  double * GetNodeDataField(const char * name);
  //! give vector to per element data
  double * GetElemDataField(const char * name);

  //! return number of read nodes
  int GetNumberNodes();
  //! return number of read elements
  int GetNumberElements();

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

protected:
  std::string prefix;
  std::string base_name;
  int dump_step;

  //! flag to produce zipped files
  bool flag_compressed;

  // //! vector of positions
  // Field<double> * position;
  // //! vector of connectivity
  // Field<int> * connec;
  // //! vector of additional per node data
  // std::map<std::string,Field<double> *> per_node_data;
  // //! vector of additional per element data
  // std::map<std::string,Field<double> *> per_element_data;

  int world_size;
  int my_rank;
  int elem_type;
  int connectivity_mode;
};

} // namespace iohelper

#endif /* IOHELPER_READER_H_ */
