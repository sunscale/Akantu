/**
 * @file   dumper_paraview.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Thu Mar 11 2010
 * @date last modification: Wed Jun 05 2013
 *
 * @brief  header for paraview dumper
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

#ifndef IOHELPER_DUMPER_PARAVIEW_H_
#define IOHELPER_DUMPER_PARAVIEW_H_
/* -------------------------------------------------------------------------- */
#include <map>
#include <string>
#include "dumper.hh"
				    //#include "field.hh"
#include "paraview_helper.hh"
/* -------------------------------------------------------------------------- */

namespace iohelper {

/** Class DumperParaview
 * Implementation of a dumper to paraview vtu files
 */
class DumperParaview : public Dumper {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

public:

  DumperParaview(const std::string & prefix = std::string("./"));
  ~DumperParaview() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

  //! dump to file
  void dump(const std::string & name = std::string(""),
            UInt count = UInt(-1)) override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */

  //! set mode for file creation : TEXT, BASE64, COMPRESSED
  void setMode(Int mode) override;

  //! set the sub directory where to store the vtu files
  void setVTUSubDirectory(const std::string & sub);

  //! push a single field nodal with templated type
  template <typename T> void pushNodeField(ParaviewHelper & paraHelper, Field<T> & f);
  //! push a single element field with templated type
  template <typename T> void pushElemField(ParaviewHelper & paraHelper, Field<T> & f);

  //! give vector to connectivity
  void setConnectivity(int * connectivity, ElemType element_type, UInt nb_elem,
                       int mode) override;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

private:

  //! flag to produce zipped files
  bool flag_compressed;

  //! offset to compute connectivities
  Int mode_offset;

  //! subdirectory where to put the vtu files
  //std::string vtu_subdirectory; //became sub_folder in dumper.hh
  std::vector< std::pair<Real, std::string> > pvtu_file_names;
};


/* -------------------------------------------------------------------------- */

}



#endif /* IOHELPER_DUMPER_PARAVIEW_H_ */
