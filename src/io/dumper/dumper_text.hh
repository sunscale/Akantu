/**
 * @file   dumper_text.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Fri May 17 2013
 * @date last modification: Fri Sep 05 2014
 *
 * @brief  to dump into a text file
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "dumper_iohelper.hh"
/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_DUMPER_TEXT_HH__
#define __AKANTU_DUMPER_TEXT_HH__
/* -------------------------------------------------------------------------- */
#include <io_helper.hh>
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class DumperText : public DumperIOHelper {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  DumperText(const std::string & basename = "dumper_text",
	     iohelper::TextDumpMode mode = iohelper::_tdm_space,
	     bool parallel = true);
  virtual ~DumperText() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void registerMesh(const Mesh & mesh, UInt spatial_dimension = _all_dimensions,
			    const GhostType & ghost_type = _not_ghost,
			    const ElementKind & element_kind = _ek_not_defined);

  virtual void registerFilteredMesh(const Mesh & mesh,
				    const ElementTypeMapArray<UInt> & elements_filter,
				    const Array<UInt> & nodes_filter,
				    UInt spatial_dimension = _all_dimensions,
				    const GhostType & ghost_type = _not_ghost,
				    const ElementKind & element_kind = _ek_not_defined);

  virtual void setBaseName(const std::string & basename);

private:
  void registerNodeTypeField();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  void setPrecision(UInt prec);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
};

__END_AKANTU__

#endif /* __AKANTU_DUMPER_TEXT_HH__ */
