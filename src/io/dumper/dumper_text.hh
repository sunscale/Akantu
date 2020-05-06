/**
 * @file   dumper_text.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  to dump into a text file
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "dumper_iohelper.hh"
/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_DUMPER_TEXT_HH__
#define __AKANTU_DUMPER_TEXT_HH__
/* -------------------------------------------------------------------------- */
#include <io_helper.hh>
/* -------------------------------------------------------------------------- */

namespace akantu {

class DumperText : public DumperIOHelper {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  DumperText(const std::string & basename = "dumper_text",
             iohelper::TextDumpMode mode = iohelper::_tdm_space,
             bool parallel = true);
  ~DumperText() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void
  registerMesh(const Mesh & mesh, UInt spatial_dimension = _all_dimensions,
               const GhostType & ghost_type = _not_ghost,
               const ElementKind & element_kind = _ek_not_defined) override;

  void registerFilteredMesh(
      const Mesh & mesh, const ElementTypeMapArray<UInt> & elements_filter,
      const Array<UInt> & nodes_filter,
      UInt spatial_dimension = _all_dimensions,
      const GhostType & ghost_type = _not_ghost,
      const ElementKind & element_kind = _ek_not_defined) override;

  void setBaseName(const std::string & basename) override;

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

} // namespace akantu

#endif /* __AKANTU_DUMPER_TEXT_HH__ */
