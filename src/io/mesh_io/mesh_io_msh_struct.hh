/**
 * @file   mesh_io_msh_struct.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Fri Jan 26 2018
 *
 * @brief  Read/Write for MSH files
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
#ifndef AKANTU_MESH_IO_MSH_STRUCT_HH_
#define AKANTU_MESH_IO_MSH_STRUCT_HH_
/* -------------------------------------------------------------------------- */
#include "mesh_io.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

class MeshIOMSHStruct : public MeshIOMSH {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MeshIOMSHStruct();

  /// read a mesh from the file
  void read(const std::string & filename, Mesh & mesh) override;
};

} // namespace akantu

#endif /* AKANTU_MESH_IO_MSH_STRUCT_HH_ */
