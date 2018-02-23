/**
 * @file   mesh_io_abaqus.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sun Sep 26 2010
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  read a mesh from an abaqus input file
 *
 * @section LICENSE
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
#include "mesh_accessor.hh"
#include "mesh_io.hh"

#ifndef __AKANTU_MESH_IO_ABAQUS_HH__
#define __AKANTU_MESH_IO_ABAQUS_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
class MeshIOAbaqus : public MeshIO {
public:
  MeshIOAbaqus();
  ~MeshIOAbaqus() override;

  /// read a mesh from the file
  void read(const std::string & filename, Mesh & mesh) override;

  /// write a mesh to a file
  //  virtual void write(const std::string & filename, const Mesh & mesh);

private:
  /// correspondence between msh element types and akantu element types
  std::map<std::string, ElementType> _abaqus_to_akantu_element_types;
};

} // akantu

#endif /* __AKANTU_MESH_IO_ABAQUS_HH__ */
