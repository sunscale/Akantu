/**
 * @file   mesh_io.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Aug 09 2017
 *
 * @brief  interface of a mesh io class, reader and writer
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
#ifndef AKANTU_MESH_IO_HH_
#define AKANTU_MESH_IO_HH_

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_accessor.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

class MeshIO {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MeshIO();

  virtual ~MeshIO();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  static void read(const std::string & filename, Mesh & mesh,
                   const MeshIOType & type);
  static void write(const std::string & filename, Mesh & mesh,
                    const MeshIOType & type);

  /// read a mesh from the file
  virtual void read(const std::string & /*filename*/, Mesh & /*mesh*/) {}

  /// write a mesh to a file
  virtual void write(const std::string & /*filename*/, const Mesh & /*mesh*/) {}

  /// function to request the manual construction of the physical names maps
  virtual void constructPhysicalNames(const std::string & tag_name,
                                      Mesh & mesh);

  /// method to permit to be printed to a generic stream
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /// static contruction of a meshio object
  static std::unique_ptr<MeshIO> getMeshIO(const std::string & filename,
                                           const MeshIOType & type);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  auto & getPhysicalNames() { return this->physical_names; }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  bool canReadSurface{false};
  bool canReadExtendedData{false};

  /// correspondance between a tag and physical names (if applicable)
  std::map<int, std::string> physical_names;
};

/* -------------------------------------------------------------------------- */

inline std::ostream & operator<<(std::ostream & stream, const MeshIO & _this) {
  _this.printself(stream);
  return stream;
}

/* -------------------------------------------------------------------------- */

} // namespace akantu

#include "mesh_io_diana.hh"
#include "mesh_io_msh.hh"

#if defined(AKANTU_STRUCTURAL_MECHANICS)
#include "mesh_io_msh_struct.hh"
#endif

#endif /* AKANTU_MESH_IO_HH_ */
