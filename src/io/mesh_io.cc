/**
 * @file   mesh_io.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Mon Jun 01 2015
 *
 * @brief  common part for all mesh io classes
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "aka_common.hh"
#include "mesh_io.hh"

/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
MeshIO::MeshIO() {
  canReadSurface = false;
  canReadExtendedData = false;
}

/* -------------------------------------------------------------------------- */
MeshIO::~MeshIO() {}

/* -------------------------------------------------------------------------- */
MeshIO * MeshIO::getMeshIO(const std::string & filename,
                           const MeshIOType & type) {
  MeshIOType t = type;
  if (type == _miot_auto) {
    std::string::size_type idx = filename.rfind('.');
    std::string ext;
    if (idx != std::string::npos) {
      ext = filename.substr(idx + 1);
    }

    if (ext == "msh") {
      t = _miot_gmsh;
    } else if (ext == "diana") {
      t = _miot_diana;
    } else if (ext == "inp") {
      t = _miot_abaqus;
    } else
      AKANTU_EXCEPTION("Cannot guess the type of file of "
                       << filename << " (ext " << ext << "). "
                       << "Please provide the MeshIOType to the read function");
  }

  switch (t) {
  case _miot_gmsh:
    return new MeshIOMSH();
#if defined(AKANTU_STRUCTURAL_MECHANICS)
  case _miot_gmsh_struct:
    return new MeshIOMSHStruct();
#endif
  case _miot_diana:
    return new MeshIODiana();
  case _miot_abaqus:
    return new MeshIOAbaqus();
  default:
    return NULL;
  }
}

/* -------------------------------------------------------------------------- */
void MeshIO::read(const std::string & filename, Mesh & mesh,
                  const MeshIOType & type) {
  MeshIO * mesh_io = getMeshIO(filename, type);
  mesh_io->read(filename, mesh);
  delete mesh_io;
}

/* -------------------------------------------------------------------------- */
void MeshIO::write(const std::string & filename, Mesh & mesh,
                   const MeshIOType & type) {
  MeshIO * mesh_io = getMeshIO(filename, type);
  mesh_io->write(filename, mesh);
  delete mesh_io;
}

/* -------------------------------------------------------------------------- */
void MeshIO::constructPhysicalNames(const std::string & tag_name, Mesh & mesh) {

  if (!phys_name_map.empty()) {
    for (Mesh::type_iterator type_it = mesh.firstType();
         type_it != mesh.lastType(); ++type_it) {

      Array<std::string> * name_vec =
          mesh.getDataPointer<std::string>("physical_names", *type_it);

      const Array<UInt> & tags_vec = mesh.getData<UInt>(tag_name, *type_it);

      for (UInt i(0); i < tags_vec.getSize(); i++) {
        std::map<UInt, std::string>::const_iterator map_it =
            phys_name_map.find(tags_vec(i));

        if (map_it == phys_name_map.end()) {
          std::stringstream sstm;
          sstm << tags_vec(i);
          name_vec->operator()(i) = sstm.str();
        } else {
          name_vec->operator()(i) = map_it->second;
        }
      }
    }
  }
}

/* -------------------------------------------------------------------------- */

void MeshIO::printself(std::ostream & stream, int indent) const {

  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
    ;

  if (phys_name_map.size()) {

    stream << space << "Physical map:" << std::endl;

    std::map<UInt, std::string>::const_iterator it = phys_name_map.begin();
    std::map<UInt, std::string>::const_iterator end = phys_name_map.end();

    for (; it != end; ++it) {
      stream << space << it->first << ": " << it->second << std::endl;
    }
  }
}

/* -------------------------------------------------------------------------- */

} // akantu
