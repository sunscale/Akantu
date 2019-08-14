/**
 * @file   mesh_io.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Thu Feb 01 2018
 *
 * @brief  common part for all mesh io classes
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
#include "mesh_io.hh"
#include "aka_common.hh"
#include "aka_iterators.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
MeshIO::MeshIO() {
  canReadSurface = false;
  canReadExtendedData = false;
}

/* -------------------------------------------------------------------------- */
MeshIO::~MeshIO() = default;

/* -------------------------------------------------------------------------- */
std::unique_ptr<MeshIO> MeshIO::getMeshIO(const std::string & filename,
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
    } else
      AKANTU_EXCEPTION("Cannot guess the type of file of "
                       << filename << " (ext " << ext << "). "
                       << "Please provide the MeshIOType to the read function");
  }

  switch (t) {
  case _miot_gmsh:
    return std::make_unique<MeshIOMSH>();
#if defined(AKANTU_STRUCTURAL_MECHANICS)
  case _miot_gmsh_struct:
    return std::make_unique<MeshIOMSHStruct>();
#endif
  case _miot_diana:
    return std::make_unique<MeshIODiana>();
  default:
    return nullptr;
  }
}

/* -------------------------------------------------------------------------- */
void MeshIO::read(const std::string & filename, Mesh & mesh,
                  const MeshIOType & type) {
  std::unique_ptr<MeshIO> mesh_io = getMeshIO(filename, type);
  mesh_io->read(filename, mesh);
}

/* -------------------------------------------------------------------------- */
void MeshIO::write(const std::string & filename, Mesh & mesh,
                   const MeshIOType & type) {
  std::unique_ptr<MeshIO> mesh_io = getMeshIO(filename, type);
  mesh_io->write(filename, mesh);
}

/* -------------------------------------------------------------------------- */
void MeshIO::constructPhysicalNames(const std::string & tag_name, Mesh & mesh) {
  if (!physical_names.empty()) {
    for (auto type : mesh.elementTypes()) {
      auto & name_vec =
          mesh.getDataPointer<std::string>("physical_names", type);

      const auto & tags_vec = mesh.getData<UInt>(tag_name, type);

      for (auto pair : zip(tags_vec, name_vec)) {
        auto tag = std::get<0>(pair);
        auto & name = std::get<1>(pair);
        auto map_it = physical_names.find(tag);

        if (map_it == physical_names.end()) {
          std::stringstream sstm;
          sstm << tag;

          name = sstm.str();
        } else {
          name = map_it->second;
        }
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void MeshIO::printself(std::ostream & stream, int indent) const {
  std::string space(AKANTU_INDENT, indent);

  if (physical_names.size()) {
    stream << space << "Physical map:" << std::endl;
    for (auto & pair : physical_names) {
      stream << space << pair.first << ": " << pair.second << std::endl;
    }
  }
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
