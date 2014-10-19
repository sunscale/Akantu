/**
 * @file   mesh_io.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Jul 15 2010
 * @date last modification: Fri Jun 13 2014
 *
 * @brief  common part for all mesh io classes
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_common.hh"
#include "mesh_io.hh"

/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
MeshIO::MeshIO() {
  canReadSurface      = false;
  canReadExtendedData = false;
}

/* -------------------------------------------------------------------------- */
MeshIO::~MeshIO() {

}

/* -------------------------------------------------------------------------- */
MeshIO * MeshIO::getMeshIO(const std::string & filename, const MeshIOType & type) {
  MeshIOType t = type;
  if(type == _miot_auto) {
    std::string::size_type idx = filename.rfind('.');
    std::string ext;
    if(idx != std::string::npos) {
      ext = filename.substr(idx+1);
    }

    if(ext == "msh") { t = _miot_gmsh;
    } else if(ext == "diana") { t = _miot_diana;
    } else if(ext == "inp")   { t = _miot_abaqus;
    } else AKANTU_EXCEPTION("Cannot guess the type of file of "
			    << filename << " (ext "<< ext <<"). "
			    << "Please provide the MeshIOType to the read function");
  }

  switch(t) {
  case _miot_gmsh         : return new MeshIOMSH();
#if defined(AKANTU_STRUCTURAL_MECHANICS)
  case _miot_gmsh_struct  : return new MeshIOMSHStruct();
#endif
  case _miot_diana        : return new MeshIODiana();
  case _miot_abaqus       : return new MeshIOAbaqus();
  default:
    return NULL;
  }
}

/* -------------------------------------------------------------------------- */
void MeshIO::read(const std::string & filename, Mesh & mesh, const MeshIOType & type) {
  MeshIO * mesh_io = getMeshIO(filename, type);
  mesh_io->read(filename, mesh);
  delete mesh_io;
}

/* -------------------------------------------------------------------------- */
void MeshIO::write(const std::string & filename, Mesh & mesh, const MeshIOType & type) {
  MeshIO * mesh_io = getMeshIO(filename, type);
  mesh_io->write(filename, mesh);
  delete mesh_io;
}

__END_AKANTU__
