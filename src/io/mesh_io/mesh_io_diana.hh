/**
 * @file   mesh_io_diana.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Alodie Schneuwly <alodie.schneuwly@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sat Mar 26 2011
 * @date last modification: Mon Aug 19 2013
 *
 * @brief  diana mesh reader description
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

#ifndef __AKANTU_MESH_IO_DIANA_HH__
#define __AKANTU_MESH_IO_DIANA_HH__

/* -------------------------------------------------------------------------- */
#include "mesh_io.hh"

/* -------------------------------------------------------------------------- */
#include <vector>

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class MeshIODiana : public MeshIO {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MeshIODiana();
  virtual ~MeshIODiana();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// read a mesh from the file
  virtual void read(const std::string & filename, Mesh & mesh);

  /// write a mesh to a file
  virtual void write(const std::string & filename, const Mesh & mesh);

private:
  std::string readCoordinates(std::ifstream & infile, Mesh & mesh,
                              UInt & first_node_number);

  std::string readElements(std::ifstream & infile, Mesh & mesh,
                           UInt first_node_number);

  std::string readGroups(std::ifstream & infile, Mesh & mesh, UInt first_node_number);

  std::string readConnectivity(std::ifstream & infile, Mesh & mesh,
                               UInt first_node_number);

  std::string readMaterialElement(std::ifstream & infile, Mesh & mesh);

  std::string readMaterial(std::ifstream & infile,
                           const std::string & filename);

  UInt readInterval(std::stringstream & line, std::set<UInt> & interval);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  std::map<std::string, ElementType> _diana_to_akantu_element_types;
  std::map<std::string, std::string> _diana_to_akantu_mat_prop;

  /// order in witch element as to be read, akantu_node_order = _read_order[diana_node_order]
  std::map<ElementType, UInt *> _read_order;

  std::map<UInt, Element> diana_element_number_to_elements;
  std::map<Element, UInt> akantu_number_to_diana_number;
};

__END_AKANTU__

#endif /* __AKANTU_MESH_IO_DIANA_HH__ */
