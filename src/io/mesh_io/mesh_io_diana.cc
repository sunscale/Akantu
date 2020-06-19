/**
 * @file   mesh_io_diana.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Alodie Schneuwly <alodie.schneuwly@epfl.ch>
 *
 * @date creation: Sat Mar 26 2011
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  handles diana meshes
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

/* -------------------------------------------------------------------------- */
#include <fstream>
#include <iostream>

/* -------------------------------------------------------------------------- */
#include "element_group.hh"
#include "mesh_io_diana.hh"
#include "mesh_utils.hh"
/* -------------------------------------------------------------------------- */
#include <string.h>
/* -------------------------------------------------------------------------- */
#include <stdio.h>

namespace akantu {

/* -------------------------------------------------------------------------- */
/*   Methods Implentations                                                    */
/* -------------------------------------------------------------------------- */

MeshIODiana::MeshIODiana() {
  canReadSurface = true;
  canReadExtendedData = true;
  _diana_to_akantu_element_types["T9TM"] = _triangle_3;
  _diana_to_akantu_element_types["CT6CM"] = _triangle_6;
  _diana_to_akantu_element_types["Q12TM"] = _quadrangle_4;
  _diana_to_akantu_element_types["CQ8CM"] = _quadrangle_8;
  _diana_to_akantu_element_types["TP18L"] = _pentahedron_6;
  _diana_to_akantu_element_types["CTP45"] = _pentahedron_15;
  _diana_to_akantu_element_types["TE12L"] = _tetrahedron_4;
  _diana_to_akantu_element_types["HX24L"] = _hexahedron_8;
  _diana_to_akantu_element_types["CHX60"] = _hexahedron_20;
  _diana_to_akantu_mat_prop["YOUNG"] = "E";
  _diana_to_akantu_mat_prop["DENSIT"] = "rho";
  _diana_to_akantu_mat_prop["POISON"] = "nu";

  std::map<std::string, ElementType>::iterator it;
  for (it = _diana_to_akantu_element_types.begin();
       it != _diana_to_akantu_element_types.end(); ++it) {
    UInt nb_nodes = Mesh::getNbNodesPerElement(it->second);

    auto * tmp = new UInt[nb_nodes];
    for (UInt i = 0; i < nb_nodes; ++i) {
      tmp[i] = i;
    }

    switch (it->second) {
    case _tetrahedron_10:
      tmp[8] = 9;
      tmp[9] = 8;
      break;
    case _pentahedron_15:
      tmp[0] = 2;
      tmp[1] = 8;
      tmp[2] = 0;
      tmp[3] = 6;
      tmp[4] = 1;
      tmp[5] = 7;
      tmp[6] = 11;
      tmp[7] = 9;
      tmp[8] = 10;
      tmp[9] = 5;
      tmp[10] = 14;
      tmp[11] = 3;
      tmp[12] = 12;
      tmp[13] = 4;
      tmp[14] = 13;
      break;
    case _hexahedron_20:
      tmp[0] = 5;
      tmp[1] = 16;
      tmp[2] = 4;
      tmp[3] = 19;
      tmp[4] = 7;
      tmp[5] = 18;
      tmp[6] = 6;
      tmp[7] = 17;
      tmp[8] = 13;
      tmp[9] = 12;
      tmp[10] = 15;
      tmp[11] = 14;
      tmp[12] = 1;
      tmp[13] = 8;
      tmp[14] = 0;
      tmp[15] = 11;
      tmp[16] = 3;
      tmp[17] = 10;
      tmp[18] = 2;
      tmp[19] = 9;
      break;
    default:
      // nothing to change
      break;
    }
    _read_order[it->second] = tmp;
  }
}

/* -------------------------------------------------------------------------- */
MeshIODiana::~MeshIODiana() = default;

/* -------------------------------------------------------------------------- */
inline void my_getline(std::ifstream & infile, std::string & line) {
  std::getline(infile, line);   // read the line
  size_t pos = line.find('\r'); /// remove the extra \\r if needed
  line = line.substr(0, pos);
}

/* -------------------------------------------------------------------------- */
void MeshIODiana::read(const std::string & filename, Mesh & mesh) {
  AKANTU_DEBUG_IN();

  MeshAccessor mesh_accessor(mesh);

  std::ifstream infile;
  infile.open(filename.c_str());

  std::string line;
  UInt first_node_number = std::numeric_limits<UInt>::max();
  diana_element_number_to_elements.clear();

  if (!infile.good()) {
    AKANTU_ERROR("Cannot open file " << filename);
  }

  while (infile.good()) {
    my_getline(infile, line);

    /// read all nodes
    if (line == "'COORDINATES'") {
      line = readCoordinates(infile, mesh, first_node_number);
    }

    /// read all elements
    if (line == "'ELEMENTS'") {
      line = readElements(infile, mesh, first_node_number);
    }

    /// read the material properties and write a .dat file
    if (line == "'MATERIALS'") {
      line = readMaterial(infile, filename);
    }

    /// read the material properties and write a .dat file
    if (line == "'GROUPS'") {
      line = readGroups(infile, mesh, first_node_number);
    }
  }
  infile.close();

  mesh_accessor.setNbGlobalNodes(mesh.getNbNodes());

  MeshUtils::fillElementToSubElementsData(mesh);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshIODiana::write(__attribute__((unused)) const std::string & filename,
                        __attribute__((unused)) const Mesh & mesh) {
  AKANTU_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
std::string MeshIODiana::readCoordinates(std::ifstream & infile, Mesh & mesh,
                                         UInt & first_node_number) {
  AKANTU_DEBUG_IN();

  MeshAccessor mesh_accessor(mesh);

  Array<Real> & nodes = mesh_accessor.getNodes();

  std::string line;

  UInt index;
  Vector<Real> coord(3);

  do {
    my_getline(infile, line);
    if ("'ELEMENTS'" == line)
      break;

    std::stringstream sstr_node(line);
    sstr_node >> index >> coord(0) >> coord(1) >> coord(2);

    first_node_number = first_node_number < index ? first_node_number : index;

    nodes.push_back(coord);
  } while (true);

  AKANTU_DEBUG_OUT();
  return line;
}

/* -------------------------------------------------------------------------- */
UInt MeshIODiana::readInterval(std::stringstream & line,
                               std::set<UInt> & interval) {
  UInt first;
  line >> first;
  if (line.fail()) {
    return 0;
  }
  interval.insert(first);

  UInt second;
  int dash;
  dash = line.get();
  if (dash == '-') {
    line >> second;
    interval.insert(second);
    return 2;
  }

  if (line.fail())
    line.clear(std::ios::eofbit); // in case of get at end of the line
  else
    line.unget();
  return 1;
}

/* -------------------------------------------------------------------------- */
std::string MeshIODiana::readGroups(std::ifstream & infile, Mesh & mesh,
                                    UInt first_node_number) {
  AKANTU_DEBUG_IN();

  std::string line;
  my_getline(infile, line);

  bool reading_nodes_group = false;

  while (line != "'SUPPORTS'") {
    if (line == "NODES") {
      reading_nodes_group = true;
      my_getline(infile, line);
    }

    if (line == "ELEMEN") {
      reading_nodes_group = false;
      my_getline(infile, line);
    }

    auto * str = new std::stringstream(line);

    UInt id;
    std::string name;
    char c;
    *str >> id >> name >> c;

    auto * list_ids = new Array<UInt>(0, 1, name);

    UInt s = 1;
    bool end = false;
    while (!end) {
      while (!str->eof() && s != 0) {
        std::set<UInt> interval;
        s = readInterval(*str, interval);
        auto it = interval.begin();
        if (s == 1)
          list_ids->push_back(*it);
        if (s == 2) {
          UInt first = *it;
          ++it;
          UInt second = *it;
          for (UInt i = first; i <= second; ++i) {
            list_ids->push_back(i);
          }
        }
      }
      if (str->fail())
        end = true;
      else {
        my_getline(infile, line);
        delete str;
        str = new std::stringstream(line);
      }
    }

    delete str;

    if (reading_nodes_group) {
      NodeGroup & ng = mesh.createNodeGroup(name);
      for (UInt i = 0; i < list_ids->size(); ++i) {
        UInt node = (*list_ids)(i)-first_node_number;
        ng.add(node, false);
      }
      delete list_ids;

    } else {
      ElementGroup & eg = mesh.createElementGroup(name);
      for (UInt i = 0; i < list_ids->size(); ++i) {
        Element & elem = diana_element_number_to_elements[(*list_ids)(i)];
        if (elem.type != _not_defined)
          eg.add(elem, false, false);
      }

      eg.optimize();
      delete list_ids;
    }

    my_getline(infile, line);
  }

  AKANTU_DEBUG_OUT();
  return line;
}

/* -------------------------------------------------------------------------- */
std::string MeshIODiana::readElements(std::ifstream & infile, Mesh & mesh,
                                      UInt first_node_number) {
  AKANTU_DEBUG_IN();

  std::string line;
  my_getline(infile, line);

  if ("CONNECTIVITY" == line) {
    line = readConnectivity(infile, mesh, first_node_number);
  }

  /// read the line corresponding to the materials
  if ("MATERIALS" == line) {
    line = readMaterialElement(infile, mesh);
  }

  AKANTU_DEBUG_OUT();
  return line;
}

/* -------------------------------------------------------------------------- */
std::string MeshIODiana::readConnectivity(std::ifstream & infile, Mesh & mesh,
                                          UInt first_node_number) {
  AKANTU_DEBUG_IN();

  MeshAccessor mesh_accessor(mesh);
  Int index;
  std::string lline;

  std::string diana_type;
  ElementType akantu_type, akantu_type_old = _not_defined;
  Array<UInt> * connectivity = nullptr;
  UInt node_per_element = 0;
  Element elem;
  UInt * read_order = nullptr;

  while (true) {
    my_getline(infile, lline);
    //    std::cerr << lline << std::endl;
    std::stringstream sstr_elem(lline);
    if (lline == "MATERIALS")
      break;

    /// traiter les coordonnees
    sstr_elem >> index;
    sstr_elem >> diana_type;

    akantu_type = _diana_to_akantu_element_types[diana_type];

    if (akantu_type == _not_defined)
      continue;

    if (akantu_type != akantu_type_old) {
      connectivity = &(mesh_accessor.getConnectivity(akantu_type));

      node_per_element = connectivity->getNbComponent();
      akantu_type_old = akantu_type;
      read_order = _read_order[akantu_type];
    }

    Vector<UInt> local_connect(node_per_element);

    // used if element is written on two lines
    UInt j_last = 0;

    for (UInt j = 0; j < node_per_element; ++j) {
      UInt node_index;
      sstr_elem >> node_index;

      // check s'il y a pas plus rien après un certain point

      if (sstr_elem.fail()) {
        sstr_elem.clear();
        sstr_elem.ignore();
        break;
      }

      node_index -= first_node_number;
      local_connect(read_order[j]) = node_index;
      j_last = j;
    }

    // check if element is written in two lines
    if (j_last != (node_per_element - 1)) {

      // if this is the case, read on more line
      my_getline(infile, lline);
      std::stringstream sstr_elem(lline);

      for (UInt j = (j_last + 1); j < node_per_element; ++j) {

        UInt node_index;
        sstr_elem >> node_index;

        node_index -= first_node_number;
        local_connect(read_order[j]) = node_index;
      }
    }

    connectivity->push_back(local_connect);

    elem.type = akantu_type;
    elem.element = connectivity->size() - 1;

    diana_element_number_to_elements[index] = elem;
    akantu_number_to_diana_number[elem] = index;
  }

  AKANTU_DEBUG_OUT();
  return lline;
}

/* -------------------------------------------------------------------------- */
std::string MeshIODiana::readMaterialElement(std::ifstream & infile,
                                             Mesh & mesh) {
  AKANTU_DEBUG_IN();

  std::string line;

  for (auto type : mesh.elementTypes()) {
    UInt nb_element = mesh.getNbElement(type);
    mesh.getDataPointer<UInt>("material", type, _not_ghost, 1)
        .resize(nb_element);
  }

  my_getline(infile, line);
  while (line != "'MATERIALS'") {
    line =
        line.substr(line.find('/') + 1,
                    std::string::npos); // erase the first slash / of the line
    char tutu[250] = {'\0'};
    strncpy(tutu, line.c_str(), 249);

    AKANTU_DEBUG_WARNING("reading line " << line);
    Array<UInt> temp_id(0, 2);
    UInt mat;
    while (true) {
      std::stringstream sstr_intervals_elements(line);
      Vector<UInt> id(2);
      char temp;
      while (sstr_intervals_elements.good()) {
        sstr_intervals_elements >> id(0) >> temp >> id(1); // >> "/" >> mat;
        if (!sstr_intervals_elements.fail())
          temp_id.push_back(id);
      }
      if (sstr_intervals_elements.fail()) {
        sstr_intervals_elements.clear();
        sstr_intervals_elements.ignore();
        sstr_intervals_elements >> mat;
        break;
      }
      my_getline(infile, line);
    }

    // loop over elements
    //    UInt * temp_id_val = temp_id.storage();
    for (UInt i = 0; i < temp_id.size(); ++i)
      for (UInt j = temp_id(i, 0); j <= temp_id(i, 1); ++j) {
        Element & element = diana_element_number_to_elements[j];
        if (element.type == _not_defined)
          continue;
        UInt elem = element.element;
        ElementType type = element.type;
        Array<UInt> & data =
            mesh.getDataPointer<UInt>("material", type, _not_ghost);
        data(elem) = mat;
      }

    my_getline(infile, line);
  }

  AKANTU_DEBUG_OUT();
  return line;
}

/* -------------------------------------------------------------------------- */
std::string MeshIODiana::readMaterial(std::ifstream & infile,
                                      const std::string & filename) {
  AKANTU_DEBUG_IN();

  std::stringstream mat_file_name;

  mat_file_name << "material_" << filename;

  std::ofstream material_file;
  material_file.open(mat_file_name.str().c_str()); // mat_file_name.str());

  UInt mat_index;
  std::string line;

  bool first_mat = true;
  bool end = false;

  UInt mat_id = 0;

  using MatProp = std::map<std::string, Real>;
  MatProp mat_prop;
  do {
    my_getline(infile, line);
    std::stringstream sstr_material(line);
    if (("'GROUPS'" == line) || ("'END'" == line)) {
      if (!mat_prop.empty()) {
        material_file << "material elastic [" << std::endl;
        material_file << "\tname = material" << ++mat_id << std::endl;
        for (auto it = mat_prop.begin(); it != mat_prop.end(); ++it)
          material_file << "\t" << it->first << " = " << it->second
                        << std::endl;
        material_file << "]" << std::endl;
        mat_prop.clear();
      }
      end = true;
    } else {
      /// traiter les caractéristiques des matériaux
      sstr_material >> mat_index;

      if (!sstr_material.fail()) {
        if (!first_mat) {
          if (!mat_prop.empty()) {
            material_file << "material elastic [" << std::endl;
            material_file << "\tname = material" << ++mat_id << std::endl;
            for (auto it = mat_prop.begin(); it != mat_prop.end(); ++it)
              material_file << "\t" << it->first << " = " << it->second
                            << std::endl;
            material_file << "]" << std::endl;
            mat_prop.clear();
          }
        }
        first_mat = false;
      } else {
        sstr_material.clear();
      }

      std::string prop_name;
      sstr_material >> prop_name;

      std::map<std::string, std::string>::iterator it;
      it = _diana_to_akantu_mat_prop.find(prop_name);

      if (it != _diana_to_akantu_mat_prop.end()) {
        Real value;
        sstr_material >> value;
        mat_prop[it->second] = value;
      } else {
        AKANTU_DEBUG_INFO("In material reader, property " << prop_name
                                                          << "not recognized");
      }
    }
  } while (!end);

  AKANTU_DEBUG_OUT();
  return line;
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
