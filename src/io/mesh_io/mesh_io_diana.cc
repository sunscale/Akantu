/**
 * @file   mesh_io_diana.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Alodie Schneuwly <alodie.schneuwly@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sat Mar 26 2011
 * @date last modification: Thu Mar 27 2014
 *
 * @brief  handles diana meshes
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


/* -------------------------------------------------------------------------- */
#include <fstream>
#include <iostream>

/* -------------------------------------------------------------------------- */
#include "mesh_io_diana.hh"

/* -------------------------------------------------------------------------- */
#include <string.h>
/* -------------------------------------------------------------------------- */
#include <stdio.h>

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/*   Methods Implentations                                                    */
/* -------------------------------------------------------------------------- */

MeshIODiana::MeshIODiana() {
  canReadSurface      = true;
  canReadExtendedData = true;
  _diana_to_akantu_element_types["T9TM"] = _triangle_3;
  _diana_to_akantu_element_types["Q12TM"] = _quadrangle_4;
  _diana_to_akantu_element_types["TP18L"] = _pentahedron_6;
  _diana_to_akantu_element_types["TE12L"] = _tetrahedron_4;
  _diana_to_akantu_element_types["HX24L"] = _hexahedron_8;
  _diana_to_akantu_mat_prop["YOUNG"] = "E";
  _diana_to_akantu_mat_prop["DENSIT"] = "rho";
  _diana_to_akantu_mat_prop["POISON"] = "nu";
}

/* -------------------------------------------------------------------------- */
MeshIODiana::~MeshIODiana() {
  std::map<std::string, Array<UInt> *>::iterator ng_it;
  std::map<std::string, std::vector<Element> *>::iterator eg_it;

  for (ng_it = node_groups.begin(); ng_it != node_groups.end(); ++ng_it) {
    delete ng_it->second;
  }

  for (eg_it = element_groups.begin(); eg_it != element_groups.end(); ++eg_it) {
    delete eg_it->second;
  }

}

/* -------------------------------------------------------------------------- */
inline void my_getline(std::ifstream & infile, std::string & line) {
  std::getline(infile, line); //read the line
  size_t pos = line.find("\r"); /// remove the extra \r if needed
  line = line.substr(0, pos);
}


/* -------------------------------------------------------------------------- */
void MeshIODiana::read(const std::string & filename, Mesh & mesh) {
  AKANTU_DEBUG_IN();

  std::ifstream infile;
  infile.open(filename.c_str());

  std::string line;
  UInt first_node_number = std::numeric_limits<UInt>::max();
  std::vector<Element> global_to_local_index;

  if(!infile.good()) {
    AKANTU_DEBUG_ERROR("Cannot open file " << filename);
  }

  while(infile.good()) {
    my_getline(infile, line);

    /// read all nodes
    if(line == "'COORDINATES'") {
      line = readCoordinates(infile, mesh, first_node_number);
    }

    /// read all elements
    if (line == "'ELEMENTS'") {
      line = readElements(infile, mesh, global_to_local_index, first_node_number);
    }

    /// read the material properties and write a .dat file
    if (line == "'MATERIALS'") {
      line = readMaterial(infile, filename);
    }

    /// read the material properties and write a .dat file
    if (line == "'GROUPS'") {
      line = readGroups(infile, global_to_local_index, first_node_number);
    }

  }
  infile.close();

  mesh.nb_global_nodes = mesh.nodes->getSize();

  mesh.registerEventHandler(*this);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshIODiana::write(__attribute__((unused)) const std::string & filename,
			__attribute__((unused)) const Mesh & mesh) {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
std::string MeshIODiana::readCoordinates(std::ifstream & infile, Mesh & mesh, UInt & first_node_number) {
  AKANTU_DEBUG_IN();

  Array<Real> & nodes = const_cast<Array<Real> &>(mesh.getNodes());

  std::string line;

  UInt index;
  Real coord[3];

  do {
    my_getline(infile, line);
    if("'ELEMENTS'" == line)
      break;
    //end = true;
    //else {
    /// for each node, read the coordinates

    std::stringstream sstr_node(line);

    sstr_node >> index >> coord[0] >> coord[1] >> coord[2];

    //if (!sstr_node.fail())
    //break;

    first_node_number = first_node_number < index ? first_node_number : index;

    nodes.push_back(coord);
    // }
  } while(true);//!end);

  AKANTU_DEBUG_OUT();
  return line;
}

/* -------------------------------------------------------------------------- */
UInt MeshIODiana::readInterval(std::stringstream & line,
			       std::set<UInt> & interval) {
  UInt first;
  line >> first;
  if(line.fail()) { return 0; }
  interval.insert(first);

  UInt second;
  int dash;
  dash = line.get();
  if(dash == '-') {
    line >> second;
    interval.insert(second);
    return 2;
  }

  if(line.fail())
    line.clear(std::ios::eofbit);  // in case of get at end of the line
  else line.unget();
  return 1;
}

/* -------------------------------------------------------------------------- */
std::string MeshIODiana::readGroups(std::ifstream & infile,
				      std::vector<Element> & global_to_local_index,
				      UInt first_node_number) {
  AKANTU_DEBUG_IN();

  std::string line;
  my_getline(infile, line);

  bool reading_nodes_group = false;

  while(line != "'SUPPORTS'") {
    if(line == "NODES") {
      reading_nodes_group   = true;
      my_getline(infile, line);
    }

    if(line == "ELEMEN") {
      reading_nodes_group   = false;
      my_getline(infile, line);
    }

    std::stringstream *str = new std::stringstream(line);

    UInt id;
    std::string name;
    char c;
    *str >> id >> name >> c;

    Array<UInt> * list_ids = new Array<UInt>(0, 1, name);

    UInt s = 1; bool end = false;
    while(!end) {
      while(!str->eof() && s != 0) {
	std::set<UInt> interval;
	s = readInterval(*str, interval);
	std::set<UInt>::iterator it = interval.begin();
	if(s == 1) list_ids->push_back(*it);
	if(s == 2) {
	  UInt first = *it;
	  ++it;
	  UInt second = *it;
	  for(UInt i = first; i <= second; ++i) {
	    list_ids->push_back(i);
	  }
	}
      }
      if(str->fail()) end = true;
      else {
	my_getline(infile, line);
	delete str;
	str = new std::stringstream(line);
      }
    }

    delete str;

    if(reading_nodes_group) {
      for (UInt i = 0; i < list_ids->getSize(); ++i) {
	(*list_ids)(i) -= first_node_number;
      }
      node_groups[name] = list_ids;
    } else {
      std::vector<Element> * elem = new std::vector<Element>;
      elem->reserve(list_ids->getSize());
      for (UInt i = 0; i < list_ids->getSize(); ++i) {
	Element & e = global_to_local_index[(*list_ids)(i)-1];
	if(e.type != _not_defined)
	  elem->push_back(e);
      }

      element_groups[name] = elem;
      delete list_ids;
    }

    my_getline(infile, line);
  }

  AKANTU_DEBUG_OUT();
  return line;
}

/* -------------------------------------------------------------------------- */
std::string MeshIODiana::readElements(std::ifstream & infile,
				      Mesh & mesh,
				      std::vector<Element> & global_to_local_index,
				      UInt first_node_number) {
  AKANTU_DEBUG_IN();

  std::string line;
  my_getline(infile, line);

  if("CONNECTIVITY" == line) {
    line = readConnectivity(infile, mesh, global_to_local_index, first_node_number);
  }

  /// read the line corresponding to the materials
  if ("MATERIALS" == line) {
    line = readMaterialElement(infile, mesh, global_to_local_index);
  }

  AKANTU_DEBUG_OUT();
  return line;
}


/* -------------------------------------------------------------------------- */
std::string MeshIODiana::readConnectivity(std::ifstream & infile,
					  Mesh & mesh,
					  std::vector<Element> & global_to_local_index,
					  UInt first_node_number) {
  AKANTU_DEBUG_IN();

  Int index;
  std::string lline;

  std::string diana_type;
  ElementType akantu_type, akantu_type_old = _not_defined;
  Array<UInt> *connectivity = NULL;
  UInt node_per_element = 0;
  Element elem;

  bool end = false;
  do {
    my_getline(infile, lline);
    std::stringstream sstr_elem(lline);
    if(lline == "MATERIALS") end = true;
    else {
      /// traiter les coordonnees
      sstr_elem >> index;
      sstr_elem >> diana_type;

      akantu_type = _diana_to_akantu_element_types[diana_type];

      if(akantu_type != _not_defined) {
	if(akantu_type != akantu_type_old) {
	  connectivity = mesh.getConnectivityPointer(akantu_type);

	  node_per_element = connectivity->getNbComponent();
	  akantu_type_old = akantu_type;
	}

	UInt local_connect[node_per_element];
	for(UInt j = 0; j < node_per_element; ++j) {
	  UInt node_index;
	  sstr_elem >> node_index;

	  node_index -= first_node_number;
	  local_connect[j] = node_index;
	}
	connectivity->push_back(local_connect);

	elem.type = akantu_type;
	elem.element = connectivity->getSize() - 1;
      } else {
	elem.type = _not_defined;
	elem.element = UInt(-1);
      }

      global_to_local_index.push_back(elem);
    }
  }
  while(!end);

  AKANTU_DEBUG_OUT();
  return lline;
}

/* -------------------------------------------------------------------------- */
std::string MeshIODiana::readMaterialElement(std::ifstream & infile,
					     Mesh & mesh,
					     std::vector<Element> & global_to_local_index) {
  AKANTU_DEBUG_IN();


  std::string line;
  std::stringstream sstr_tag_name; sstr_tag_name << "tag_" << 0;

  Mesh::type_iterator it  = mesh.firstType();
  Mesh::type_iterator end = mesh.lastType();
  for(; it != end; ++it) {
    UInt nb_element = mesh.getNbElement(*it);
    mesh.getDataPointer<UInt>("material", *it, _not_ghost, 1)->resize(nb_element);
  }

  my_getline(infile, line);
  while(line != "'MATERIALS'") {
    line = line.substr(line.find('/') + 1, std::string::npos); // erase the first slash / of the line
    char tutu[250];
    strcpy(tutu, line.c_str());
    Array<UInt> temp_id(0, 2);
    UInt mat;
    while(true){
      std::stringstream sstr_intervals_elements(line);
      UInt id[2];
      char temp;
      while(sstr_intervals_elements.good()) {
	sstr_intervals_elements >> id[0] >> temp >> id[1]; // >> "/" >> mat;
	if(!sstr_intervals_elements.fail())
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
    for (UInt i = 0; i < temp_id.getSize(); ++i)
      for (UInt j = temp_id(i, 0); j <= temp_id(i, 1); ++j) {
	Element & element = global_to_local_index[j - 1];
	if(element.type == _not_defined) continue;
	UInt elem = element.element;
	ElementType type = element.type;
	Array<UInt> & data = *(mesh.getDataPointer<UInt>("material", type, _not_ghost));
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
  material_file.open(mat_file_name.str().c_str());//mat_file_name.str());

  UInt mat_index;
  std::string line;

  bool first_mat = true;
  bool end = false;

  UInt mat_id = 0;

  typedef std::map<std::string, Real> MatProp;
  MatProp mat_prop;
  do{
    my_getline(infile, line);
    std::stringstream sstr_material(line);
    if("'GROUPS'" == line) {
      if(!mat_prop.empty()) {
	material_file << "material elastic [" << std::endl;
	material_file << "\tname = material" << ++mat_id << std::endl;
	for(MatProp::iterator it = mat_prop.begin();
	    it != mat_prop.end(); ++it)
	  material_file << "\t" << it->first << " = " << it->second << std::endl;
	material_file << "]" << std::endl;
	mat_prop.clear();
      }
      end = true;
    }
    else {
      /// traiter les caractéristiques des matériaux
      sstr_material >> mat_index;

      if(!sstr_material.fail()) {
	if(!first_mat) {
	  if(!mat_prop.empty()) {
	    material_file << "material elastic [" << std::endl;
	    material_file << "\tname = material" << ++mat_id << std::endl;
	    for(MatProp::iterator it = mat_prop.begin();
		it != mat_prop.end(); ++it)
	      material_file << "\t" << it->first << " = " << it->second << std::endl;
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

      if(it != _diana_to_akantu_mat_prop.end()) {
	Real value;
	sstr_material >> value;
	mat_prop[it->second] = value;
      } else {
	AKANTU_DEBUG_INFO("In material reader, property " << prop_name << "not recognized");
      }
    }
  } while (!end);

  AKANTU_DEBUG_OUT();
  return line;
}


void MeshIODiana::onNodesRemoved(const Array<UInt> & element_list,
				 const Array<UInt> & new_numbering,
				 const RemovedNodesEvent & event) {
  std::map<std::string, Array<UInt> *>::iterator it = node_groups.begin();
  std::map<std::string, Array<UInt> *>::iterator end = node_groups.end();
  for (; it != end; ++it) {
    Array<UInt> & group = *(it->second);
    Array<UInt>::iterator<> git = group.begin();
    Array<UInt>::iterator<> gend = group.end();
    for (; git != gend; ++git) {
      UInt new_id = new_numbering(*git);
      AKANTU_DEBUG_ASSERT(new_id != UInt(-1), "Argh " << *git << " was suppressed!");
      *git = new_id;
    }
  }
}


__END_AKANTU__
