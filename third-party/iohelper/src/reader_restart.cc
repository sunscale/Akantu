/**
 * @file   reader_restart.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Thu Mar 11 2010
 * @date last modification: Thu Nov 01 2012
 *
 * @brief  implementation for the restart reader
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * IOHelper is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * IOHelper is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with IOHelper. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include <iomanip>

#include "reader_restart.hh"
#include "file_manager.hh"
/* -------------------------------------------------------------------------- */

namespace iohelper {

void ReaderRestart::Read(){

  stringstream temp;
  File file;
  string line;
  stringstream buffer;
  UInt tmp;

  /* node coordinates */
  temp << prefix << "/" << base_name << "-coordinates-" << std::setfill('0') << std::setw(4) << dump_step << ".restart";
  //  file.open(temp.str(),std::ios_base::in | std::ios_base::binary,flag_compressed);
  file >> tmp;
  position->setNbDof(tmp);
  file >> tmp;
  position->setDim(tmp);
  position->AllocateData();
  double * ptr_pos = position->getData();
  file.read((char*)ptr_pos,position->getNbDof()*position->getDim()*sizeof(double));
  file.close();

  /* nodal data */
  {
    std::map<std::string,Field<double> *>::iterator it = per_node_data.begin();
    std::map<std::string,Field<double> *>::iterator end = per_node_data.end();
    while (it != end){
      temp.str(std::string());
      temp << prefix << "/" << base_name << "-" << (*it).second->getName() << "-" << std::setfill('0') << std::setw(4) << dump_step << ".restart";
      //      file.open(temp.str(),fstream::in | fstream::binary ,flag_compressed);
      Field<double> * ptr_field = (*it).second;
      file >> tmp;
      ptr_field->setNbDof(tmp);
      file >> tmp;
      ptr_field->setDim(tmp);
      if (ptr_field->getNbDof() != position->getNbDof()) 
	DUMP("number of degree of freedom of " << ptr_field->name 
	     << " does not match that one of positions field");
      ptr_field->AllocateData();
      double * ptr = ptr_field->getData();
      file.read((char*)ptr,ptr_field->getNbDof()*ptr_field->getDim()*sizeof(double));
      file.close();
      ++it;
    }
  }
  /* if only collection of points do not dump any element data */
  if (elem_type == POINT_SET) return;

  /* connectivity */
  temp.str(std::string());
  temp << prefix << "/" << base_name << "-connectivity-" << std::setfill('0') << std::setw(4) << dump_step << ".restart";
  //  file.open(temp.str(),fstream::in | fstream::binary,flag_compressed);
  unsigned int offset;

  offset = nb_node_per_elem[elem_type];

  file >> tmp;
  connec->setNbDof(tmp);
  file >> tmp;
  connec->setDim(tmp);
  if (offset != connec->getDim()) 
    FATAL("It appears that you are trying to reload a different kind of element !");
  connec->AllocateData();
  file.read((char*)connec->getData(),connec->getNbDof()*connec->getDim()*sizeof(int));
  file.close();

  /* element data */
  {
      std::map<std::string,Field<double> *>::iterator it = per_element_data.begin();
      std::map<std::string,Field<double> *>::iterator end = per_element_data.end();
      while (it != end){
	temp.str(std::string());
	temp << prefix << "/" << base_name << "-" << (*it).second->getName() << "-" << std::setfill('0') << std::setw(4) << dump_step << ".restart";
	//	file.open(temp.str(),fstream::in | fstream::binary,flag_compressed);
	Field<double> * ptr_field = (*it).second;
	file >> tmp;
	ptr_field->setNbDof(tmp);
	file >> tmp;
	ptr_field->setDim(tmp);
	if (ptr_field->getNbDof() != position->getNbDof()) 
	  DUMP("number of degree of freedom of " << ptr_field->name 
	       << " does not match that one of connectivity field");
	ptr_field->AllocateData();
	double * ptr = ptr_field->getData();
	file.read((char*)ptr,ptr_field->getNbDof()*ptr_field->getDim()*sizeof(double));
	file.close();
	++it;
      }
  }
  ++dump_step;
}



}
