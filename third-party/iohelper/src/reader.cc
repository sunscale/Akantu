/**
 * @file   reader.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Thu Mar 11 2010
 * @date last modification: Thu Nov 01 2012
 *
 * @brief  reader implementation
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

#include "reader.hh"
/* -------------------------------------------------------------------------- */
namespace iohelper {
/* -------------------------------------------------------------------------- */


Reader::Reader(){
  position = NULL;
  connec = NULL;
  dump_step = 0;
  elem_type = -1;
  base_name[0] = '\0';
  prefix[0] = '\0';
}
/* -------------------------------------------------------------------------- */


Reader::~Reader(){
  if (position)
    delete position;
  if (connec)
    delete connec;

  {
    std::map<std::string,Field<double> *>::iterator it = per_node_data.begin();
    std::map<std::string,Field<double> *>::iterator end = per_node_data.end();
    while (it != end){
      delete (*it).second;
      ++it;
    }
  }
  {
    std::map<std::string,Field<double> *>::iterator it = per_element_data.begin();
    std::map<std::string,Field<double> *>::iterator end = per_element_data.end();
    while (it != end){
      delete (*it).second;
      ++it;
    }
  }
}
/* -------------------------------------------------------------------------- */


void Reader::Init(){
}
/* -------------------------------------------------------------------------- */


void Reader::AddNodeDataField(const string & name){
  if (per_node_data.count(name) != 0) delete(per_node_data[name]);
  Field<double> * temp = new Field<double>();
  temp->setName(name);
  per_node_data[name] = temp;
}
/* -------------------------------------------------------------------------- */


void Reader::AddElemDataField(const char * name){
  if (per_element_data.count(name) != 0) delete(per_element_data[name]);
  if (connec == NULL) FATAL("connectivity should be provided before elemental fields ! Please use SetConnectivity function before AddElemDataField");
  Field<double> * temp = new Field<double>();
  temp->setName(name);
  per_element_data[name] = temp;
}
/* -------------------------------------------------------------------------- */


void Reader::SetPoints(const std::string & n){
  position = new Field<double>();
  base_name = n;
}
/* -------------------------------------------------------------------------- */


void Reader::SetConnectivity(int elem_type){
  connec = new Field<int>();
  this->elem_type = elem_type;
}
/* -------------------------------------------------------------------------- */


void Reader::SetPrefix(const std::string & dir){
  prefix = dir;
}
/* -------------------------------------------------------------------------- */


void Reader::SetParallelContext(int me, int wld_size){
  my_rank = me;
  world_size = wld_size;
}
/* -------------------------------------------------------------------------- */


double * Reader::GetPoints(){
  return position->getData();
}
/* -------------------------------------------------------------------------- */


int * Reader::GetConnectivity(){
  return connec->getData();
}
/* -------------------------------------------------------------------------- */


double * Reader::GetNodeDataField(const char * name){
  if (per_node_data.count(name) == 0) FATAL("node data field named " << name << " was not reloaded");
  return per_node_data[name]->getData();
}
/* -------------------------------------------------------------------------- */


double * Reader::GetElemDataField(const char * name){
  if (per_element_data.count(name) == 0) FATAL("elem data field named " << name << " was not reloaded");
  return per_element_data[name]->getData();
}
/* -------------------------------------------------------------------------- */


int Reader::GetNumberNodes(){
  return position->getNbDof();
}
/* -------------------------------------------------------------------------- */


int Reader::GetNumberElements(){
  return connec->getNbDof();
}
/* -------------------------------------------------------------------------- */



}
