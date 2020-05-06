/**
 * @file   reader_C_wrapper.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Thu Mar 11 2010
 * @date last modification: Thu Nov 01 2012
 *
 * @brief  wrapper for C interface
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

#include <string>
#include "iohelper_common.hh"

#include "reader.hh"
#include "reader_restart.hh"

extern "C" {
#include "reader_C_wrapper.h"
}

extern "C" ReadHelper * getNewReaderHandle(int reader_style){
  ReadHelper * box;
  switch (reader_style){
  case RESTART:
    {
      ReaderRestart * obj = new ReaderRestart();
      box = (ReadHelper*)malloc(sizeof(ReadHelper));
      box->object_ptr = obj;
    }
    break;
  default:
    FATAL("unknown reader style");
  }
  return box;
}
extern "C" void Read(ReadHelper * pH){
  Reader * ptr = (Reader *)pH->object_ptr;
  ptr->Read();
}
extern "C" void ReaderInit(ReadHelper * pH){
  Reader * ptr = (Reader *)pH->object_ptr;
  ptr->Init();
}
extern "C" void ReaderSetPoints(ReadHelper * pH,const char * name){
  Reader * ptr = (Reader *)pH->object_ptr;
  ptr->SetPoints(name);
}
extern "C" void ReaderSetConnectivity(ReadHelper * pH,int element_type){
  Reader * ptr = (Reader *)pH->object_ptr;
  ptr->SetConnectivity(element_type);
}
extern "C" void ReaderAddNodeDataField(ReadHelper * pH,const char * name){
  Reader * ptr = (Reader *)pH->object_ptr;
  ptr->AddNodeDataField(name);
}
extern "C" void ReaderAddElemDataField(ReadHelper * pH,const char * name){
  Reader * ptr = (Reader *)pH->object_ptr;
  ptr->AddElemDataField(name);
}
extern "C" void ReaderSetPrefix(ReadHelper * pH,const char * dir){
  Reader * ptr = (Reader *)pH->object_ptr;
  ptr->SetPrefix(dir);
}
extern "C" void ReaderSetParallelContext(ReadHelper * pH,const int me,const int wld_size){
  Reader * ptr = (Reader *)pH->object_ptr;
  ptr->SetParallelContext(me,wld_size);
}
extern "C" void ReaderFreeHandle(ReadHelper * pH){
  Reader * ptr = (Reader *)pH->object_ptr;
  delete(ptr);
  free(pH);
}
extern "C" double * ReaderGetPoints(ReadHelper * pH){
  Reader * ptr = (Reader *)pH->object_ptr;
  return ptr->GetPoints();
}
extern "C" int * ReaderGetConnectivity(ReadHelper * pH){
  Reader * ptr = (Reader *)pH->object_ptr;
  return ptr->GetConnectivity();
}
extern "C" double * ReaderGetNodeDataField(ReadHelper * pH,const char * name){
  Reader * ptr = (Reader *)pH->object_ptr;
  return ptr->GetNodeDataField(name);
}
extern "C" double * ReaderGetElemDataField(ReadHelper * pH,const char * name){
  Reader * ptr = (Reader *)pH->object_ptr;
  return ptr->GetElemDataField(name);
}
extern "C" int ReaderGetNumberNodes(ReadHelper * pH){
  Reader * ptr = (Reader *)pH->object_ptr;
  return ptr->GetNumberNodes();
}
extern "C" int ReaderGetNumberElements(ReadHelper * pH){
  Reader * ptr = (Reader *)pH->object_ptr;
  return ptr->GetNumberElements();
}
extern "C" void ReaderSetMode(ReadHelper * pH,int mode){
  Reader * ptr = (Reader *)pH->object_ptr;
  ptr->SetMode(mode);
}
