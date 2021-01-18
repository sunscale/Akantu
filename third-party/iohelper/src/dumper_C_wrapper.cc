/**
 * @file   dumper_C_wrapper.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Thu Mar 11 2010
 * @date last modification: Thu Nov 01 2012
 *
 * @brief  dumper C api implementation
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

#include "dumper.hh"
#include "dumper_paraview.hh"
#include "dumper_restart.hh"
/* -------------------------------------------------------------------------- */

extern "C" {
#include "dumper_C_wrapper.h"
}

/* -------------------------------------------------------------------------- */

namespace iohelper {

extern "C" DumpHelper * getNewDumperHandle(int dumper_style){
  DumpHelper * box;
  switch (dumper_style){
  case PARAVIEW:
    {
      DumperParaview * obj = new DumperParaview();
      box = (DumpHelper*)malloc(sizeof(DumpHelper));
      box->object_ptr = obj;
    }
    break;
  case RESTART:
    {
      DumperRestart * obj = new DumperRestart();
      box = (DumpHelper*)malloc(sizeof(DumpHelper));
      box->object_ptr = obj;
    }
    break;
  default:
    FATAL("unknown dumper style");
  }
  return box;
}

/* -------------------------------------------------------------------------- */


extern "C" void Dump(DumpHelper * pH){
  Dumper * ptr = (Dumper *)pH->object_ptr;
  ptr->dump();
}
extern "C" void DumperInit(DumpHelper * pH){
  Dumper * ptr = (Dumper *)pH->object_ptr;
  ptr->init();
}
extern "C" void DumperSetPoints(DumpHelper * pH,double * points,int dimension,int nb,const char * name){
  Dumper * ptr = (Dumper *)pH->object_ptr;
  ptr->setPoints(points,dimension,nb,name);
}
extern "C" void DumperSetConnectivity(DumpHelper * pH,int * connectivity,int element_type,int nb_elem,int mode){
  Dumper * ptr = (Dumper *)pH->object_ptr;
  ptr->setConnectivity(connectivity,static_cast<iohelper::ElemType>(element_type),nb_elem,mode);
}
extern "C" void DumperAddNodeDataField(DumpHelper * pH,double * data,int dimension,const char * name){
  Dumper * ptr = (Dumper *)pH->object_ptr;
  ptr->addNodeDataField(data,dimension,name);
}
extern "C" void DumperAddElemDataField(DumpHelper * pH,double * data,int dimension,const char * name){
  Dumper * ptr = (Dumper *)pH->object_ptr;
  ptr->addElemDataField(data,dimension,name);
}
extern "C" void DumperSetMode(DumpHelper * pH,int mode){
  Dumper * ptr = (Dumper *)pH->object_ptr;
  ptr->setMode(mode);
}
extern "C" void DumperSetEmbeddedValue(DumpHelper * pH,const char * name,int value){
  Dumper * ptr = (Dumper *)pH->object_ptr;
  ptr->setEmbeddedValue(name,value);
}
extern "C" void DumperSetPrefix(DumpHelper * pH,const char * dir){
  Dumper * ptr = (Dumper *)pH->object_ptr;
  ptr->setPrefix(dir);
}
extern "C" void DumperSetParallelContext(DumpHelper * pH,const int me,const int wld_size){
  Dumper * ptr = (Dumper *)pH->object_ptr;
  ptr->setParallelContext(me,wld_size);
}

extern "C" void DumperFreeHandle(DumpHelper * pH){
  Dumper * ptr = (Dumper *)pH->object_ptr;
  delete(ptr);
  free(pH);
}

//! set number of filtered elements
extern "C" void DumperSetNumberFilteredElements(DumpHelper * pH,int nb_filtered){
  Dumper * ptr = (Dumper *)pH->object_ptr;
  ptr->setNumberFilteredElements(nb_filtered);
}

//! set internal dumper count
extern "C" void DumperSetDumpStep(DumpHelper * pH,int step){
  Dumper * ptr = (Dumper *)pH->object_ptr;
  ptr->setDumpStep(step);
}



}
