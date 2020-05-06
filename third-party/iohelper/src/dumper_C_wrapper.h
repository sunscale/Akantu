/**
 * @file   dumper_C_wrapper.h
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date   Thu Mar 11 12:45:06 2010
 *
 * @brief  header for C api
 *
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef DUMPER_C_WRAPPER
#define DUMPER_C_WRAPPER

#define C_MODE 0
#define FORTRAN_MODE 1
#define DIM 3

#define TRIANGLE1 0
#define TRIANGLE2 1
#define TETRA1    2
#define TETRA2    3
#define POINT_SET 4
#define LINE1     5
#define LINE2     6
#define QUAD1     7
#define QUAD2     8
#define HEX1      9
#define BEAM2     10
#define HEX2      11

#define COMPRESSED 2
#define BASE64 1
#define TEXT 0

#define PARAVIEW 0
#define RESTART 1

typedef struct DumpHelper_{
  void * object_ptr;
}DumpHelper;

extern DumpHelper * getNewDumperHandle(int dumper_style);
//! dump to file
extern void Dump(DumpHelper * pH);
//! initialisation of the dumper
extern void DumperInit(DumpHelper * pH);
//! give vector with coordinates
extern void DumperSetPoints(DumpHelper * pH,double * points,int dimension,int nb,const char * name);
//! give vector to connectivity
extern void DumperSetConnectivity(DumpHelper * pH,int * connectivity,int element_type,int nb_elem,int mode);
//! give vector to per node data
extern void DumperAddNodeDataField(DumpHelper * pH,double * data,int dimension,const char * name);
//! give vector to per element data
extern void DumperAddElemDataField(DumpHelper * pH,double * data,int dimension,const char * name);
//! set mode for file creation : TEXT, BASE64, COMPRESSED
extern void DumperSetMode(DumpHelper * pH,int mode);
//! set flag embedded dimension. This allow paraview to complement 1D and 2D to obtain 3D fields
extern void DumperSetEmbeddedValue(DumpHelper * ph,const char * name,int flag);
//! set prefix directory
extern void DumperSetPrefix(DumpHelper * ph,const char * dir);
//! set rank and world size params for parallel treatment
extern void DumperSetParallelContext(DumpHelper * pH,const int me,const int wld_size);
//! free memory
extern void DumperFreeHandle(DumpHelper * pH);
//! set number of filtered elements
extern void DumperSetNumberFilteredElements(DumpHelper * pH,int nb_filtered);
//! set internal dumper count
extern void DumperSetDumpStep(DumpHelper * pH,int step);
#endif //DUMPER_C_WRAPPER
