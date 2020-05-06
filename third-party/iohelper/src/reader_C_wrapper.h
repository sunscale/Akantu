/**
 * @file   reader_C_wrapper.h
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date   Thu Mar 11 12:45:06 2010
 *
 * @brief  Reader C api header
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

#ifndef READER_C_WRAPPER
#define READER_C_WRAPPER

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

#define RESTART 1

typedef struct ReadHelper_{
  void * object_ptr;
}ReadHelper;

extern ReadHelper * getNewReaderHandle(int reader_style);
//! read to file
extern void Read(ReadHelper * pH);
//! initialisation of the reader
extern void ReaderInit(ReadHelper * pH);
//! give vector with coordinates
extern void ReaderSetPoints(ReadHelper * pH,const char * name);
//! give vector to connectivity
extern void ReaderSetConnectivity(ReadHelper * pH,int element_type);
//! give vector to per node data
extern void ReaderAddNodeDataField(ReadHelper * pH,const char * name);
//! give vector to per element data
extern void ReaderAddElemDataField(ReadHelper * pH,const char * name);
//! set prefix directory
extern void ReaderSetPrefix(ReadHelper * ph,const char * dir);
//! set rank and world size params for parallel treatment
extern void ReaderSetParallelContext(ReadHelper * pH,const int me,const int wld_size);
//! free memory
extern void ReaderFreeHandle(ReadHelper * pH);

//! give vector with coordinates
extern double * ReaderGetPoints(ReadHelper* pH);
//! give vector to connectivity
extern int * ReaderGetConnectivity(ReadHelper * pH);
//! give vector to per node data
extern double * ReaderGetNodeDataField(ReadHelper * pH, const char *name);
//! give vector to per element data
extern double * ReaderGetElemDataField(ReadHelper * pH,const char *name);
//! give number of nodes
extern int ReaderGetNumberNodes(ReadHelper * pH);
//! give number of elements
extern int ReaderGetNumberElements(ReadHelper * pH);
//! set mode for file creation : TEXT, BASE64, COMPRESSED
extern void ReaderSetMode(ReadHelper * pH,int mode);
#endif //READERPARAVIEW_C_WRAPPER
