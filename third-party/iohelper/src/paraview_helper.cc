/**
 * @file   paraview_helper.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Oct 12 2012
 * @date last modification: Wed Jun 05 2013
 *
 * @brief  implementation of paraview helper
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

#include "paraview_helper.hh"
#include <fstream>


/* -------------------------------------------------------------------------- */
#if defined(__INTEL_COMPILER)
/// remark #981: operands are evaluated in unspecified order
#pragma warning ( disable : 981 )
/// remark #383: value copied to temporary, reference to temporary used
#pragma warning ( disable : 383 )
#endif //defined(__INTEL_COMPILER)

namespace iohelper {

/* -------------------------------------------------------------------------- */
ParaviewHelper::ParaviewHelper(File & f, UInt mode):
  b64(f), file(f), position_flag(false) {
  bflag = BASE64;
  compteur = 0;
  setMode(mode);

  this->paraview_code_type[TRIANGLE1] = VTK_TRIANGLE;
  this->paraview_code_type[TRIANGLE2] = VTK_QUADRATIC_TRIANGLE;
  this->paraview_code_type[TETRA1   ] = VTK_TETRA;
  this->paraview_code_type[TETRA2   ] = VTK_QUADRATIC_TETRA;
  this->paraview_code_type[POINT_SET] = VTK_POLY_VERTEX;
  this->paraview_code_type[LINE1    ] = VTK_LINE;
  this->paraview_code_type[LINE2    ] = VTK_QUADRATIC_EDGE;
  this->paraview_code_type[QUAD1    ] = VTK_QUAD;
  this->paraview_code_type[QUAD2    ] = VTK_QUADRATIC_QUAD;
  this->paraview_code_type[HEX1     ] = VTK_HEXAHEDRON;
  this->paraview_code_type[HEX2     ] = VTK_QUADRATIC_HEXAHEDRON;
  this->paraview_code_type[BEAM2    ] = VTK_LINE;
  this->paraview_code_type[BEAM3    ] = VTK_LINE;
  this->paraview_code_type[PRISM1   ] = VTK_WEDGE;
  this->paraview_code_type[PRISM2   ] = VTK_QUADRATIC_WEDGE;
  this->paraview_code_type[COH1D2   ] = VTK_LINE;
  this->paraview_code_type[COH2D4   ] = VTK_POLYGON;
  this->paraview_code_type[COH2D6   ] = VTK_POLYGON;
  this->paraview_code_type[COH3D6   ] = VTK_WEDGE;
  this->paraview_code_type[COH3D12  ] = VTK_QUADRATIC_LINEAR_WEDGE;
  this->paraview_code_type[COH3D8   ] = VTK_HEXAHEDRON;


  std::map<ElemType, VTKCellType>::iterator it;
  for(it = paraview_code_type.begin();
      it != paraview_code_type.end(); ++it) {
    UInt nb_nodes = nb_node_per_elem[it->first];

    UInt * tmp = new UInt[nb_nodes];
    for (UInt i = 0; i < nb_nodes; ++i) {
      tmp[i] = i;
    }

    switch(it->first) {
    case COH3D12:
      tmp[ 0] = 0;
      tmp[ 1] = 1;
      tmp[ 2] = 2;
      tmp[ 3] = 6;
      tmp[ 4] = 7;
      tmp[ 5] = 8;
      tmp[ 6] = 3;
      tmp[ 7] = 4;
      tmp[ 8] = 5;
      tmp[ 9] = 9;
      tmp[10] = 10;
      tmp[11] = 11;
      break;
    case COH2D6:
      tmp[0] = 0;
      tmp[1] = 2;
      tmp[2] = 1;
      tmp[3] = 4;
      tmp[4] = 5;
      tmp[5] = 3;
      break;
    case COH2D4:
      tmp[0] = 0;
      tmp[1] = 1;
      tmp[2] = 3;
      tmp[3] = 2;
      break;
    case HEX2:
      tmp[12] = 16;
      tmp[13] = 17;
      tmp[14] = 18;
      tmp[15] = 19;
      tmp[16] = 12;
      tmp[17] = 13;
      tmp[18] = 14;
      tmp[19] = 15;
      break;
    case PRISM2:
      tmp[ 0] =  0;
      tmp[ 1] =  1;
      tmp[ 2] =  2;
      
      tmp[ 3] =  3;
      tmp[ 4] =  4;
      tmp[ 5] =  5;
      
      tmp[ 6] =  6;
      tmp[ 7] =  7;
      tmp[ 8] =  8;

      tmp[ 9] = 12;
      tmp[10] = 13;
      tmp[11] = 14;
      
      tmp[12] = 9;
      tmp[13] = 10;
      tmp[14] = 11;
      
      break;
    default:
      //nothing to change
      break;
    }
    this->write_reorder[it->first] = tmp;
  }

}

/* -------------------------------------------------------------------------- */
ParaviewHelper::~ParaviewHelper(){
  std::map<ElemType, VTKCellType>::iterator it;
  for(it = this->paraview_code_type.begin();
      it != this->paraview_code_type.end(); ++it) {
    delete [] this->write_reorder[it->first];
  }
}

/* -------------------------------------------------------------------------- */
void ParaviewHelper::writeTimePVD(const std::string & filename,
				  const std::vector< std::pair<Real, std::string> > & pvtus) {
  std::ofstream pvd_file;
  pvd_file.open(filename.c_str());

  if(!pvd_file.good()) {
    IOHELPER_THROW("DumperParaview was not able to open the file \"" << filename, _et_file_error);
  }

  pvd_file << "<?xml version=\"1.0\"?>" << std::endl
	   << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl
	   << "  <Collection>" << std::endl;

  auto it = pvtus.begin();
  auto end = pvtus.end();
  for (;it != end; ++it) {
    pvd_file << "    <DataSet timestep=\"" << it->first << "\" group=\"\" part=\"0\" file=\""
	     << it->second << "\"/>" << std::endl;
  }

  pvd_file << "  </Collection>" << std::endl
	   << "</VTKFile>" << std::endl;
  pvd_file.close();
}

/* -------------------------------------------------------------------------- */
void ParaviewHelper::writeHeader(int nb_nodes,int nb_elems){
  file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" " ;
  file << "byte_order=\"LittleEndian\">" << std::endl;
  file << " <UnstructuredGrid>" << std::endl
       << "  <Piece NumberOfPoints= \""
       << nb_nodes << "\" NumberOfCells=\""
       << nb_elems << "\">" << std::endl;
}

/* -------------------------------------------------------------------------- */
void ParaviewHelper::PDataArray(const std::string & name, int nb_components, const std::string & type){
  file << "   <PDataArray type=\"" << type << "\" NumberOfComponents=\""
       << nb_components << "\" Name=\"" << name << "\" format=\"";
  if (bflag == BASE64) {
    file << "binary";
  } else {
    file << "ascii";
  }

  file << "\"></PDataArray>" << std::endl;
}

/* -------------------------------------------------------------------------- */
void ParaviewHelper::write_conclusion(){
  file << "  </Piece>" << std::endl;
  file << " </UnstructuredGrid>" << std::endl;
  file << "</VTKFile>" << std::endl;
}

/* -------------------------------------------------------------------------- */
void ParaviewHelper::startData(const std::string & name,
			       int nb_components,
			       const std::string & type){

  file << "    <DataArray type=\"" << type << "\" ";
  if (nb_components) {
    file << "NumberOfComponents=\"" << nb_components << "\" ";
  }

  file << "Name=\"" << name << "\" format=\"";
  if (bflag == BASE64) {
    file << "binary";
  } else {
    file << "ascii";
  }
  file << "\">" << std::endl;
  if (bflag == BASE64) {
    b64.CreateHeader();
  }
}

/* -------------------------------------------------------------------------- */
void ParaviewHelper::endData(){
  if (bflag == BASE64) {
    b64.WriteHeader();
  }
  file << std::endl << "    </DataArray>" << std::endl;
}

/* -------------------------------------------------------------------------- */
void ParaviewHelper::startDofList(int dimension){
  file << "   <Points>" << std::endl;
  startData("positions", dimension, "Float64");
}

/* -------------------------------------------------------------------------- */
void ParaviewHelper::endDofList(){
  endData();
  file << "   </Points>" << std::endl;
}

/* -------------------------------------------------------------------------- */
void ParaviewHelper::startCells(){
  file << "   <Cells>" << std::endl;
}

/* -------------------------------------------------------------------------- */
void ParaviewHelper::endCells(){
  file << "   </Cells>" << std::endl;
}

/* -------------------------------------------------------------------------- */
void ParaviewHelper::startCellsConnectivityList(){
  startData("connectivity",0,"Int32");
}

/* -------------------------------------------------------------------------- */
void ParaviewHelper::endCellsConnectivityList(){
  endData();
}

/* -------------------------------------------------------------------------- */
void ParaviewHelper::startCellsoffsetsList(){
  startData("offsets",0,"Int32");
}

/* -------------------------------------------------------------------------- */
void ParaviewHelper::endCellsoffsetsList(){
  endData();
}

/* -------------------------------------------------------------------------- */
void ParaviewHelper::startCellstypesList(){
  startData("types",0,"UInt32");
}

/* -------------------------------------------------------------------------- */
void ParaviewHelper::endCellstypesList(){
  endData();
}

/* -------------------------------------------------------------------------- */
void ParaviewHelper::startPointDataList(){
  file << "   <PointData>" << std::endl;
}

/* -------------------------------------------------------------------------- */
void ParaviewHelper::endPointDataList(){
  file << "   </PointData>" << std::endl;
}

/* -------------------------------------------------------------------------- */
void ParaviewHelper::startCellDataList(){
  file << "   <CellData>" << std::endl;
}

/* -------------------------------------------------------------------------- */
void ParaviewHelper::endCellDataList(){
  file << "   </CellData>" << std::endl;
}

/* -------------------------------------------------------------------------- */

}
