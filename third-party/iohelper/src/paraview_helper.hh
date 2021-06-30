/**
 * @file   paraview_helper.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Thu Mar 11 2010
 * @date last modification: Wed Jun 05 2013
 *
 * @brief  paraview helper header
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

#ifndef IOHELPER_PARAVIEW_HELPER_H_
#define IOHELPER_PARAVIEW_HELPER_H_
/* -------------------------------------------------------------------------- */
#include "base64.hh"
#include <iomanip>
#include <map>
#include "visitor.hh"
#include "field_interface.hh"
/* -------------------------------------------------------------------------- */

namespace iohelper {

// Taken from vtkCellType.h
enum VTKCellType {
  // Linear cells
  VTK_EMPTY_CELL = 0,
  VTK_VERTEX = 1,
  VTK_POLY_VERTEX = 2,
  VTK_LINE = 3,
  VTK_POLY_LINE = 4,
  VTK_TRIANGLE = 5,
  VTK_TRIANGLE_STRIP = 6,
  VTK_POLYGON = 7,
  VTK_PIXEL = 8,
  VTK_QUAD = 9,
  VTK_TETRA = 10,
  VTK_VOXEL = 11,
  VTK_HEXAHEDRON = 12,
  VTK_WEDGE = 13,
  VTK_PYRAMID = 14,
  VTK_PENTAGONAL_PRISM = 15,
  VTK_HEXAGONAL_PRISM = 16,

  // Quadratic, isoparametric cells
  VTK_QUADRATIC_EDGE = 21,
  VTK_QUADRATIC_TRIANGLE = 22,
  VTK_QUADRATIC_QUAD = 23,
  VTK_QUADRATIC_POLYGON = 36,
  VTK_QUADRATIC_TETRA = 24,
  VTK_QUADRATIC_HEXAHEDRON = 25,
  VTK_QUADRATIC_WEDGE = 26,
  VTK_QUADRATIC_PYRAMID = 27,
  VTK_BIQUADRATIC_QUAD = 28,
  VTK_TRIQUADRATIC_HEXAHEDRON = 29,
  VTK_QUADRATIC_LINEAR_QUAD = 30,
  VTK_QUADRATIC_LINEAR_WEDGE = 31,
  VTK_BIQUADRATIC_QUADRATIC_WEDGE = 32,
  VTK_BIQUADRATIC_QUADRATIC_HEXAHEDRON = 33,
  VTK_BIQUADRATIC_TRIANGLE = 34,

  // Polyhedron cell (consisting of polygonal faces)
  VTK_POLYHEDRON = 42,

  // Higher order cells in parametric form
  VTK_PARAMETRIC_CURVE = 51,
  VTK_PARAMETRIC_SURFACE = 52,
  VTK_PARAMETRIC_TRI_SURFACE = 53,
  VTK_PARAMETRIC_QUAD_SURFACE = 54,
  VTK_PARAMETRIC_TETRA_REGION = 55,
  VTK_PARAMETRIC_HEX_REGION = 56,

  // Higher order cells
  VTK_HIGHER_ORDER_EDGE = 60,
  VTK_HIGHER_ORDER_TRIANGLE = 61,
  VTK_HIGHER_ORDER_QUAD = 62,
  VTK_HIGHER_ORDER_POLYGON = 63,
  VTK_HIGHER_ORDER_TETRAHEDRON = 64,
  VTK_HIGHER_ORDER_WEDGE = 65,
  VTK_HIGHER_ORDER_PYRAMID = 66,
  VTK_HIGHER_ORDER_HEXAHEDRON = 67,

  VTK_NUMBER_OF_CELL_TYPES
};


inline std::ostream & operator <<(std::ostream & stream, const VTKCellType & type) {
  stream << UInt(type);
  return stream;
}


class ParaviewHelper : public Visitor {

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */

  enum Stage{
    _s_writePosition,
    _s_writeFieldProperty,
    _s_writeField,
    _s_writeConnectivity,
    _s_writeElemType,
    _s_buildOffsets
  };

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

 public:

  ParaviewHelper(File & f, UInt mode);
  ~ParaviewHelper() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */


  //! write the PVTU file
  template <typename T>
  void writePVTU(T & per_node_data, T & per_elem_data,
		 const std::vector<std::string> & vtus);

  //! write the PVD file for time description
  static void writeTimePVD(const std::string & filename,
			   const std::vector< std::pair<Real, std::string> > & pvtus);

  //! write the header of a vtu file
  void writeHeader(int nb_nodes,int nb_elems);
  //! write a field
  template <typename T> void writeField(T & data);
  //! write a connectivity field
  template <typename T>  void writeConnectivity(T & data);
  //! write an element type field
  template <typename T> void writeElemType(T & data);
  //! write the field properties
  template <typename T> void writeFieldProperty(T & data);
  //! write the connectivities offset
  template <typename T> void writeOffsets(T & data);


  template <typename T>
  void pushDataFields(T & per_node_data, T & per_elem_data);

  //! push the position field to the paraview file
  void pushPosition(FieldInterface & f);
  //! push a field to the paraview file
  void pushField(FieldInterface & f);
  //! build the offset from connectivities
  void buildOffsets(FieldInterface & f);
  //! push a connectivity field
  void pushConnectivity(FieldInterface & f);
  //! push a element type field
  void pushElemType(FieldInterface & f);

  //! get the formated vtu name
  // static std::string getVTUName(const std::string & basename, UInt proc);

  //! push a small array of values
  template <template<typename T> class Cont, typename T>
  void pushData(const Cont<T> & n);

  //! push a small array of values of homogeneous values with padding to size dim
  template <template<typename T> class Cont, typename T>
  inline void pushData(const Cont<T> & n, UInt dim);

  //! pushing datum
  template <typename T> void pushDatum(const T & n, UInt size = 3);

  //! visitor system
  template <typename T> void visitField(T & visited);
private:

  void setMode(int mode);

  static std::string dataTypeToStr(DataType data_type);

  /* ------------------------------------------------------------------------ */
  /* Methods for writing control sequences in the paraview files              */
  /* ------------------------------------------------------------------------ */

public:
  void startDofList(int dimension);
  void endDofList();
  void startCells();
  void endCells();
  void startCellsConnectivityList();
  void endCellsConnectivityList();
  void startCellsoffsetsList();
  void endCellsoffsetsList();
  void startCellstypesList();
  void endCellstypesList();
  void startPointDataList();
  void endPointDataList();
  void startCellDataList();
  void endCellDataList();
  void startData(const std::string & name, int nb_components, const std::string & type);
  void PDataArray(const std::string & name, int nb_components, const std::string & type);
  void endData();
  void write_conclusion();

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

  Base64Writer b64;
  int bflag;

  File & file;
  long header_offset;

  unsigned int compteur;

  Stage current_stage;

  bool position_flag;

  //! mapping between iohelper elements and paraview elements
  std::map<ElemType, VTKCellType> paraview_code_type;

  //! mapping of the connectivities between iohelper and paraview
  std::map<ElemType, UInt *> write_reorder;
};

/* -------------------------------------------------------------------------- */
#include "paraview_helper.tcc"
/* -------------------------------------------------------------------------- */




}



#endif /* IOHELPER_PARAVIEW_HELPER_H_ */
