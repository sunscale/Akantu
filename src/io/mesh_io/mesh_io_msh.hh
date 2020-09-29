/**
 * @file   mesh_io_msh.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  Read/Write for MSH files
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

#ifndef AKANTU_MESH_IO_MSH_HH_
#define AKANTU_MESH_IO_MSH_HH_

/* -------------------------------------------------------------------------- */
#include "mesh_io.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

class MeshIOMSH : public MeshIO {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MeshIOMSH();

  ~MeshIOMSH() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// read a mesh from the file
  void read(const std::string & filename, Mesh & mesh) override;

  /// write a mesh to a file
  void write(const std::string & filename, const Mesh & mesh) override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// MSH element types
  enum MSHElementType {
    _msh_not_defined = 0,
    _msh_segment_2 = 1,       // 2-node line.
    _msh_triangle_3 = 2,      // 3-node triangle.
    _msh_quadrangle_4 = 3,    // 4-node quadrangle.
    _msh_tetrahedron_4 = 4,   // 4-node tetrahedron.
    _msh_hexahedron_8 = 5,    // 8-node hexahedron.
    _msh_prism_1 = 6,         // 6-node prism.
    _msh_pyramid_1 = 7,       // 5-node pyramid.
    _msh_segment_3 = 8,       // 3-node second order line
    _msh_triangle_6 = 9,      // 6-node second order triangle
    _msh_quadrangle_9 = 10,   // 9-node second order quadrangle
    _msh_tetrahedron_10 = 11, // 10-node second order tetrahedron
    _msh_hexahedron_27 = 12,  // 27-node second order hexahedron
    _msh_prism_18 = 13,       // 18-node second order prism
    _msh_pyramid_14 = 14,     // 14-node second order pyramid
    _msh_point = 15,          // 1-node point.
    _msh_quadrangle_8 = 16,   // 8-node second order quadrangle
    _msh_hexahedron_20 = 17,  // 20-node second order hexahedron
    _msh_prism_15 = 18        // 15-node second order prism
  };

#define MAX_NUMBER_OF_NODE_PER_ELEMENT 10 // tetrahedron of second order

  /// order in witch element as to be read
  std::map<ElementType, std::vector<UInt>> _read_order;

  /// number of nodes per msh element
  std::map<MSHElementType, UInt> _msh_nodes_per_elem;

  /// correspondence between msh element types and akantu element types
  std::map<MSHElementType, ElementType> _msh_to_akantu_element_types;

  /// correspondence between akantu element types and msh element types
  std::map<ElementType, MSHElementType> _akantu_to_msh_element_types;

protected:
  template <typename File, typename Readers>
  void populateReaders2(File & file, Readers & readers);

  template <typename File, typename Readers>
  void populateReaders4(File & file, Readers & readers);
};

} // namespace akantu

#endif /* AKANTU_MESH_IO_MSH_HH_ */
