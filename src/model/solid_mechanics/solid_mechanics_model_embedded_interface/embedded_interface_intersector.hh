/**
 * @file   embedded_interface_intersector.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Fri May 01 2015
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  Class that loads the interface from mesh and computes intersections
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef AKANTU_EMBEDDED_INTERFACE_INTERSECTOR_HH_
#define AKANTU_EMBEDDED_INTERFACE_INTERSECTOR_HH_

#include "aka_common.hh"
#include "mesh_geom_abstract.hh"
#include "mesh_geom_common.hh"
#include "mesh_segment_intersector.hh"

/* -------------------------------------------------------------------------- */

namespace akantu {

namespace {
  using K = cgal::Cartesian;
}

/**
 * @brief Computes the intersections of the reinforcements defined in the
 * primitive mesh
 *
 * The purpose of this class is to look for reinforcements in the primitive
 * mesh, which
 * should be defined by physical groups with the same names as the reinforcement
 * materials
 * in the model.
 *
 * It then constructs the CGAL primitives from the elements of those
 * reinforcements
 * and computes the intersections with the background mesh, to create an
 * `interface_mesh`,
 * which is in turn used by the EmbeddedInterfaceModel.
 *
 * @see MeshSegmentIntersector, MeshGeomAbstract
 * @see EmbeddedInterfaceModel
 */
class EmbeddedInterfaceIntersector : public MeshGeomAbstract {

public:
  /// Construct from mesh and a reinforcement mesh
  explicit EmbeddedInterfaceIntersector(Mesh & mesh,
                                        const Mesh & primitive_mesh);

  /// Destructor
  ~EmbeddedInterfaceIntersector() override = default;

public:
  /// Generate the interface mesh
  void constructData(GhostType ghost_type = _not_ghost) override;

  /// Create a segment with an element connectivity
  K::Segment_3 createSegment(const Vector<UInt> & connectivity);

  /// Getter for interface mesh
  AKANTU_GET_MACRO_NOT_CONST(InterfaceMesh, interface_mesh, Mesh &);

protected:
  /// Resulting mesh of intersection
  Mesh interface_mesh;

  /// Mesh used for primitive construction
  const Mesh & primitive_mesh;
};

} // namespace akantu

#endif // AKANTU_EMBEDDED_INTERFACE_INTERSECTOR_HH_
