/**
 * @file   mesh_geom_abstract.hh
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Thu Feb 26 2015
 * @date last modification: Fri Mar 6 2015
 *
 * @brief  Class for constructing the CGAL primitives of a mesh
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MESH_GEOM_ABSTRACT_HH__
#define __AKANTU_MESH_GEOM_ABSTRACT_HH__

#include "aka_common.hh"
#include "mesh.hh"

#include <CGAL/Cartesian.h>

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

typedef CGAL::Cartesian<Real> K;

class MeshGeomAbstract {

public:
  /// Construct from mesh
  explicit MeshGeomAbstract(const Mesh & mesh);

  /// Destructor
  virtual ~MeshGeomAbstract();

public:
  /// Construct geometric data for computational geometry algorithms
  virtual void constructData() = 0;

  /// Compute number of intersections with geometric interface
  virtual UInt numberOfIntersectionsWithInterface(const K::Segment_3 & interface) const = 0;

  /// Compute the mesh created by a linear interface
  virtual Mesh * meshOfLinearInterface(const std::pair<K::Segment_3, std::string> & interface) = 0;

protected:
  /// Mesh used to construct the primitives
  const Mesh & mesh;
};

__END_AKANTU__

#endif // __AKANTU_MESH_GEOM_ABSTRACT_HH__
