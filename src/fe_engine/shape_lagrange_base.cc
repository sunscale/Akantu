/**
 * @file   shape_lagrange_base.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Aug 09 2017
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  common par for the shape lagrange
 *
 * @section LICENSE
 *
 * Copyright (©) 2016-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "shape_lagrange_base.hh"
#include "mesh_iterators.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

ShapeLagrangeBase::ShapeLagrangeBase(const Mesh & mesh, UInt spatial_dimension,
                                     const ElementKind & kind, const ID & id,
                                     const MemoryID & memory_id)
    : ShapeFunctions(mesh, spatial_dimension, id, memory_id), _kind(kind) {}

/* -------------------------------------------------------------------------- */
ShapeLagrangeBase::~ShapeLagrangeBase() = default;

/* -------------------------------------------------------------------------- */
#define AKANTU_COMPUTE_SHAPES(type)                                            \
  _this.template computeShapesOnIntegrationPoints<type>(                       \
      nodes, integration_points, shapes, ghost_type, filter_elements)

namespace shape_lagrange {
  namespace details {
    template <ElementKind kind> struct Helper {
      template <class S>
      static void call(const S &, const Array<Real> &, const Matrix<Real> &,
                       Array<Real> &, const ElementType &, const GhostType &,
                       const Array<UInt> &) {
        AKANTU_TO_IMPLEMENT();
      }
    };

#define AKANTU_COMPUTE_SHAPES_KIND(kind)                                       \
  template <> struct Helper<kind> {                                            \
    template <class S>                                                         \
    static void call(const S & _this, const Array<Real> & nodes,               \
                     const Matrix<Real> & integration_points,                  \
                     Array<Real> & shapes, const ElementType & type,           \
                     const GhostType & ghost_type,                             \
                     const Array<UInt> & filter_elements) {                    \
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(AKANTU_COMPUTE_SHAPES, kind);           \
    }                                                                          \
  };

    AKANTU_BOOST_ALL_KIND_LIST(AKANTU_COMPUTE_SHAPES_KIND,
                               AKANTU_FE_ENGINE_LIST_LAGRANGE_BASE)

  } // namespace details
} // namespace shape_lagrange

/* -------------------------------------------------------------------------- */
void ShapeLagrangeBase::computeShapesOnIntegrationPoints(
    const Array<Real> & nodes, const Matrix<Real> & integration_points,
    Array<Real> & shapes, const ElementType & type,
    const GhostType & ghost_type, const Array<UInt> & filter_elements) const {

  auto kind = Mesh::getKind(type);

#define AKANTU_COMPUTE_SHAPES_KIND_SWITCH(kind)                                \
  shape_lagrange::details::Helper<kind>::call(                                 \
      *this, nodes, integration_points, shapes, type, ghost_type,              \
      filter_elements);

  AKANTU_BOOST_LIST_SWITCH(
      AKANTU_COMPUTE_SHAPES_KIND_SWITCH,
      BOOST_PP_LIST_TO_SEQ(AKANTU_FE_ENGINE_LIST_LAGRANGE_BASE), kind);

#undef AKANTU_COMPUTE_SHAPES
#undef AKANTU_COMPUTE_SHAPES_KIND
#undef AKANTU_COMPUTE_SHAPES_KIND_SWITCH
}

/* -------------------------------------------------------------------------- */
void ShapeLagrangeBase::onElementsAdded(const Array<Element> & new_elements) {
  AKANTU_DEBUG_IN();
  const auto & nodes = mesh.getNodes();

  for (auto elements_range : MeshElementsByTypes(new_elements)) {
    auto type = elements_range.getType();
    auto ghost_type = elements_range.getGhostType();

    if (mesh.getSpatialDimension(type) != _spatial_dimension)
      continue;

    if (mesh.getKind(type) != _kind)
      continue;

    auto & elements = elements_range.getElements();

    auto itp_type = FEEngine::getInterpolationType(type);

    if (not this->shapes_derivatives.exists(itp_type, ghost_type)) {
      auto size_of_shapesd = this->getShapeDerivativesSize(type);
      this->shapes_derivatives.alloc(0, size_of_shapesd, itp_type, ghost_type);
    }

    if (not shapes.exists(itp_type, ghost_type)) {
      auto size_of_shapes = this->getShapeSize(type);
      this->shapes.alloc(0, size_of_shapes, itp_type, ghost_type);
    }

    const auto & natural_coords = integration_points(type, ghost_type);
    computeShapesOnIntegrationPoints(nodes, natural_coords,
                                     shapes(itp_type, ghost_type), type,
                                     ghost_type, elements);

    computeShapeDerivativesOnIntegrationPoints(
        nodes, natural_coords, shapes_derivatives(itp_type, ghost_type), type,
        ghost_type, elements);
  }
#undef INIT_SHAPE_FUNCTIONS

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ShapeLagrangeBase::onElementsRemoved(
    const Array<Element> &, const ElementTypeMapArray<UInt> & new_numbering) {
  this->shapes.onElementsRemoved(new_numbering);
  this->shapes_derivatives.onElementsRemoved(new_numbering);
}

/* -------------------------------------------------------------------------- */
void ShapeLagrangeBase::printself(std::ostream & stream, int indent) const {
  std::string space(indent, AKANTU_INDENT);

  stream << space << "Shapes Lagrange [" << std::endl;
  ShapeFunctions::printself(stream, indent + 1);
  shapes.printself(stream, indent + 1);
  shapes_derivatives.printself(stream, indent + 1);
  stream << space << "]" << std::endl;
}

} // namespace akantu
