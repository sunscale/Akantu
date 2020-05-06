/**
 * @file   shape_lagrange_base_inline_impl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Aug 09 2017
 * @date last modification: Wed Oct 11 2017
 *
 * @brief  A Documented file.
 *
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
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SHAPE_LAGRANGE_BASE_INLINE_IMPL_HH__
#define __AKANTU_SHAPE_LAGRANGE_BASE_INLINE_IMPL_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeLagrangeBase::computeShapesOnIntegrationPoints(
    const Array<Real> &, const Matrix<Real> & integration_points,
    Array<Real> & shapes, const GhostType & ghost_type,
    const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  UInt nb_points = integration_points.cols();
  UInt nb_element = mesh.getConnectivity(type, ghost_type).size();

  shapes.resize(nb_element * nb_points);

#if !defined(AKANTU_NDEBUG)
  UInt size_of_shapes = ElementClass<type>::getShapeSize();
  AKANTU_DEBUG_ASSERT(shapes.getNbComponent() == size_of_shapes,
                      "The shapes array does not have the correct "
                          << "number of component");
#endif

  auto shapes_it = shapes.begin_reinterpret(
      ElementClass<type>::getNbNodesPerInterpolationElement(), nb_points,
      nb_element);
  auto shapes_begin = shapes_it;
  if (filter_elements != empty_filter) {
    nb_element = filter_elements.size();
  }

  for (UInt elem = 0; elem < nb_element; ++elem) {
    if (filter_elements != empty_filter)
      shapes_it = shapes_begin + filter_elements(elem);

    Matrix<Real> & N = *shapes_it;
    ElementClass<type>::computeShapes(integration_points, N);

    if (filter_elements == empty_filter)
      ++shapes_it;
  }

  AKANTU_DEBUG_OUT();
}

} // namespace akantu

#endif /* __AKANTU_SHAPE_LAGRANGE_BASE_INLINE_IMPL_HH__ */
