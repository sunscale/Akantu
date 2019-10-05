/**
 * @file   aka_element_classes_info_inline_impl.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Jun 18 2015
 * @date last modification: Wed Jan 10 2018
 *
 * @brief  Implementation of the streaming fonction for the element classes
 * enums
 *
 * @section LICENSE
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
#include <unordered_map>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_AKA_ELEMENT_CLASSES_INFO_INLINE_IMPL_CC__
#define __AKANTU_AKA_ELEMENT_CLASSES_INFO_INLINE_IMPL_CC__

namespace akantu {

AKANTU_ENUM_OUTPUT_STREAM(
    ElementType, AKANTU_ALL_ELEMENT_TYPE(_not_defined)(_max_element_type))
AKANTU_ENUM_INPUT_STREAM(ElementType, AKANTU_ALL_ELEMENT_TYPE)

AKANTU_ENUM_OUTPUT_STREAM(InterpolationType, AKANTU_INTERPOLATION_TYPES)
AKANTU_ENUM_INPUT_STREAM(InterpolationType, AKANTU_INTERPOLATION_TYPES)

AKANTU_ENUM_OUTPUT_STREAM(ElementKind, AKANTU_ELEMENT_KIND)
AKANTU_ENUM_INPUT_STREAM(ElementKind, AKANTU_ELEMENT_KIND)

} // namespace akantu

#endif /* __AKANTU_AKA_ELEMENT_CLASSES_INFO_INLINE_IMPL_CC__ */
