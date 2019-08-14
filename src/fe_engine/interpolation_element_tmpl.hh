/**
 * @file   interpolation_element_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Jun 06 2013
 * @date last modification: Wed Nov 29 2017
 *
 * @brief  interpolation property description
 *
 * @section LICENSE
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "element_class.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_INTERPOLATION_ELEMENT_TMPL_HH__
#define __AKANTU_INTERPOLATION_ELEMENT_TMPL_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
/* Regular Elements                                                           */
/* -------------------------------------------------------------------------- */
AKANTU_DEFINE_INTERPOLATION_TYPE_PROPERTY(_itp_not_defined, _itk_not_defined, 0,
                                          0);
AKANTU_DEFINE_INTERPOLATION_TYPE_PROPERTY(_itp_lagrange_point_1,
                                          _itk_lagrangian, 1, 0);
AKANTU_DEFINE_INTERPOLATION_TYPE_PROPERTY(_itp_lagrange_segment_2,
                                          _itk_lagrangian, 2, 1);
AKANTU_DEFINE_INTERPOLATION_TYPE_PROPERTY(_itp_lagrange_segment_3,
                                          _itk_lagrangian, 3, 1);
AKANTU_DEFINE_INTERPOLATION_TYPE_PROPERTY(_itp_lagrange_triangle_3,
                                          _itk_lagrangian, 3, 2);
AKANTU_DEFINE_INTERPOLATION_TYPE_PROPERTY(_itp_lagrange_triangle_6,
                                          _itk_lagrangian, 6, 2);
AKANTU_DEFINE_INTERPOLATION_TYPE_PROPERTY(_itp_lagrange_tetrahedron_4,
                                          _itk_lagrangian, 4, 3);
AKANTU_DEFINE_INTERPOLATION_TYPE_PROPERTY(_itp_lagrange_tetrahedron_10,
                                          _itk_lagrangian, 10, 3);
AKANTU_DEFINE_INTERPOLATION_TYPE_PROPERTY(_itp_lagrange_quadrangle_4,
                                          _itk_lagrangian, 4, 2);
AKANTU_DEFINE_INTERPOLATION_TYPE_PROPERTY(_itp_serendip_quadrangle_8,
                                          _itk_lagrangian, 8, 2);
AKANTU_DEFINE_INTERPOLATION_TYPE_PROPERTY(_itp_lagrange_hexahedron_8,
                                          _itk_lagrangian, 8, 3);
AKANTU_DEFINE_INTERPOLATION_TYPE_PROPERTY(_itp_serendip_hexahedron_20,
                                          _itk_lagrangian, 20, 3);
AKANTU_DEFINE_INTERPOLATION_TYPE_PROPERTY(_itp_lagrange_pentahedron_6,
                                          _itk_lagrangian, 6, 3);
AKANTU_DEFINE_INTERPOLATION_TYPE_PROPERTY(_itp_lagrange_pentahedron_15,
                                          _itk_lagrangian, 15, 3);

} // namespace akantu
#endif /* __AKANTU_INTERPOLATION_ELEMENT_TMPL_HH__ */
