/**
 * @file   cohesive_element.hh
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Feb 03 2012
 * @date last modification: Fri Jun 13 2014
 *
 * @brief  Cohesive element class
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "element_class.hh"

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_COHESIVE_ELEMENT_HH__
#define __AKANTU_COHESIVE_ELEMENT_HH__

__BEGIN_AKANTU__

AKANTU_DEFINE_ELEMENT_CLASS_PROPERTY(_cohesive_2d_4, _gt_cohesive_2d_4, _itp_lagrange_segment_2,
				     _ek_cohesive, 2,
				     _git_segment, 1);

AKANTU_DEFINE_ELEMENT_CLASS_PROPERTY(_cohesive_2d_6, _gt_cohesive_2d_6, _itp_lagrange_segment_3,
				     _ek_cohesive, 2,
				     _git_segment, 2);

AKANTU_DEFINE_ELEMENT_CLASS_PROPERTY(_cohesive_1d_2, _gt_cohesive_1d_2, _itp_lagrange_point_1,
				     _ek_cohesive, 1,
				     _git_point, 1);

AKANTU_DEFINE_ELEMENT_CLASS_PROPERTY(_cohesive_3d_6, _gt_cohesive_3d_6, _itp_lagrange_triangle_3,
				     _ek_cohesive, 3,
				     _git_triangle, 1);

AKANTU_DEFINE_ELEMENT_CLASS_PROPERTY(_cohesive_3d_12, _gt_cohesive_3d_12, _itp_lagrange_triangle_6,
				     _ek_cohesive, 3,
				     _git_triangle, 2);


__END_AKANTU__

#endif /* __AKANTU_COHESIVE_ELEMENT_HH__ */
