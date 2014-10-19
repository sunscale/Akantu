/**
 * @file   integration_element.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Jan 16 2013
 * @date last modification: Mon Jul 07 2014
 *
 * @brief  Definition of the intagration constants
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
#include "aka_common.hh"
#include "element_class.hh"

using std::sqrt;

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/* Points                                                                     */
/* -------------------------------------------------------------------------- */
template<> UInt GaussIntegrationTypeData<_git_point, 1>::nb_quadrature_points = 1;
template<> Real GaussIntegrationTypeData<_git_point, 1>::quad_positions[]     = {0};
template<> Real GaussIntegrationTypeData<_git_point, 1>::quad_weights[]       = {1.};

/* -------------------------------------------------------------------------- */
/* Segments                                                                   */
/* -------------------------------------------------------------------------- */
template<> UInt GaussIntegrationTypeData<_git_segment, 1>::nb_quadrature_points = 1;
template<> Real GaussIntegrationTypeData<_git_segment, 1>::quad_positions[]     = {0.};
template<> Real GaussIntegrationTypeData<_git_segment, 1>::quad_weights[]       = {2.};
/* -------------------------------------------------------------------------- */
template<> UInt GaussIntegrationTypeData<_git_segment, 2>::nb_quadrature_points = 2;
template<> Real GaussIntegrationTypeData<_git_segment, 2>::quad_positions[]     = {-1./sqrt(3.), 1./sqrt(3.)};
template<> Real GaussIntegrationTypeData<_git_segment, 2>::quad_weights[]       = {1., 1.};
/* -------------------------------------------------------------------------- */
template<> UInt GaussIntegrationTypeData<_git_segment, 3>::nb_quadrature_points = 3;
template<> Real GaussIntegrationTypeData<_git_segment, 3>::quad_positions[]     = {-sqrt(3./5.), 0., sqrt(3./5.)};
template<> Real GaussIntegrationTypeData<_git_segment, 3>::quad_weights[]       = {5./9., 8./9., 5./9.};
/* -------------------------------------------------------------------------- */
template<> UInt GaussIntegrationTypeData<_git_segment, 4>::nb_quadrature_points = 4;
template<> Real GaussIntegrationTypeData<_git_segment, 4>::quad_positions[]     = {-sqrt((3. + 2.*sqrt(6./5.))/7.),
										   -sqrt((3. - 2.*sqrt(6./5.))/7.), 
										   sqrt((3. - 2.*sqrt(6./5.))/7.), 
										   sqrt((3. + 2.*sqrt(6./5.))/7.)};
template<> Real GaussIntegrationTypeData<_git_segment, 4>::quad_weights[]       = {(18. - sqrt(30.))/36., 
										   (18. + sqrt(30.))/36., 
										   (18. + sqrt(30.))/36., 
										   (18. - sqrt(30.))/36.};
/* -------------------------------------------------------------------------- */
/* Triangles                                                                  */
/* -------------------------------------------------------------------------- */
template<> UInt GaussIntegrationTypeData<_git_triangle, 1>::nb_quadrature_points = 1;
template<> Real GaussIntegrationTypeData<_git_triangle, 1>::quad_positions[]     = {1./3., 1./3.};
template<> Real GaussIntegrationTypeData<_git_triangle, 1>::quad_weights[]       = {1./2.};
/* -------------------------------------------------------------------------- */
template<> UInt GaussIntegrationTypeData<_git_triangle, 2>::nb_quadrature_points = 3;
template<> Real GaussIntegrationTypeData<_git_triangle, 2>::quad_positions[]     = {1./6., 1./6.,
										    2./3., 1./6.,
										    1./6., 2./3.};
template<> Real GaussIntegrationTypeData<_git_triangle, 2>::quad_weights[]       = {1./6., 1./6., 1./6.};
/* -------------------------------------------------------------------------- */
template<> UInt GaussIntegrationTypeData<_git_triangle, 3>::nb_quadrature_points = 4;
template<> Real GaussIntegrationTypeData<_git_triangle, 3>::quad_positions[]     = {1./5., 1./5.,
										    3./5., 1./5.,
										    1./5., 3./5.,
										    1./3., 1./3.};
template<> Real GaussIntegrationTypeData<_git_triangle, 3>::quad_weights[]       = {25./(24.*4.), 25./(24.*4.), 25./(24.*4.), -27/(24.*4.)};

/* -------------------------------------------------------------------------- */
/* Tetrahedrons                                                               */
/* -------------------------------------------------------------------------- */
template<> UInt GaussIntegrationTypeData<_git_tetrahedron, 1>::nb_quadrature_points = 1;
template<> Real GaussIntegrationTypeData<_git_tetrahedron, 1>::quad_positions[]     = {1./4., 1./4., 1./4.};
template<> Real GaussIntegrationTypeData<_git_tetrahedron, 1>::quad_weights[]       = {1./6.};
/* -------------------------------------------------------------------------- */
static const Real tet_2_a = (5. -    std::sqrt(5.))/20.;
static const Real tet_2_b = (5. + 3.*std::sqrt(5.))/20.;
template<> UInt GaussIntegrationTypeData<_git_tetrahedron, 2>::nb_quadrature_points = 4;
template<> Real GaussIntegrationTypeData<_git_tetrahedron, 2>::quad_positions[]     = {tet_2_a, tet_2_a, tet_2_a,
										       tet_2_b, tet_2_a, tet_2_a,
										       tet_2_a, tet_2_b, tet_2_a,
										       tet_2_a, tet_2_a, tet_2_b};
template<> Real GaussIntegrationTypeData<_git_tetrahedron, 2>::quad_weights[]       = {1./24., 1./24., 1./24., 1./24.};


/* -------------------------------------------------------------------------- */
/* Tetrahedrons                                                               */
/* -------------------------------------------------------------------------- */
template<> UInt GaussIntegrationTypeData<_git_pentahedron, 1>::nb_quadrature_points = 6;
template<> Real GaussIntegrationTypeData<_git_pentahedron, 1>::quad_positions[]     = {-1./std::sqrt(3.), 0.5, 0.5,
										       -1./std::sqrt(3.), 0. , 0.5,
										       -1./std::sqrt(3.), 0.5, 0.,
  										        1./std::sqrt(3.), 0.5, 0.5,
										        1./std::sqrt(3.), 0. , 0.5,
										        1./std::sqrt(3.), 0.5 ,0.};
template<> Real GaussIntegrationTypeData<_git_pentahedron, 1>::quad_weights[]       = {1./6, 1./6, 1./6, 
										       1./6, 1./6, 1./6};

__END_AKANTU__
