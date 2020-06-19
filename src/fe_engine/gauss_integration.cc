/**
 * @file   gauss_integration.cc
 *
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 * @author Sébastien Hartmann <sebastien.hartmann@epfl.ch>
 * @author Thomas Menouillard <tmenouillard@stucky.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Feb 21 2013
 * @date last modification: Mon Jun 19 2017
 *
 * @brief  Definition of the integration constants, some of the value are taken
 * from r3.01.01 doc from Code Aster
 *
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
#include "aka_common.hh"
#include "element_class.hh"

using std::sqrt;

namespace akantu {

/* clang-format off */
/* -------------------------------------------------------------------------- */
/* Points                                                                     */
/* -------------------------------------------------------------------------- */
template<> Real GaussIntegrationTypeData<_git_point, 1>::quad_positions[]     = {0};
template<> Real GaussIntegrationTypeData<_git_point, 1>::quad_weights[]       = {1.};

/* -------------------------------------------------------------------------- */
/* Segments                                                                   */
/* -------------------------------------------------------------------------- */
template<> Real GaussIntegrationTypeData<_git_segment, 1>::quad_positions[]     = {0.};
template<> Real GaussIntegrationTypeData<_git_segment, 1>::quad_weights[]       = {2.};
/* -------------------------------------------------------------------------- */
template<> Real GaussIntegrationTypeData<_git_segment, 2>::quad_positions[]     = {-1./sqrt(3.), 1./sqrt(3.)};
template<> Real GaussIntegrationTypeData<_git_segment, 2>::quad_weights[]       = {1., 1.};
/* -------------------------------------------------------------------------- */
template<> Real GaussIntegrationTypeData<_git_segment, 3>::quad_positions[]     = {-sqrt(3./5.), 0., sqrt(3./5.)};
template<> Real GaussIntegrationTypeData<_git_segment, 3>::quad_weights[]       = {5./9., 8./9., 5./9.};
/* -------------------------------------------------------------------------- */
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
template<> Real GaussIntegrationTypeData<_git_triangle, 1>::quad_positions[]     = {1./3., 1./3.};
template<> Real GaussIntegrationTypeData<_git_triangle, 1>::quad_weights[]       = {1./2.};
/* -------------------------------------------------------------------------- */
template<> Real GaussIntegrationTypeData<_git_triangle, 3>::quad_positions[]     = {1./6., 1./6.,
                                                                                    2./3., 1./6.,
                                                                                    1./6., 2./3.};
template<> Real GaussIntegrationTypeData<_git_triangle, 3>::quad_weights[]       = {1./6., 1./6., 1./6.};
/* -------------------------------------------------------------------------- */
template<> Real GaussIntegrationTypeData<_git_triangle, 4>::quad_positions[]     = {1./5., 1./5.,
                                                                                    3./5., 1./5.,
                                                                                    1./5., 3./5.,
                                                                                    1./3., 1./3.};
template<> Real GaussIntegrationTypeData<_git_triangle, 4>::quad_weights[]       = {25./(24.*4.), 25./(24.*4.), 25./(24.*4.), -27/(24.*4.)};
/* -------------------------------------------------------------------------- */
/// Found those one in the TrigGaussRuleInfo from mathematica and matched them to the code aster values
/// http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.AppI.d/AFEM.AppI.pdf
static const Real tri_6_a = (8. - std::sqrt(10.) + std::sqrt(38.-44.*std::sqrt(2./5.)))/18.;
static const Real tri_6_b = (8. - std::sqrt(10.) - std::sqrt(38.-44.*std::sqrt(2./5.)))/18.;
static const Real tri_6_w1 = (620. - std::sqrt(213125. - 53320.*std::sqrt(10.)))/7440.;
static const Real tri_6_w2 = (620. + std::sqrt(213125. - 53320.*std::sqrt(10.)))/7440.;
template<> Real GaussIntegrationTypeData<_git_triangle, 6>::quad_positions[]     = {tri_6_b, tri_6_b,
                                                                                    1. - 2. * tri_6_b, tri_6_b,
                                                                                    tri_6_b, 1. - 2. * tri_6_b,
                                                                                    tri_6_a, 1. - 2. * tri_6_a,
                                                                                    tri_6_a, tri_6_a,
                                                                                    1. - 2. * tri_6_a, tri_6_a};
template<> Real GaussIntegrationTypeData<_git_triangle, 6>::quad_weights[]       = {tri_6_w1, tri_6_w1, tri_6_w1,
                                                                                    tri_6_w2, tri_6_w2, tri_6_w2};
/* -------------------------------------------------------------------------- */
static const Real tri_7_a  = (6. + std::sqrt(15.)) / 21.;
static const Real tri_7_b  = (6. - std::sqrt(15.)) / 21.;
static const Real tri_7_w1 = (155. + std::sqrt(15.))/2400.;
static const Real tri_7_w2 = (155. - std::sqrt(15.))/2400.;
template<> Real GaussIntegrationTypeData<_git_triangle, 7>::quad_positions[]     = {          1./3.,           1./3.,
                                                                                            tri_7_a,         tri_7_a,
                                                                                    1. - 2.*tri_7_a,         tri_7_a,
                                                                                            tri_7_a, 1. - 2.*tri_7_a,
                                                                                            tri_7_b,         tri_7_b,
                                                                                    1. - 2.*tri_7_b,         tri_7_b,
                                                                                            tri_7_b, 1. - 2.*tri_7_b};
template<> Real GaussIntegrationTypeData<_git_triangle, 7>::quad_weights[]       = {9./80.,
                                                                                    tri_7_w1, tri_7_w1, tri_7_w1,
                                                                                    tri_7_w2, tri_7_w2, tri_7_w2};

/* -------------------------------------------------------------------------- */


/* -------------------------------------------------------------------------- */
/* Tetrahedrons                                                               */
/* -------------------------------------------------------------------------- */
template<> Real GaussIntegrationTypeData<_git_tetrahedron, 1>::quad_positions[]     = {1./4., 1./4., 1./4.};
template<> Real GaussIntegrationTypeData<_git_tetrahedron, 1>::quad_weights[]       = {1./6.};
/* -------------------------------------------------------------------------- */
static const Real tet_4_a = (5. -    std::sqrt(5.))/20.;
static const Real tet_4_b = (5. + 3.*std::sqrt(5.))/20.;
template<> Real GaussIntegrationTypeData<_git_tetrahedron, 4>::quad_positions[]     = {tet_4_a, tet_4_a, tet_4_a,
                                                                                       tet_4_b, tet_4_a, tet_4_a,
                                                                                       tet_4_a, tet_4_b, tet_4_a,
                                                                                       tet_4_a, tet_4_a, tet_4_b};
template<> Real GaussIntegrationTypeData<_git_tetrahedron, 4>::quad_weights[]       = {1./24., 1./24., 1./24., 1./24.};
/* -------------------------------------------------------------------------- */
template<> Real GaussIntegrationTypeData<_git_tetrahedron, 5>::quad_positions[]     = {1./4., 1./4., 1./4.,
                                                                                       1./6., 1./6., 1./6.,
                                                                                       1./6., 1./6., 1./2.,
                                                                                       1./6., 1./2., 1./6.,
                                                                                       1./2., 1./6., 1./6.,};
template<> Real GaussIntegrationTypeData<_git_tetrahedron, 5>::quad_weights[]       = {-2./15., 3./40.,
                                                                                       3./40., 3./40.,
                                                                                       3./40.};
/* -------------------------------------------------------------------------- */
static const Real tet_15_a = (7. + std::sqrt(15.))/34.;
static const Real tet_15_b = (13. - 3. * std::sqrt(15.))/34.;
static const Real tet_15_c = (7. - std::sqrt(15.))/34.;
static const Real tet_15_d = (13. + 3. * std::sqrt(15.))/34.;
static const Real tet_15_e = (5. - std::sqrt(15.))/20.;
static const Real tet_15_f = (5. + std::sqrt(15.))/20.;
static const Real tet_15_w1 = (2665. - 14. * std::sqrt(15.))/226800.;
static const Real tet_15_w2 = (2665. + 14. * std::sqrt(15.))/226800.;
static const Real tet_15_w3 = 5./567.;
template<> Real GaussIntegrationTypeData<_git_tetrahedron, 15>::quad_positions[]    = {1./4., 1./4., 1./4.,
                                                                                       tet_15_a, tet_15_a, tet_15_a,
                                                                                       tet_15_a, tet_15_a, tet_15_b,
                                                                                       tet_15_a, tet_15_b, tet_15_a,
                                                                                       tet_15_b, tet_15_a, tet_15_a,

                                                                                       tet_15_c, tet_15_c, tet_15_c,
                                                                                       tet_15_c, tet_15_c, tet_15_d,
                                                                                       tet_15_c, tet_15_d, tet_15_c,
                                                                                       tet_15_d, tet_15_c, tet_15_c,

                                                                                       tet_15_e, tet_15_e, tet_15_f,
                                                                                       tet_15_e, tet_15_f, tet_15_e,
                                                                                       tet_15_f, tet_15_e, tet_15_e,

                                                                                       tet_15_e, tet_15_f, tet_15_f,
                                                                                       tet_15_f, tet_15_e, tet_15_f,
                                                                                       tet_15_f, tet_15_f, tet_15_e};
template<> Real GaussIntegrationTypeData<_git_tetrahedron, 15>::quad_weights[]      = {8./405.,
                                                                                       tet_15_w1, tet_15_w1, tet_15_w1, tet_15_w1,
                                                                                       tet_15_w2, tet_15_w2, tet_15_w2, tet_15_w2,
                                                                                       tet_15_w3, tet_15_w3, tet_15_w3, tet_15_w3, tet_15_w3, tet_15_w3};

/* -------------------------------------------------------------------------- */
/* Pentahedrons                                                               */
/* -------------------------------------------------------------------------- */
template<> Real GaussIntegrationTypeData<_git_pentahedron, 6>::quad_positions[]     = {-1./sqrt(3.), 0.5, 0.5,
                                                                                       -1./sqrt(3.), 0. , 0.5,
                                                                                       -1./sqrt(3.), 0.5, 0.,
                                                                                        1./sqrt(3.), 0.5, 0.5,
                                                                                        1./sqrt(3.), 0. , 0.5,
                                                                                        1./sqrt(3.), 0.5 ,0.};
template<> Real GaussIntegrationTypeData<_git_pentahedron, 6>::quad_weights[]       = {1./6., 1./6., 1./6.,
                                                                                       1./6., 1./6., 1./6.};
/* -------------------------------------------------------------------------- */
template<> Real GaussIntegrationTypeData<_git_pentahedron, 8>::quad_positions[]     = {-sqrt(3.)/3., 1./3., 1./3.,
                                                                                       -sqrt(3.)/3.,   0.6,   0.2,
                                                                                       -sqrt(3.)/3.,   0.2,   0.6,
                                                                                       -sqrt(3.)/3.,   0.2,   0.2,
                                                                                        sqrt(3.)/3., 1./3., 1./3.,
                                                                                        sqrt(3.)/3.,   0.6,   0.2,
                                                                                        sqrt(3.)/3.,   0.2,   0.6,
                                                                                        sqrt(3.)/3.,   0.2,   0.2};
template<> Real GaussIntegrationTypeData<_git_pentahedron, 8>::quad_weights[]       = {-27./96., 25./96., 25./96., 25./96.,
                                                                                       -27./96., 25./96., 25./96., 25./96.};
/* -------------------------------------------------------------------------- */
static const Real pent_21_x = std::sqrt(3./5.);
static const Real pent_21_a = (6. + std::sqrt(15.)) / 21.;
static const Real pent_21_b = (6. - std::sqrt(15.)) / 21.;

static const Real pent_21_w1_1 = 5./9.;
static const Real pent_21_w2_1 = 8./9.;
static const Real pent_21_w1_2 = (155. + std::sqrt(15.))/2400.;
static const Real pent_21_w2_2 = (155. - std::sqrt(15.))/2400.;
template<> Real GaussIntegrationTypeData<_git_pentahedron, 21>::quad_positions[]    = {- pent_21_x,             1./3.,             1./3.,
                                                                                       - pent_21_x,         pent_21_a,         pent_21_a,
                                                                                       - pent_21_x, 1. - 2.*pent_21_a,         pent_21_a,
                                                                                       - pent_21_x,         pent_21_a, 1. - 2.*pent_21_a,
                                                                                       - pent_21_x,         pent_21_b,         pent_21_b,
                                                                                       - pent_21_x, 1. - 2.*pent_21_b,         pent_21_b,
                                                                                       - pent_21_x,         pent_21_b, 1. - 2.*pent_21_b,
                                                                                                0.,             1./3.,             1./3.,
                                                                                                0.,         pent_21_a,         pent_21_a,
                                                                                                0., 1. - 2.*pent_21_a,         pent_21_a,
                                                                                                0.,         pent_21_a, 1. - 2.*pent_21_a,
                                                                                                0.,         pent_21_b,         pent_21_b,
                                                                                                0., 1. - 2.*pent_21_b,         pent_21_b,
                                                                                                0.,         pent_21_b, 1. - 2.*pent_21_b,
                                                                                         pent_21_x,             1./3.,             1./3.,
                                                                                         pent_21_x,         pent_21_a,         pent_21_a,
                                                                                         pent_21_x, 1. - 2.*pent_21_a,         pent_21_a,
                                                                                         pent_21_x,         pent_21_a, 1. - 2.*pent_21_a,
                                                                                         pent_21_x,         pent_21_b,         pent_21_b,
                                                                                         pent_21_x, 1. - 2.*pent_21_b,         pent_21_b,
                                                                                         pent_21_x,         pent_21_b, 1. - 2.*pent_21_b};
template<> Real GaussIntegrationTypeData<_git_pentahedron, 21>::quad_weights[]      = {pent_21_w1_1 * 9. / 80.,
                                                                                       pent_21_w1_1*pent_21_w1_2, pent_21_w1_1*pent_21_w1_2, pent_21_w1_1*pent_21_w1_2,
                                                                                       pent_21_w1_1*pent_21_w2_2, pent_21_w1_1*pent_21_w2_2, pent_21_w1_1*pent_21_w2_2,
                                                                                       pent_21_w2_1 * 9. / 80.,
                                                                                       pent_21_w1_2*pent_21_w1_2, pent_21_w1_2*pent_21_w1_2, pent_21_w1_2*pent_21_w1_2,
                                                                                       pent_21_w1_2*pent_21_w2_2, pent_21_w1_2*pent_21_w2_2, pent_21_w1_2*pent_21_w2_2,
                                                                                       pent_21_w1_1 * 9. / 80.,
                                                                                       pent_21_w1_1*pent_21_w1_2, pent_21_w1_1*pent_21_w1_2, pent_21_w1_1*pent_21_w1_2,
                                                                                       pent_21_w1_1*pent_21_w2_2, pent_21_w1_1*pent_21_w2_2, pent_21_w1_1*pent_21_w2_2};


} // akantu
