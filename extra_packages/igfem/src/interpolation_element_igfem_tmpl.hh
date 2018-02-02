/**
 * @file   element_class_igfem.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 *
 * @brief  Interpolation property description for IGFEM
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/* IGFEM elements                                                             */
/* -------------------------------------------------------------------------- */
#if defined(AKANTU_IGFEM)
AKANTU_DEFINE_IGFEM_INTERPOLATION_TYPE_PROPERTY(_itp_igfem_segment_3,         _itk_igfem,       3, 1, _itp_lagrange_segment_2, _itp_lagrange_segment_2, _itp_lagrange_segment_2);
AKANTU_DEFINE_IGFEM_INTERPOLATION_TYPE_PROPERTY(_itp_igfem_triangle_4,        _itk_igfem,       4, 2, _itp_lagrange_triangle_3, _itp_lagrange_triangle_3, _itp_lagrange_triangle_3);
AKANTU_DEFINE_IGFEM_INTERPOLATION_TYPE_PROPERTY(_itp_igfem_triangle_5,        _itk_igfem,       5, 2, _itp_lagrange_triangle_3, _itp_lagrange_triangle_3, _itp_lagrange_quadrangle_4);
#endif
