/**
 * @file   weight_function_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Apr 13 2012
 * @date last modification: Thu Jun 05 2014
 *
 * @brief  implementation of the weight function classes
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

/* -------------------------------------------------------------------------- */
/* Stress based weight function                                               */
/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
StressBasedWeightFunction<spatial_dimension>::StressBasedWeightFunction(Material & material) :
  BaseWeightFunction<spatial_dimension>(material, "stress_based"),
  stress_diag("stress_diag", material), selected_stress_diag(NULL),
  stress_base("stress_base", material), selected_stress_base(NULL),
  characteristic_size("lc", material),  selected_characteristic_size(NULL) {

  this->registerParam("ft", this->ft, 0., _pat_parsable, "Tensile strength");
  stress_diag.initialize(spatial_dimension);
  stress_base.initialize(spatial_dimension * spatial_dimension);
  characteristic_size.initialize(1);
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void StressBasedWeightFunction<spatial_dimension>::init() {
  const Mesh & mesh = this->material.getModel().getFEEngine().getMesh();
  for (UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = GhostType(g);
    Mesh::type_iterator it = mesh.firstType(spatial_dimension, gt);
    Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, gt);
    for(; it != last_type; ++it) {
      UInt nb_quadrature_points =
	this->material.getModel().getFEEngine().getNbQuadraturePoints(*it, gt);
      const Array<UInt> & element_filter = this->material.getElementFilter(*it, gt);
      UInt nb_element = element_filter.getSize();

      Array<Real> ones(nb_element*nb_quadrature_points, 1, 1.);
      Array<Real> & lc = characteristic_size(*it, gt);
      this->material.getModel().getFEEngine().integrateOnQuadraturePoints(ones,
								     lc,
								     1,
								     *it,
								     gt,
								     element_filter);

      for (UInt q = 0;  q < nb_quadrature_points * nb_element; q++) {
	lc(q) = pow(lc(q), 1./ Real(spatial_dimension));
      }
    }
  }
}


/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void StressBasedWeightFunction<spatial_dimension>::updatePrincipalStress(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  const Mesh & mesh = this->material.getModel().getFEEngine().getMesh();

  Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type);
  for(; it != last_type; ++it) {
    Array<Real>::const_matrix_iterator sigma =
      this->material.getStress(*it, ghost_type).begin(spatial_dimension, spatial_dimension);
    Array<Real>::vector_iterator eigenvalues =
      stress_diag(*it, ghost_type).begin(spatial_dimension);
    Array<Real>::vector_iterator eigenvalues_end =
      stress_diag(*it, ghost_type).end(spatial_dimension);
    Array<Real>::matrix_iterator eigenvector =
      stress_base(*it, ghost_type).begin(spatial_dimension, spatial_dimension);

#ifndef __trick__
    Array<Real>::iterator<Real> cl = characteristic_size(*it, ghost_type).begin();
#endif
    UInt q = 0;
    for(;eigenvalues != eigenvalues_end; ++sigma, ++eigenvalues, ++eigenvector, ++cl, ++q) {
      sigma->eig(*eigenvalues, *eigenvector);
      *eigenvalues /= ft;
#ifndef __trick__
      // specify a lower bound for principal stress based on the size of the element
      for (UInt i = 0; i < spatial_dimension; ++i) {
        (*eigenvalues)(i) = std::max(*cl / this->R, (*eigenvalues)(i));
      }
#endif
    }
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline void StressBasedWeightFunction<spatial_dimension>::selectType(ElementType type1,
								     GhostType ghost_type1,
								     ElementType type2,
								     GhostType ghost_type2) {
  selected_stress_diag = &stress_diag(type2, ghost_type2);
  selected_stress_base = &stress_base(type2, ghost_type2);

  selected_characteristic_size = &characteristic_size(type1, ghost_type1);
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline Real StressBasedWeightFunction<spatial_dimension>::operator()(Real r,
								     QuadraturePoint & q1,
								     QuadraturePoint & q2) {
  Real zero = std::numeric_limits<Real>::epsilon();

  if(r < zero) return 1.; // means x and s are the same points

  const Vector<Real> & x = q1.getPosition();
  const Vector<Real> & s = q2.getPosition();

  Vector<Real> eigs =
    selected_stress_diag->begin(spatial_dimension)[q2.global_num];

  Matrix<Real> eigenvects =
    selected_stress_base->begin(spatial_dimension, spatial_dimension)[q2.global_num];

  Real min_rho_lc = selected_characteristic_size->begin()[q1.global_num];

  Vector<Real> x_s(spatial_dimension);
  x_s  = x;
  x_s -= s;

  Real rho_2 = computeRhoSquare(r, eigs, eigenvects, x_s);

  Real rho_lc_2 = std::max(this->R2 * rho_2, min_rho_lc*min_rho_lc);

  // Real w = std::max(0., 1. - r*r / rho_lc_2);
  // w = w*w;
  Real w = exp(- 2*2*r*r / rho_lc_2);
  return w;
}

/* -------------------------------------------------------------------------- */
template<>
inline Real StressBasedWeightFunction<1>::computeRhoSquare(__attribute__ ((unused)) Real r,
							   Vector<Real> & eigs,
							   __attribute__ ((unused)) Matrix<Real> & eigenvects,
							   __attribute__ ((unused)) Vector<Real> & x_s) {
  return eigs[0];
}

/* -------------------------------------------------------------------------- */
template<>
inline Real StressBasedWeightFunction<2>::computeRhoSquare(__attribute__ ((unused)) Real r,
							   Vector<Real> & eigs,
							   Matrix<Real> & eigenvects,
							   Vector<Real> & x_s) {
  Vector<Real> u1(eigenvects.storage(), 2);
  Real cos_t = x_s.dot(u1) / (x_s.norm() * u1.norm());

  Real cos_t_2;
  Real sin_t_2;

  Real sigma1_2 = eigs[0]*eigs[0];
  Real sigma2_2 = eigs[1]*eigs[1];

#ifdef __trick__
  Real zero = std::numeric_limits<Real>::epsilon();
  if(std::abs(cos_t) < zero) {
    cos_t_2 = 0;
    sin_t_2 = 1;
  } else {
    cos_t_2 = cos_t * cos_t;
    sin_t_2 = (1 - cos_t_2);
  }

  Real rhop1 = std::max(0., cos_t_2 / sigma1_2);
  Real rhop2 = std::max(0., sin_t_2 / sigma2_2);
#else
  cos_t_2 = cos_t * cos_t;
  sin_t_2 = (1 - cos_t_2);

  Real rhop1 = cos_t_2 / sigma1_2;
  Real rhop2 = sin_t_2 / sigma2_2;
#endif

  return 1./ (rhop1 + rhop2);
}

/* -------------------------------------------------------------------------- */
template<>
inline Real StressBasedWeightFunction<3>::computeRhoSquare(Real r,
							   Vector<Real> & eigs,
							   Matrix<Real> & eigenvects,
							   Vector<Real> & x_s) {
  Vector<Real> u1(eigenvects.storage() + 0*3, 3);
//Vector<Real> u2(eigenvects.storage() + 1*3, 3);
  Vector<Real> u3(eigenvects.storage() + 2*3, 3);

  Real zero = std::numeric_limits<Real>::epsilon();

  Vector<Real> tmp(3);
  tmp.crossProduct(x_s, u3);

  Vector<Real> u3_C_x_s_C_u3(3);
  u3_C_x_s_C_u3.crossProduct(u3, tmp);

  Real norm_u3_C_x_s_C_u3 = u3_C_x_s_C_u3.norm();
  Real cos_t = 0.;
  if(std::abs(norm_u3_C_x_s_C_u3) > zero) {
    Real inv_norm_u3_C_x_s_C_u3 = 1. / norm_u3_C_x_s_C_u3;
    cos_t = u1.dot(u3_C_x_s_C_u3) * inv_norm_u3_C_x_s_C_u3;
  }

  Real cos_p = u3.dot(x_s) / r;

  Real cos_t_2;
  Real sin_t_2;
  Real cos_p_2;
  Real sin_p_2;

  Real sigma1_2 = eigs[0]*eigs[0];
  Real sigma2_2 = eigs[1]*eigs[1];
  Real sigma3_2 = eigs[2]*eigs[2];

#ifdef __trick__
  if(std::abs(cos_t) < zero) {
    cos_t_2 = 0;
    sin_t_2 = 1;
  } else {
    cos_t_2 = cos_t * cos_t;
    sin_t_2 = (1 - cos_t_2);
  }

  if(std::abs(cos_p) < zero) {
    cos_p_2 = 0;
    sin_p_2 = 1;
  } else {
    cos_p_2 = cos_p * cos_p;
    sin_p_2 = (1 - cos_p_2);
  }

  Real rhop1 = std::max(0., sin_p_2 * cos_t_2 / sigma1_2);
  Real rhop2 = std::max(0., sin_p_2 * sin_t_2 / sigma2_2);
  Real rhop3 = std::max(0., cos_p_2 / sigma3_2);
#else
  cos_t_2 = cos_t * cos_t;
  sin_t_2 = (1 - cos_t_2);

  cos_p_2 = cos_p * cos_p;
  sin_p_2 = (1 - cos_p_2);

  Real rhop1 = sin_p_2 * cos_t_2 / sigma1_2;
  Real rhop2 = sin_p_2 * sin_t_2 / sigma2_2;
  Real rhop3 = cos_p_2 / sigma3_2;
#endif

  return 1./ (rhop1 + rhop2 + rhop3);
}
