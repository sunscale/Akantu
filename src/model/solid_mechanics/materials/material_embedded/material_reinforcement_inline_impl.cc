/**
 * @file   material_reinforcement_inline_impl.cc
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Thu Mar 12 2015
 * @date last modification: Thu Mar 12 2015
 *
 * @brief  Reinforcement material
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

/**
 * The structure of the directing cosines matrix is :
 * \f{eqnarray*}{
 *  C_{1,\cdot} & = & (l^2, m^2, n^2, lm, mn, ln) \\
 *  C_{i,j} & = & 0
 * \f}
 *
 * with :
 * \f[
 * (l, m, n) = \frac{1}{\|\frac{\mathrm{d}\vec{r}(s)}{\mathrm{d}s}\|} \cdot \frac{\mathrm{d}\vec{r}(s)}{\mathrm{d}s}
 * \f]
 */
template<UInt dim>
inline void MaterialReinforcement<dim>::computeDirectingCosinesOnQuad(const Matrix<Real> & nodes,
                                                                 Matrix<Real> & cosines) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(nodes.cols() == 2, "Higher order reinforcement elements not implemented");

  const Vector<Real> & a = nodes(0), b = nodes(1);

  cosines.clear();

  Real sq_length = 0.;
  for (UInt i = 0 ; i < dim ; i++) {
    sq_length += Math::pow<2, Real>(b(i) - a(i));
  }

  // Fill the first row of cosine matrix
  for (UInt i = 0 ; i < dim ; i++) {
    cosines(0, i) = Math::pow<2, Real>(b(i) - a(i)) / sq_length;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template<UInt dim>
inline void MaterialReinforcement<dim>::stressTensorToVoigtVector(const Matrix<Real> & tensor,
                                                                  Vector<Real> & vector) {
  AKANTU_DEBUG_IN();
  
  for (UInt i = 0; i < dim; i++) {
    vector(i) = tensor(i, i);
  }

  if (dim == 2) {
    vector(2) = tensor(0, 1);
  } else if (dim == 3) {
    vector(3) = tensor(0, 1);
    vector(4) = tensor(0, 2);
    vector(5) = tensor(1, 2);
  }
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template<UInt dim>
inline void MaterialReinforcement<dim>::strainTensorToVoigtVector(const Matrix<Real> & tensor,
                                                                  Vector<Real> & vector) {
  AKANTU_DEBUG_IN();
  
  for (UInt i = 0; i < dim; i++) {
    vector(i) = tensor(i, i);
  }

  if (dim == 2) {
    vector(2) = 2 * tensor(0, 1);
  } else if (dim == 3) {
    vector(3) = 2 * tensor(0, 1);
    vector(4) = 2 * tensor(0, 2);
    vector(5) = 2 * tensor(1, 2);
  }
  
  AKANTU_DEBUG_OUT();
}
