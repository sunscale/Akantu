/**
 * @file   material_elastic_linear_anisotropic.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Till Junge <till.junge@epfl.ch>
 * @author Enrico Milanese <enrico.milanese@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Sep 25 2013
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Anisotropic elastic material
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

/* -------------------------------------------------------------------------- */
#include "material_elastic_linear_anisotropic.hh"
#include "solid_mechanics_model.hh"
#include <algorithm>
#include <sstream>

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt dim>
MaterialElasticLinearAnisotropic<dim>::MaterialElasticLinearAnisotropic(
    SolidMechanicsModel & model, const ID & id, bool symmetric)
    : Material(model, id), rot_mat(dim, dim), Cprime(dim * dim, dim * dim),
      C(voigt_h::size, voigt_h::size), eigC(voigt_h::size),
      symmetric(symmetric), alpha(0), was_stiffness_assembled(false) {
  AKANTU_DEBUG_IN();

  this->dir_vecs.push_back(std::make_unique<Vector<Real>>(dim));
  (*this->dir_vecs.back())[0] = 1.;
  this->registerParam("n1", *(this->dir_vecs.back()), _pat_parsmod,
                      "Direction of main material axis");

  if (dim > 1) {
    this->dir_vecs.push_back(std::make_unique<Vector<Real>>(dim));
    (*this->dir_vecs.back())[1] = 1.;
    this->registerParam("n2", *(this->dir_vecs.back()), _pat_parsmod,
                        "Direction of secondary material axis");
  }

  if (dim > 2) {
    this->dir_vecs.push_back(std::make_unique<Vector<Real>>(dim));
    (*this->dir_vecs.back())[2] = 1.;
    this->registerParam("n3", *(this->dir_vecs.back()), _pat_parsmod,
                        "Direction of tertiary material axis");
  }

  for (UInt i = 0; i < voigt_h::size; ++i) {
    UInt start = 0;
    if (this->symmetric) {
      start = i;
    }
    for (UInt j = start; j < voigt_h::size; ++j) {
      std::stringstream param("C");
      param << "C" << i + 1 << j + 1;
      this->registerParam(param.str(), this->Cprime(i, j), Real(0.),
                          _pat_parsmod, "Coefficient " + param.str());
    }
  }

  this->registerParam("alpha", this->alpha, _pat_parsmod,
                      "Proportion of viscous stress");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim> void MaterialElasticLinearAnisotropic<dim>::initMaterial() {
  AKANTU_DEBUG_IN();
  Material::initMaterial();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void MaterialElasticLinearAnisotropic<dim>::updateInternalParameters() {

  Material::updateInternalParameters();
  if (this->symmetric) {
    for (UInt i = 0; i < voigt_h::size; ++i) {
      for (UInt j = i + 1; j < voigt_h::size; ++j) {
        this->Cprime(j, i) = this->Cprime(i, j);
      }
    }
  }
  this->rotateCprime();
  this->C.eig(this->eigC);

  this->was_stiffness_assembled = false;
}

/* -------------------------------------------------------------------------- */
template <UInt Dim> void MaterialElasticLinearAnisotropic<Dim>::rotateCprime() {
  // start by filling the empty parts fo Cprime
  UInt diff = Dim * Dim - voigt_h::size;
  for (UInt i = voigt_h::size; i < Dim * Dim; ++i) {
    for (UInt j = 0; j < Dim * Dim; ++j) {
      this->Cprime(i, j) = this->Cprime(i - diff, j);
    }
  }
  for (UInt i = 0; i < Dim * Dim; ++i) {
    for (UInt j = voigt_h::size; j < Dim * Dim; ++j) {
      this->Cprime(i, j) = this->Cprime(i, j - diff);
    }
  }
  // construction of rotator tensor
  // normalise rotation matrix
  for (UInt j = 0; j < Dim; ++j) {
    Vector<Real> rot_vec = this->rot_mat(j);
    rot_vec = *this->dir_vecs[j];
    rot_vec.normalize();
  }

  // make sure the vectors form a right-handed base
  Vector<Real> test_axis(3);
  Vector<Real> v1(3);
  Vector<Real> v2(3);
  Vector<Real> v3(3, 0.);

  if (Dim == 2) {
    for (UInt i = 0; i < Dim; ++i) {
      v1[i] = this->rot_mat(0, i);
      v2[i] = this->rot_mat(1, i);
    }

    v3.crossProduct(v1, v2);
    if (v3.norm() < 8 * std::numeric_limits<Real>::epsilon()) {
      AKANTU_ERROR("The axis vectors parallel.");
    }

    v3.normalize();
  } else if (Dim == 3) {
    v1 = this->rot_mat(0);
    v2 = this->rot_mat(1);
    v3 = this->rot_mat(2);
  }

  test_axis.crossProduct(v1, v2);
  test_axis -= v3;
  if (test_axis.norm() > 8 * std::numeric_limits<Real>::epsilon()) {
    AKANTU_ERROR("The axis vectors do not form a right-handed coordinate "
                 << "system. I. e., ||n1 x n2 - n3|| should be zero, but "
                 << "it is " << test_axis.norm() << ".");
  }

  // create the rotator and the reverse rotator
  Matrix<Real> rotator(Dim * Dim, Dim * Dim);
  Matrix<Real> revrotor(Dim * Dim, Dim * Dim);
  for (UInt i = 0; i < Dim; ++i) {
    for (UInt j = 0; j < Dim; ++j) {
      for (UInt k = 0; k < Dim; ++k) {
        for (UInt l = 0; l < Dim; ++l) {
          UInt I = voigt_h::mat[i][j];
          UInt J = voigt_h::mat[k][l];
          rotator(I, J) = this->rot_mat(k, i) * this->rot_mat(l, j);
          revrotor(I, J) = this->rot_mat(i, k) * this->rot_mat(j, l);
        }
      }
    }
  }

  // create the full rotated matrix
  Matrix<Real> Cfull(Dim * Dim, Dim * Dim);
  Cfull = rotator * Cprime * revrotor;

  for (UInt i = 0; i < voigt_h::size; ++i) {
    for (UInt j = 0; j < voigt_h::size; ++j) {
      this->C(i, j) = Cfull(i, j);
    }
  }
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void MaterialElasticLinearAnisotropic<dim>::computeStress(
    ElementType el_type, GhostType ghost_type) {
  // Wikipedia convention:
  // 2*eps_ij (i!=j) = voigt_eps_I
  // http://en.wikipedia.org/wiki/Voigt_notation
  AKANTU_DEBUG_IN();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  this->computeStressOnQuad(grad_u, sigma);

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void MaterialElasticLinearAnisotropic<dim>::computeTangentModuli(
    ElementType el_type, Array<Real> & tangent_matrix,
    GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_matrix);

  this->computeTangentModuliOnQuad(tangent);

  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END;

  this->was_stiffness_assembled = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void MaterialElasticLinearAnisotropic<dim>::computePotentialEnergy(
    ElementType el_type) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(!this->finite_deformation,
                      "finite deformation not possible in material anisotropic "
                      "(TO BE IMPLEMENTED)");

  Array<Real>::scalar_iterator epot =
      this->potential_energy(el_type, _not_ghost).begin();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, _not_ghost);

  computePotentialEnergyOnQuad(grad_u, sigma, *epot);
  ++epot;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
Real MaterialElasticLinearAnisotropic<dim>::getCelerity(
    __attribute__((unused)) const Element & element) const {
  return std::sqrt(this->eigC(0) / rho);
}

/* -------------------------------------------------------------------------- */

INSTANTIATE_MATERIAL(elastic_anisotropic, MaterialElasticLinearAnisotropic);

} // namespace akantu
