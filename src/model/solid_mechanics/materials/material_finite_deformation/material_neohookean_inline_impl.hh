/**
 * @file   material_neohookean_inline_impl.hh
 *
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 *
 * @date creation: Mon Apr 08 2013
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  Implementation of the inline functions of the material elastic
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
#include "material_neohookean.hh"
/* -------------------------------------------------------------------------- */
#include <cmath>
#include <iostream>
#include <utility>
/* -------------------------------------------------------------------------- */

namespace akantu {
/* -------------------------------------------------------------------------- */
template <UInt dim>
inline void MaterialNeohookean<dim>::computeDeltaStressOnQuad(
    __attribute__((unused)) const Matrix<Real> & grad_u,
    __attribute__((unused)) const Matrix<Real> & grad_delta_u,
    __attribute__((unused)) Matrix<Real> & delta_S) {}

//! computes the second piola kirchhoff stress, called S
template <UInt dim>
inline void MaterialNeohookean<dim>::computeStressOnQuad(Matrix<Real> & grad_u,
                                                         Matrix<Real> & S,
                                                         const Real & C33) {
  // Neo hookean book
  Matrix<Real> F(dim, dim);
  Matrix<Real> C(dim, dim);      // Right green
  Matrix<Real> Cminus(dim, dim); // Right green

  this->template gradUToF<dim>(grad_u, F);
  this->rightCauchy(F, C);
  Real J = F.det() * sqrt(C33); // the term  sqrt(C33) corresponds to the off
                                // plane strain (2D plane stress)
  //  std::cout<<"det(F) -> "<<J<<std::endl;
  Cminus.inverse(C);

  for (UInt i = 0; i < dim; ++i)
    for (UInt j = 0; j < dim; ++j)
      S(i, j) = (i == j) * mu + (lambda * log(J) - mu) * Cminus(i, j);
}

/* -------------------------------------------------------------------------- */
class C33_NR : public Math::NewtonRaphsonFunctor {
public:
  C33_NR(std::string name, const Real & lambda, const Real & mu,
         const Matrix<Real> & C)
      : NewtonRaphsonFunctor(std::move(name)), lambda(lambda), mu(mu), C(C) {}

  inline Real f(Real x) const override {
    return (this->lambda / 2. *
                (std::log(x) + std::log(this->C(0, 0) * this->C(1, 1) -
                                        Math::pow<2>(this->C(0, 1)))) +
            this->mu * (x - 1.));
  }

  inline Real f_prime(Real x) const override {
    AKANTU_DEBUG_ASSERT(std::abs(x) > Math::getTolerance(),
                        "x is zero (x should be the off plane right Cauchy"
                            << " measure in this function so you made a mistake"
                            << " somewhere else that lead to a zero here!!!");
    return (this->lambda / (2. * x) + this->mu);
  }

private:
  const Real & lambda;
  const Real & mu;
  const Matrix<Real> & C;
};

/* -------------------------------------------------------------------------- */
template <UInt dim>
inline void MaterialNeohookean<dim>::computeThirdAxisDeformationOnQuad(
    Matrix<Real> & grad_u, Real & c33_value) {
  // Neo hookean book
  Matrix<Real> F(dim, dim);
  Matrix<Real> C(dim, dim); // Right green

  this->template gradUToF<dim>(grad_u, F);
  this->rightCauchy(F, C);

  Math::NewtonRaphson nr(1e-5, 100);
  c33_value = nr.solve(
      C33_NR("Neohookean_plan_stress", this->lambda, this->mu, C), c33_value);
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
inline void
MaterialNeohookean<dim>::computePiolaKirchhoffOnQuad(const Matrix<Real> & E,
                                                     Matrix<Real> & S) {

  Real trace = E.trace(); /// \f$ trace = (\nabla u)_{kk} \f$

  /// \f$ \sigma_{ij} = \lambda * (\nabla u)_{kk} * \delta_{ij} + \mu * (\nabla
  /// u_{ij} + \nabla u_{ji}) \f$
  for (UInt i = 0; i < dim; ++i)
    for (UInt j = 0; j < dim; ++j)
      S(i, j) = (i == j) * lambda * trace + 2.0 * mu * E(i, j);
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
inline void MaterialNeohookean<dim>::computeFirstPiolaKirchhoffOnQuad(
    const Matrix<Real> & grad_u, const Matrix<Real> & S, Matrix<Real> & P) {

  Matrix<Real> F(dim, dim);
  Matrix<Real> C(dim, dim); // Right green

  this->template gradUToF<dim>(grad_u, F);

  // first Piola-Kirchhoff is computed as the product of the deformation
  // gracient
  // tensor and the second Piola-Kirchhoff stress tensor

  P = F * S;
}

/**************************************************************************************/
/*  Computation of the potential energy for a this neo hookean material */
template <UInt dim>
inline void MaterialNeohookean<dim>::computePotentialEnergyOnQuad(
    const Matrix<Real> & grad_u, Real & epot) {
  Matrix<Real> F(dim, dim);
  Matrix<Real> C(dim, dim); // Right green

  this->template gradUToF<dim>(grad_u, F);
  this->rightCauchy(F, C);
  Real J = F.det();
  //  std::cout<<"det(F) -> "<<J<<std::endl;

  epot =
      0.5 * lambda * pow(log(J), 2.) + mu * (-log(J) + 0.5 * (C.trace() - dim));
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
inline void MaterialNeohookean<dim>::computeTangentModuliOnQuad(
    Matrix<Real> & tangent, Matrix<Real> & grad_u, const Real & C33) {

  // Neo hookean book
  UInt cols = tangent.cols();
  UInt rows = tangent.rows();
  Matrix<Real> F(dim, dim);
  Matrix<Real> C(dim, dim);
  Matrix<Real> Cminus(dim, dim);
  this->template gradUToF<dim>(grad_u, F);
  this->rightCauchy(F, C);
  Real J = F.det() * sqrt(C33);
  //  std::cout<<"det(F) -> "<<J<<std::endl;
  Cminus.inverse(C);

  for (UInt m = 0; m < rows; m++) {
    UInt i = VoigtHelper<dim>::vec[m][0];
    UInt j = VoigtHelper<dim>::vec[m][1];
    for (UInt n = 0; n < cols; n++) {
      UInt k = VoigtHelper<dim>::vec[n][0];
      UInt l = VoigtHelper<dim>::vec[n][1];

      // book belytchko
      tangent(m, n) = lambda * Cminus(i, j) * Cminus(k, l) +
                      (mu - lambda * log(J)) * (Cminus(i, k) * Cminus(j, l) +
                                                Cminus(i, l) * Cminus(k, j));
    }
  }
}

/* -------------------------------------------------------------------------- */
} // namespace akantu
