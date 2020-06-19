/**
 * @file   material_igfem_elastic_inline_impl.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  Implementation of the inline functions of the material igfem elastic
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
inline void MaterialIGFEMElastic<spatial_dimension>::computeStressOnQuad(
    const Matrix<Real> & grad_u, Matrix<Real> & sigma, const Real lambda,
    const Real mu) const {
  Real trace = grad_u.trace(); // trace = (\nabla u)_{kk}

  // \sigma_{ij} = \lambda * (\nabla u)_{kk} * \delta_{ij} + \mu * (\nabla
  // u_{ij} + \nabla u_{ji})
  for (UInt i = 0; i < spatial_dimension; ++i) {
    for (UInt j = 0; j < spatial_dimension; ++j) {
      sigma(i, j) =
          (i == j) * lambda * trace + mu * (grad_u(i, j) + grad_u(j, i));
    }
  }
}

/* -------------------------------------------------------------------------- */
template <>
inline void MaterialIGFEMElastic<1>::computeStressOnQuad(
    const Matrix<Real> & grad_u, Matrix<Real> & sigma, const Real lambda,
    const Real mu) const {
  sigma(0, 0) = 2 * mu * grad_u(0, 0);
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
inline void MaterialIGFEMElastic<spatial_dimension>::computeTangentModuliOnQuad(
    Matrix<Real> & tangent, const Real lambda, const Real mu) const {
  UInt n = tangent.cols();

  // Real Ep = E/((1+nu)*(1-2*nu));
  Real Miiii = lambda + 2 * mu;
  Real Miijj = lambda;
  Real Mijij = mu;

  if (spatial_dimension == 1)
    tangent(0, 0) = 3 * mu;
  else
    tangent(0, 0) = Miiii;

  // test of dimension should by optimized out by the compiler due to the
  // template
  if (spatial_dimension >= 2) {
    tangent(1, 1) = Miiii;
    tangent(0, 1) = Miijj;
    tangent(1, 0) = Miijj;

    tangent(n - 1, n - 1) = Mijij;
  }

  if (spatial_dimension == 3) {
    tangent(2, 2) = Miiii;
    tangent(0, 2) = Miijj;
    tangent(1, 2) = Miijj;
    tangent(2, 0) = Miijj;
    tangent(2, 1) = Miijj;

    tangent(3, 3) = Mijij;
    tangent(4, 4) = Mijij;
  }
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
inline void MaterialIGFEMElastic<dim>::computePotentialEnergyOnQuad(
    const Matrix<Real> & grad_u, const Matrix<Real> & sigma, Real & epot) {
  AKANTU_DEBUG_TO_IMPLEMENT();
  ///  epot = .5 * sigma.doubleDot(grad_u);
}

/* -------------------------------------------------------------------------- */
template <>
inline void MaterialIGFEMElastic<1>::computeTangentModuliOnQuad(
    Matrix<Real> & tangent, const Real lambda, const Real mu) const {
  tangent(0, 0) = 2 * mu;
}
