/**
 * @file   material_elastic_orthotropic.cc
 *
 * @author Till Junge <till.junge@epfl.ch>
 * @author Enrico Milanese <enrico.milanese@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Orthotropic elastic material
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_elastic_orthotropic.hh"
#include "solid_mechanics_model.hh"
#include <algorithm>

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt Dim>
MaterialElasticOrthotropic<Dim>::MaterialElasticOrthotropic(
    SolidMechanicsModel & model, const ID & id)
    : MaterialElasticLinearAnisotropic<Dim>(model, id) {
  AKANTU_DEBUG_IN();
  this->registerParam("E1", E1, Real(0.), _pat_parsmod, "Young's modulus (n1)");
  this->registerParam("E2", E2, Real(0.), _pat_parsmod, "Young's modulus (n2)");
  this->registerParam("nu12", nu12, Real(0.), _pat_parsmod,
                      "Poisson's ratio (12)");
  this->registerParam("G12", G12, Real(0.), _pat_parsmod, "Shear modulus (12)");

  if (Dim > 2) {
    this->registerParam("E3", E3, Real(0.), _pat_parsmod,
                        "Young's modulus (n3)");
    this->registerParam("nu13", nu13, Real(0.), _pat_parsmod,
                        "Poisson's ratio (13)");
    this->registerParam("nu23", nu23, Real(0.), _pat_parsmod,
                        "Poisson's ratio (23)");
    this->registerParam("G13", G13, Real(0.), _pat_parsmod,
                        "Shear modulus (13)");
    this->registerParam("G23", G23, Real(0.), _pat_parsmod,
                        "Shear modulus (23)");
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt Dim> void MaterialElasticOrthotropic<Dim>::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialElasticLinearAnisotropic<Dim>::initMaterial();

  AKANTU_DEBUG_ASSERT(not this->finite_deformation,
                      "finite deformation not possible in material orthotropic "
                      "(TO BE IMPLEMENTED)");

  updateInternalParameters();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt Dim>
void MaterialElasticOrthotropic<Dim>::updateInternalParameters() {

  this->C.zero();
  this->Cprime.zero();

  /* 1) construction of temporary material frame stiffness tensor------------ */
  // http://solidmechanics.org/Text/Chapter3_2/Chapter3_2.php#Sect3_2_13
  Real nu21 = nu12 * E2 / E1;
  Real nu31 = nu13 * E3 / E1;
  Real nu32 = nu23 * E3 / E2;

  // Full (i.e. dim^2 by dim^2) stiffness tensor in material frame
  if (Dim == 1) {
    AKANTU_ERROR("Dimensions 1 not implemented: makes no sense to have "
                 "orthotropy for 1D");
  }

  Real Gamma;

  if (Dim == 3) {
    Gamma = 1 / (1 - nu12 * nu21 - nu23 * nu32 - nu31 * nu13 -
                 2 * nu21 * nu32 * nu13);
  }

  if (Dim == 2) {
    Gamma = 1 / (1 - nu12 * nu21);
  }

  // Lamé's first parameters
  this->Cprime(0, 0) = E1 * (1 - nu23 * nu32) * Gamma;
  this->Cprime(1, 1) = E2 * (1 - nu13 * nu31) * Gamma;
  if (Dim == 3) {
    this->Cprime(2, 2) = E3 * (1 - nu12 * nu21) * Gamma;
  }

  // normalised poisson's ratio's
  this->Cprime(1, 0) = this->Cprime(0, 1) = E1 * (nu21 + nu31 * nu23) * Gamma;
  if (Dim == 3) {
    this->Cprime(2, 0) = this->Cprime(0, 2) = E1 * (nu31 + nu21 * nu32) * Gamma;
    this->Cprime(2, 1) = this->Cprime(1, 2) = E2 * (nu32 + nu12 * nu31) * Gamma;
  }

  // Lamé's second parameters (shear moduli)
  if (Dim == 3) {
    this->Cprime(3, 3) = G23;
    this->Cprime(4, 4) = G13;
    this->Cprime(5, 5) = G12;
  } else {
    this->Cprime(2, 2) = G12;
  }

  /* 1) rotation of C into the global frame */
  this->rotateCprime();
  this->C.eig(this->eigC);
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialElasticOrthotropic<spatial_dimension>::
    computePotentialEnergyByElement(ElementType type, UInt index,
                                    Vector<Real> & epot_on_quad_points) {

  Array<Real>::matrix_iterator gradu_it =
      this->gradu(type).begin(spatial_dimension, spatial_dimension);
  Array<Real>::matrix_iterator gradu_end =
      this->gradu(type).begin(spatial_dimension, spatial_dimension);
  Array<Real>::matrix_iterator stress_it =
      this->stress(type).begin(spatial_dimension, spatial_dimension);

  UInt nb_quadrature_points = this->fem.getNbIntegrationPoints(type);

  gradu_it += index * nb_quadrature_points;
  gradu_end += (index + 1) * nb_quadrature_points;
  stress_it += index * nb_quadrature_points;

  Real * epot_quad = epot_on_quad_points.storage();

  Matrix<Real> grad_u(spatial_dimension, spatial_dimension);

  for (; gradu_it != gradu_end; ++gradu_it, ++stress_it, ++epot_quad) {
    grad_u.copy(*gradu_it);

    this->computePotentialEnergyOnQuad(grad_u, *stress_it, *epot_quad);
  }
}

/* -------------------------------------------------------------------------- */
INSTANTIATE_MATERIAL(elastic_orthotropic, MaterialElasticOrthotropic);

} // namespace akantu
