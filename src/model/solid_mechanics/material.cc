/**
 * @file   material.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Tue Jul 27 2010
 * @date last modification: Wed Feb 21 2018
 *
 * @brief  Implementation of the common part of the material class
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
#include "material.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
Material::Material(SolidMechanicsModel & model, const ID & id)
    : Memory(id, model.getMemoryID()), Parsable(ParserType::_material, id),
      is_init(false), fem(model.getFEEngine()), finite_deformation(false),
      name(""), model(model),
      spatial_dimension(this->model.getSpatialDimension()),
      element_filter("element_filter", id, this->memory_id),
      stress("stress", *this), eigengradu("eigen_grad_u", *this),
      gradu("grad_u", *this), green_strain("green_strain", *this),
      piola_kirchhoff_2("piola_kirchhoff_2", *this),
      potential_energy("potential_energy", *this), is_non_local(false),
      use_previous_stress(false), use_previous_gradu(false),
      interpolation_inverse_coordinates("interpolation inverse coordinates",
                                        *this),
      interpolation_points_matrices("interpolation points matrices", *this),
      eigen_grad_u(model.getSpatialDimension(), model.getSpatialDimension(),
                   0.) {
  AKANTU_DEBUG_IN();

  this->registerParam("eigen_grad_u", eigen_grad_u, _pat_parsable,
                      "EigenGradU");

  /// for each connectivity types allocate the element filer array of the
  /// material
  element_filter.initialize(model.getMesh(),
                            _spatial_dimension = spatial_dimension,
                            _element_kind = _ek_regular);
  this->initialize();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Material::Material(SolidMechanicsModel & model, UInt dim, const Mesh & mesh,
                   FEEngine & fe_engine, const ID & id)
    : Memory(id, model.getMemoryID()), Parsable(ParserType::_material, id),
      is_init(false), fem(fe_engine), finite_deformation(false), name(""),
      model(model), spatial_dimension(dim),
      element_filter("element_filter", id, this->memory_id),
      stress("stress", *this, dim, fe_engine, this->element_filter),
      eigengradu("eigen_grad_u", *this, dim, fe_engine, this->element_filter),
      gradu("gradu", *this, dim, fe_engine, this->element_filter),
      green_strain("green_strain", *this, dim, fe_engine, this->element_filter),
      piola_kirchhoff_2("piola_kirchhoff_2", *this, dim, fe_engine,
                        this->element_filter),
      potential_energy("potential_energy", *this, dim, fe_engine,
                       this->element_filter),
      is_non_local(false), use_previous_stress(false),
      use_previous_gradu(false),
      interpolation_inverse_coordinates("interpolation inverse_coordinates",
                                        *this, dim, fe_engine,
                                        this->element_filter),
      interpolation_points_matrices("interpolation points matrices", *this, dim,
                                    fe_engine, this->element_filter) {

  AKANTU_DEBUG_IN();

  element_filter.initialize(mesh, _spatial_dimension = spatial_dimension,
                            _element_kind = _ek_regular);

  this->initialize();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Material::~Material() = default;

/* -------------------------------------------------------------------------- */
void Material::initialize() {
  registerParam("rho", rho, Real(0.), _pat_parsable | _pat_modifiable,
                "Density");
  registerParam("name", name, std::string(), _pat_parsable | _pat_readable);
  registerParam("finite_deformation", finite_deformation, false,
                _pat_parsable | _pat_readable, "Is finite deformation");
  registerParam("inelastic_deformation", inelastic_deformation, false,
                _pat_internal, "Is inelastic deformation");

  /// allocate gradu stress for local elements
  eigengradu.initialize(spatial_dimension * spatial_dimension);
  gradu.initialize(spatial_dimension * spatial_dimension);
  stress.initialize(spatial_dimension * spatial_dimension);

  potential_energy.initialize(1);

  this->model.registerEventHandler(*this);
}

/* -------------------------------------------------------------------------- */
void Material::initMaterial() {
  AKANTU_DEBUG_IN();

  if (finite_deformation) {
    this->piola_kirchhoff_2.initialize(spatial_dimension * spatial_dimension);
    if (use_previous_stress)
      this->piola_kirchhoff_2.initializeHistory();
    this->green_strain.initialize(spatial_dimension * spatial_dimension);
  }

  if (use_previous_stress)
    this->stress.initializeHistory();
  if (use_previous_gradu)
    this->gradu.initializeHistory();

  this->resizeInternals();

  auto dim = model.getSpatialDimension();
  for(const auto & type : element_filter.elementTypes()) {
    for(auto eigen_gradu : make_view(eigengradu(type), dim, dim)) {
      eigen_gradu = eigen_grad_u;
    }
  }

  is_init = true;

  updateInternalParameters();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::savePreviousState() {
  AKANTU_DEBUG_IN();

  for (auto pair : internal_vectors_real)
    if (pair.second->hasHistory())
      pair.second->saveCurrentValues();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::restorePreviousState() {
  AKANTU_DEBUG_IN();

  for (auto pair : internal_vectors_real)
    if (pair.second->hasHistory())
      pair.second->restorePreviousValues();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Compute the internal forces by assembling @f$\int_{e} \sigma_e \frac{\partial
 * \varphi}{\partial X} dX @f$
 *
 * @param[in] ghost_type compute the internal forces for _ghost or _not_ghost
 * element
 */
void Material::assembleInternalForces(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = model.getSpatialDimension();

  if (!finite_deformation) {

    auto & internal_force = const_cast<Array<Real> &>(model.getInternalForce());

    // Mesh & mesh = fem.getMesh();
    for (auto && type :
         element_filter.elementTypes(spatial_dimension, ghost_type)) {
      Array<UInt> & elem_filter = element_filter(type, ghost_type);
      UInt nb_element = elem_filter.size();

      if (nb_element == 0)
        continue;

      const Array<Real> & shapes_derivatives =
          fem.getShapesDerivatives(type, ghost_type);

      UInt size_of_shapes_derivatives = shapes_derivatives.getNbComponent();
      UInt nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);
      UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

      /// compute @f$\sigma \frac{\partial \varphi}{\partial X}@f$ by
      /// @f$\mathbf{B}^t \mathbf{\sigma}_q@f$
      Array<Real> * sigma_dphi_dx =
          new Array<Real>(nb_element * nb_quadrature_points,
                          size_of_shapes_derivatives, "sigma_x_dphi_/_dX");

      fem.computeBtD(stress(type, ghost_type), *sigma_dphi_dx, type, ghost_type,
                     elem_filter);

      /**
       * compute @f$\int \sigma  * \frac{\partial \varphi}{\partial X}dX@f$ by
       * @f$ \sum_q \mathbf{B}^t
       * \mathbf{\sigma}_q \overline w_q J_q@f$
       */
      Array<Real> * int_sigma_dphi_dx =
          new Array<Real>(nb_element, nb_nodes_per_element * spatial_dimension,
                          "int_sigma_x_dphi_/_dX");

      fem.integrate(*sigma_dphi_dx, *int_sigma_dphi_dx,
                    size_of_shapes_derivatives, type, ghost_type, elem_filter);
      delete sigma_dphi_dx;

      /// assemble
      model.getDOFManager().assembleElementalArrayLocalArray(
          *int_sigma_dphi_dx, internal_force, type, ghost_type, -1,
          elem_filter);
      delete int_sigma_dphi_dx;
    }
  } else {
    switch (spatial_dimension) {
    case 1:
      this->assembleInternalForces<1>(ghost_type);
      break;
    case 2:
      this->assembleInternalForces<2>(ghost_type);
      break;
    case 3:
      this->assembleInternalForces<3>(ghost_type);
      break;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Compute  the  stress from the gradu
 *
 * @param[in] ghost_type compute the residual for _ghost or _not_ghost element
 */
void Material::computeAllStresses(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = model.getSpatialDimension();

  for (const auto & type :
       element_filter.elementTypes(spatial_dimension, ghost_type)) {
    Array<UInt> & elem_filter = element_filter(type, ghost_type);

    if (elem_filter.size() == 0)
      continue;
    Array<Real> & gradu_vect = gradu(type, ghost_type);

    /// compute @f$\nabla u@f$
    fem.gradientOnIntegrationPoints(model.getDisplacement(), gradu_vect,
                                    spatial_dimension, type, ghost_type,
                                    elem_filter);

    gradu_vect -= eigengradu(type, ghost_type);

    /// compute @f$\mathbf{\sigma}_q@f$ from @f$\nabla u@f$
    computeStress(type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::computeAllCauchyStresses(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(finite_deformation, "The Cauchy stress can only be "
                                          "computed if you are working in "
                                          "finite deformation.");

  for (auto type : element_filter.elementTypes(spatial_dimension, ghost_type)) {
    switch (spatial_dimension) {
    case 1:
      this->StoCauchy<1>(type, ghost_type);
      break;
    case 2:
      this->StoCauchy<2>(type, ghost_type);
      break;
    case 3:
      this->StoCauchy<3>(type, ghost_type);
      break;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void Material::StoCauchy(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto gradu_it = this->gradu(el_type, ghost_type).begin(dim, dim);

  auto gradu_end = this->gradu(el_type, ghost_type).end(dim, dim);

  auto piola_it = this->piola_kirchhoff_2(el_type, ghost_type).begin(dim, dim);

  auto stress_it = this->stress(el_type, ghost_type).begin(dim, dim);

  for (; gradu_it != gradu_end; ++gradu_it, ++piola_it, ++stress_it) {
    Matrix<Real> & grad_u = *gradu_it;
    Matrix<Real> & piola = *piola_it;
    Matrix<Real> & sigma = *stress_it;

    auto F_tensor = gradUToF<dim>(grad_u);
    this->StoCauchy<dim>(F_tensor, piola, sigma);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::setToSteadyState(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  const Array<Real> & displacement = model.getDisplacement();

  // resizeInternalArray(gradu);

  UInt spatial_dimension = model.getSpatialDimension();

  for (auto type : element_filter.elementTypes(spatial_dimension, ghost_type)) {
    Array<UInt> & elem_filter = element_filter(type, ghost_type);
    Array<Real> & gradu_vect = gradu(type, ghost_type);

    /// compute @f$\nabla u@f$
    fem.gradientOnIntegrationPoints(displacement, gradu_vect, spatial_dimension,
                                    type, ghost_type, elem_filter);

    setToSteadyState(type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Compute  the stiffness  matrix by  assembling @f$\int_{\omega}  B^t  \times D
 * \times B d\omega @f$
 *
 * @param[in] ghost_type compute the residual for _ghost or _not_ghost element
 */
void Material::assembleStiffnessMatrix(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = model.getSpatialDimension();

  for (auto type : element_filter.elementTypes(spatial_dimension, ghost_type)) {
    if (finite_deformation) {
      switch (spatial_dimension) {
      case 1: {
        assembleStiffnessMatrixNL<1>(type, ghost_type);
        assembleStiffnessMatrixL2<1>(type, ghost_type);
        break;
      }
      case 2: {
        assembleStiffnessMatrixNL<2>(type, ghost_type);
        assembleStiffnessMatrixL2<2>(type, ghost_type);
        break;
      }
      case 3: {
        assembleStiffnessMatrixNL<3>(type, ghost_type);
        assembleStiffnessMatrixL2<3>(type, ghost_type);
        break;
      }
      }
    } else {
      switch (spatial_dimension) {
      case 1: {
        assembleStiffnessMatrix<1>(type, ghost_type);
        break;
      }
      case 2: {
        assembleStiffnessMatrix<2>(type, ghost_type);
        break;
      }
      case 3: {
        assembleStiffnessMatrix<3>(type, ghost_type);
        break;
      }
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void Material::assembleStiffnessMatrix(const ElementType & type,
                                       GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Array<UInt> & elem_filter = element_filter(type, ghost_type);
  if (elem_filter.size() == 0) {
    AKANTU_DEBUG_OUT();
    return;
  }

  // const Array<Real> & shapes_derivatives =
  //     fem.getShapesDerivatives(type, ghost_type);

  Array<Real> & gradu_vect = gradu(type, ghost_type);

  UInt nb_element = elem_filter.size();
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);

  gradu_vect.resize(nb_quadrature_points * nb_element);

  fem.gradientOnIntegrationPoints(model.getDisplacement(), gradu_vect, dim,
                                  type, ghost_type, elem_filter);

  UInt tangent_size = getTangentStiffnessVoigtSize(dim);

  Array<Real> * tangent_stiffness_matrix =
      new Array<Real>(nb_element * nb_quadrature_points,
                      tangent_size * tangent_size, "tangent_stiffness_matrix");

  tangent_stiffness_matrix->clear();

  computeTangentModuli(type, *tangent_stiffness_matrix, ghost_type);

  /// compute @f$\mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  UInt bt_d_b_size = dim * nb_nodes_per_element;

  Array<Real> * bt_d_b = new Array<Real>(nb_element * nb_quadrature_points,
                                         bt_d_b_size * bt_d_b_size, "B^t*D*B");

  fem.computeBtDB(*tangent_stiffness_matrix, *bt_d_b, 4, type, ghost_type,
                  elem_filter);

  delete tangent_stiffness_matrix;

  /// compute @f$ k_e = \int_e \mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  Array<Real> * K_e =
      new Array<Real>(nb_element, bt_d_b_size * bt_d_b_size, "K_e");

  fem.integrate(*bt_d_b, *K_e, bt_d_b_size * bt_d_b_size, type, ghost_type,
                elem_filter);

  delete bt_d_b;

  model.getDOFManager().assembleElementalMatricesToMatrix(
      "K", "displacement", *K_e, type, ghost_type, _symmetric, elem_filter);
  delete K_e;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void Material::assembleStiffnessMatrixNL(const ElementType & type,
                                         GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  const Array<Real> & shapes_derivatives =
      fem.getShapesDerivatives(type, ghost_type);

  Array<UInt> & elem_filter = element_filter(type, ghost_type);
  // Array<Real> & gradu_vect = delta_gradu(type, ghost_type);

  UInt nb_element = elem_filter.size();
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);

  Array<Real> * shapes_derivatives_filtered = new Array<Real>(
      nb_element * nb_quadrature_points, dim * nb_nodes_per_element,
      "shapes derivatives filtered");

  fem.filterElementalData(fem.getMesh(), shapes_derivatives,
                          *shapes_derivatives_filtered, type, ghost_type,
                          elem_filter);

  /// compute @f$\mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  UInt bt_s_b_size = dim * nb_nodes_per_element;

  Array<Real> * bt_s_b = new Array<Real>(nb_element * nb_quadrature_points,
                                         bt_s_b_size * bt_s_b_size, "B^t*D*B");

  UInt piola_matrix_size = getCauchyStressMatrixSize(dim);

  Matrix<Real> B(piola_matrix_size, bt_s_b_size);
  Matrix<Real> Bt_S(bt_s_b_size, piola_matrix_size);
  Matrix<Real> S(piola_matrix_size, piola_matrix_size);

  auto shapes_derivatives_filtered_it = shapes_derivatives_filtered->begin(
      spatial_dimension, nb_nodes_per_element);

  auto Bt_S_B_it = bt_s_b->begin(bt_s_b_size, bt_s_b_size);
  auto Bt_S_B_end = bt_s_b->end(bt_s_b_size, bt_s_b_size);
  auto piola_it = piola_kirchhoff_2(type, ghost_type).begin(dim, dim);

  for (; Bt_S_B_it != Bt_S_B_end;
       ++Bt_S_B_it, ++shapes_derivatives_filtered_it, ++piola_it) {
    auto & Bt_S_B = *Bt_S_B_it;
    const auto & Piola_kirchhoff_matrix = *piola_it;

    setCauchyStressMatrix<dim>(Piola_kirchhoff_matrix, S);
    VoigtHelper<dim>::transferBMatrixToBNL(*shapes_derivatives_filtered_it, B,
                                           nb_nodes_per_element);
    Bt_S.template mul<true, false>(B, S);
    Bt_S_B.template mul<false, false>(Bt_S, B);
  }

  delete shapes_derivatives_filtered;

  /// compute @f$ k_e = \int_e \mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  Array<Real> * K_e =
      new Array<Real>(nb_element, bt_s_b_size * bt_s_b_size, "K_e");

  fem.integrate(*bt_s_b, *K_e, bt_s_b_size * bt_s_b_size, type, ghost_type,
                elem_filter);

  delete bt_s_b;

  model.getDOFManager().assembleElementalMatricesToMatrix(
      "K", "displacement", *K_e, type, ghost_type, _symmetric, elem_filter);

  delete K_e;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void Material::assembleStiffnessMatrixL2(const ElementType & type,
                                         GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  const Array<Real> & shapes_derivatives =
      fem.getShapesDerivatives(type, ghost_type);

  Array<UInt> & elem_filter = element_filter(type, ghost_type);
  Array<Real> & gradu_vect = gradu(type, ghost_type);

  UInt nb_element = elem_filter.size();
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);

  gradu_vect.resize(nb_quadrature_points * nb_element);

  fem.gradientOnIntegrationPoints(model.getDisplacement(), gradu_vect, dim,
                                  type, ghost_type, elem_filter);

  UInt tangent_size = getTangentStiffnessVoigtSize(dim);

  Array<Real> * tangent_stiffness_matrix =
      new Array<Real>(nb_element * nb_quadrature_points,
                      tangent_size * tangent_size, "tangent_stiffness_matrix");

  tangent_stiffness_matrix->clear();

  computeTangentModuli(type, *tangent_stiffness_matrix, ghost_type);

  Array<Real> * shapes_derivatives_filtered = new Array<Real>(
      nb_element * nb_quadrature_points, dim * nb_nodes_per_element,
      "shapes derivatives filtered");
  fem.filterElementalData(fem.getMesh(), shapes_derivatives,
                          *shapes_derivatives_filtered, type, ghost_type,
                          elem_filter);

  /// compute @f$\mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  UInt bt_d_b_size = dim * nb_nodes_per_element;

  Array<Real> * bt_d_b = new Array<Real>(nb_element * nb_quadrature_points,
                                         bt_d_b_size * bt_d_b_size, "B^t*D*B");

  Matrix<Real> B(tangent_size, dim * nb_nodes_per_element);
  Matrix<Real> B2(tangent_size, dim * nb_nodes_per_element);
  Matrix<Real> Bt_D(dim * nb_nodes_per_element, tangent_size);

  auto shapes_derivatives_filtered_it = shapes_derivatives_filtered->begin(
      spatial_dimension, nb_nodes_per_element);

  auto Bt_D_B_it = bt_d_b->begin(bt_d_b_size, bt_d_b_size);
  auto grad_u_it = gradu_vect.begin(dim, dim);
  auto D_it = tangent_stiffness_matrix->begin(tangent_size, tangent_size);
  auto D_end = tangent_stiffness_matrix->end(tangent_size, tangent_size);

  for (; D_it != D_end;
       ++D_it, ++Bt_D_B_it, ++shapes_derivatives_filtered_it, ++grad_u_it) {
    const auto & grad_u = *grad_u_it;
    const auto & D = *D_it;
    auto & Bt_D_B = *Bt_D_B_it;

    // transferBMatrixToBL1<dim > (*shapes_derivatives_filtered_it, B,
    // nb_nodes_per_element);
    VoigtHelper<dim>::transferBMatrixToSymVoigtBMatrix(
        *shapes_derivatives_filtered_it, B, nb_nodes_per_element);
    VoigtHelper<dim>::transferBMatrixToBL2(*shapes_derivatives_filtered_it,
                                           grad_u, B2, nb_nodes_per_element);
    B += B2;
    Bt_D.template mul<true, false>(B, D);
    Bt_D_B.template mul<false, false>(Bt_D, B);
  }

  delete tangent_stiffness_matrix;
  delete shapes_derivatives_filtered;

  /// compute @f$ k_e = \int_e \mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  Array<Real> * K_e =
      new Array<Real>(nb_element, bt_d_b_size * bt_d_b_size, "K_e");

  fem.integrate(*bt_d_b, *K_e, bt_d_b_size * bt_d_b_size, type, ghost_type,
                elem_filter);

  delete bt_d_b;

  model.getDOFManager().assembleElementalMatricesToMatrix(
      "K", "displacement", *K_e, type, ghost_type, _symmetric, elem_filter);
  delete K_e;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void Material::assembleInternalForces(GhostType ghost_type) {

  AKANTU_DEBUG_IN();

  Array<Real> & internal_force = model.getInternalForce();

  Mesh & mesh = fem.getMesh();
  for (auto type : element_filter.elementTypes(_ghost_type = ghost_type)) {
    const Array<Real> & shapes_derivatives =
        fem.getShapesDerivatives(type, ghost_type);

    Array<UInt> & elem_filter = element_filter(type, ghost_type);
    if (elem_filter.size() == 0)
      continue;
    UInt size_of_shapes_derivatives = shapes_derivatives.getNbComponent();
    UInt nb_element = elem_filter.size();
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    UInt nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);

    Array<Real> * shapesd_filtered = new Array<Real>(
        nb_element, size_of_shapes_derivatives, "filtered shapesd");

    fem.filterElementalData(mesh, shapes_derivatives, *shapesd_filtered, type,
                            ghost_type, elem_filter);

    Array<Real>::matrix_iterator shapes_derivatives_filtered_it =
        shapesd_filtered->begin(dim, nb_nodes_per_element);

    // Set stress vectors
    UInt stress_size = getTangentStiffnessVoigtSize(dim);

    // Set matrices B and BNL*
    UInt bt_s_size = dim * nb_nodes_per_element;

    auto * bt_s =
        new Array<Real>(nb_element * nb_quadrature_points, bt_s_size, "B^t*S");

    auto grad_u_it = this->gradu(type, ghost_type).begin(dim, dim);
    auto grad_u_end = this->gradu(type, ghost_type).end(dim, dim);
    auto stress_it = this->piola_kirchhoff_2(type, ghost_type).begin(dim, dim);
    shapes_derivatives_filtered_it =
        shapesd_filtered->begin(dim, nb_nodes_per_element);

    Array<Real>::matrix_iterator bt_s_it = bt_s->begin(bt_s_size, 1);

    Matrix<Real> B_tensor(stress_size, bt_s_size);
    Matrix<Real> B2_tensor(stress_size, bt_s_size);

    for (; grad_u_it != grad_u_end; ++grad_u_it, ++stress_it,
                                    ++shapes_derivatives_filtered_it,
                                    ++bt_s_it) {
      auto & grad_u = *grad_u_it;
      auto & r = *bt_s_it;
      auto & S = *stress_it;

      VoigtHelper<dim>::transferBMatrixToSymVoigtBMatrix(
          *shapes_derivatives_filtered_it, B_tensor, nb_nodes_per_element);

      VoigtHelper<dim>::transferBMatrixToBL2(*shapes_derivatives_filtered_it,
                                             grad_u, B2_tensor,
                                             nb_nodes_per_element);

      B_tensor += B2_tensor;

      auto S_vect = Material::stressToVoigt<dim>(S);
      Matrix<Real> S_voigt(S_vect.storage(), stress_size, 1);

      r.template mul<true, false>(B_tensor, S_voigt);
    }

    delete shapesd_filtered;

    /// compute @f$ k_e = \int_e \mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
    Array<Real> * r_e = new Array<Real>(nb_element, bt_s_size, "r_e");

    fem.integrate(*bt_s, *r_e, bt_s_size, type, ghost_type, elem_filter);

    delete bt_s;

    model.getDOFManager().assembleElementalArrayLocalArray(
        *r_e, internal_force, type, ghost_type, -1., elem_filter);
    delete r_e;
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::computePotentialEnergyByElements() {
  AKANTU_DEBUG_IN();

  for (auto type : element_filter.elementTypes(spatial_dimension, _not_ghost)) {
    computePotentialEnergy(type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::computePotentialEnergy(ElementType) {
  AKANTU_DEBUG_IN();
  AKANTU_TO_IMPLEMENT();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Real Material::getPotentialEnergy() {
  AKANTU_DEBUG_IN();
  Real epot = 0.;

  computePotentialEnergyByElements();

  /// integrate the potential energy for each type of elements
  for (auto type : element_filter.elementTypes(spatial_dimension, _not_ghost)) {
    epot += fem.integrate(potential_energy(type, _not_ghost), type, _not_ghost,
                          element_filter(type, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return epot;
}

/* -------------------------------------------------------------------------- */
Real Material::getPotentialEnergy(ElementType & type, UInt index) {
  AKANTU_DEBUG_IN();
  Real epot = 0.;

  Vector<Real> epot_on_quad_points(fem.getNbIntegrationPoints(type));

  computePotentialEnergyByElement(type, index, epot_on_quad_points);

  epot = fem.integrate(epot_on_quad_points, type, element_filter(type)(index));

  AKANTU_DEBUG_OUT();
  return epot;
}

/* -------------------------------------------------------------------------- */
Real Material::getEnergy(const std::string & type) {
  AKANTU_DEBUG_IN();
  if (type == "potential")
    return getPotentialEnergy();
  AKANTU_DEBUG_OUT();
  return 0.;
}

/* -------------------------------------------------------------------------- */
Real Material::getEnergy(const std::string & energy_id, ElementType type,
                         UInt index) {
  AKANTU_DEBUG_IN();
  if (energy_id == "potential")
    return getPotentialEnergy(type, index);
  AKANTU_DEBUG_OUT();
  return 0.;
}

/* -------------------------------------------------------------------------- */
void Material::initElementalFieldInterpolation(
    const ElementTypeMapArray<Real> & interpolation_points_coordinates) {
  AKANTU_DEBUG_IN();

  this->fem.initElementalFieldInterpolationFromIntegrationPoints(
      interpolation_points_coordinates, this->interpolation_points_matrices,
      this->interpolation_inverse_coordinates, &(this->element_filter));

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::interpolateStress(ElementTypeMapArray<Real> & result,
                                 const GhostType ghost_type) {

  this->fem.interpolateElementalFieldFromIntegrationPoints(
      this->stress, this->interpolation_points_matrices,
      this->interpolation_inverse_coordinates, result, ghost_type,
      &(this->element_filter));
}

/* -------------------------------------------------------------------------- */
void Material::interpolateStressOnFacets(
    ElementTypeMapArray<Real> & result,
    ElementTypeMapArray<Real> & by_elem_result, const GhostType ghost_type) {

  interpolateStress(by_elem_result, ghost_type);

  UInt stress_size = this->stress.getNbComponent();

  const Mesh & mesh = this->model.getMesh();
  const Mesh & mesh_facets = mesh.getMeshFacets();

  for (auto type : element_filter.elementTypes(spatial_dimension, ghost_type)) {
    Array<UInt> & elem_fil = element_filter(type, ghost_type);
    Array<Real> & by_elem_res = by_elem_result(type, ghost_type);
    UInt nb_element = elem_fil.size();
    UInt nb_element_full = this->model.getMesh().getNbElement(type, ghost_type);
    UInt nb_interpolation_points_per_elem =
        by_elem_res.size() / nb_element_full;

    const Array<Element> & facet_to_element =
        mesh_facets.getSubelementToElement(type, ghost_type);
    ElementType type_facet = Mesh::getFacetType(type);
    UInt nb_facet_per_elem = facet_to_element.getNbComponent();
    UInt nb_quad_per_facet =
        nb_interpolation_points_per_elem / nb_facet_per_elem;
    Element element_for_comparison{type, 0, ghost_type};
    const Array<std::vector<Element>> * element_to_facet = nullptr;
    GhostType current_ghost_type = _casper;
    Array<Real> * result_vec = nullptr;

    Array<Real>::const_matrix_iterator result_it =
        by_elem_res.begin_reinterpret(
            stress_size, nb_interpolation_points_per_elem, nb_element_full);

    for (UInt el = 0; el < nb_element; ++el) {
      UInt global_el = elem_fil(el);
      element_for_comparison.element = global_el;
      for (UInt f = 0; f < nb_facet_per_elem; ++f) {

        Element facet_elem = facet_to_element(global_el, f);
        UInt global_facet = facet_elem.element;
        if (facet_elem.ghost_type != current_ghost_type) {
          current_ghost_type = facet_elem.ghost_type;
          element_to_facet = &mesh_facets.getElementToSubelement(
              type_facet, current_ghost_type);
          result_vec = &result(type_facet, current_ghost_type);
        }

        bool is_second_element =
            (*element_to_facet)(global_facet)[0] != element_for_comparison;

        for (UInt q = 0; q < nb_quad_per_facet; ++q) {
          Vector<Real> result_local(result_vec->storage() +
                                        (global_facet * nb_quad_per_facet + q) *
                                            result_vec->getNbComponent() +
                                        is_second_element * stress_size,
                                    stress_size);

          const Matrix<Real> & result_tmp(result_it[global_el]);
          result_local = result_tmp(f * nb_quad_per_facet + q);
        }
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
template <typename T>
const Array<T> & Material::getArray(const ID & /*vect_id*/,
                                    const ElementType & /*type*/,
                                    const GhostType & /*ghost_type*/) const {
  AKANTU_TO_IMPLEMENT();
  return NULL;
}

/* -------------------------------------------------------------------------- */
template <typename T>
Array<T> & Material::getArray(const ID & /*vect_id*/,
                              const ElementType & /*type*/,
                              const GhostType & /*ghost_type*/) {
  AKANTU_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template <>
const Array<Real> & Material::getArray(const ID & vect_id,
                                       const ElementType & type,
                                       const GhostType & ghost_type) const {
  std::stringstream sstr;
  std::string ghost_id = "";
  if (ghost_type == _ghost)
    ghost_id = ":ghost";
  sstr << getID() << ":" << vect_id << ":" << type << ghost_id;

  ID fvect_id = sstr.str();
  try {
    return Memory::getArray<Real>(fvect_id);
  } catch (debug::Exception & e) {
    AKANTU_SILENT_EXCEPTION("The material " << name << "(" << getID()
                                            << ") does not contain a vector "
                                            << vect_id << " (" << fvect_id
                                            << ") [" << e << "]");
  }
}

/* -------------------------------------------------------------------------- */
template <>
Array<Real> & Material::getArray(const ID & vect_id, const ElementType & type,
                                 const GhostType & ghost_type) {
  std::stringstream sstr;
  std::string ghost_id = "";
  if (ghost_type == _ghost)
    ghost_id = ":ghost";
  sstr << getID() << ":" << vect_id << ":" << type << ghost_id;

  ID fvect_id = sstr.str();
  try {
    return Memory::getArray<Real>(fvect_id);
  } catch (debug::Exception & e) {
    AKANTU_SILENT_EXCEPTION("The material " << name << "(" << getID()
                                            << ") does not contain a vector "
                                            << vect_id << " (" << fvect_id
                                            << ") [" << e << "]");
  }
}

/* -------------------------------------------------------------------------- */
template <>
const Array<UInt> & Material::getArray(const ID & vect_id,
                                       const ElementType & type,
                                       const GhostType & ghost_type) const {
  std::stringstream sstr;
  std::string ghost_id = "";
  if (ghost_type == _ghost)
    ghost_id = ":ghost";
  sstr << getID() << ":" << vect_id << ":" << type << ghost_id;

  ID fvect_id = sstr.str();
  try {
    return Memory::getArray<UInt>(fvect_id);
  } catch (debug::Exception & e) {
    AKANTU_SILENT_EXCEPTION("The material " << name << "(" << getID()
                                            << ") does not contain a vector "
                                            << vect_id << " (" << fvect_id
                                            << ") [" << e << "]");
  }
}

/* -------------------------------------------------------------------------- */
template <>
Array<UInt> & Material::getArray(const ID & vect_id, const ElementType & type,
                                 const GhostType & ghost_type) {
  std::stringstream sstr;
  std::string ghost_id = "";
  if (ghost_type == _ghost)
    ghost_id = ":ghost";
  sstr << getID() << ":" << vect_id << ":" << type << ghost_id;

  ID fvect_id = sstr.str();
  try {
    return Memory::getArray<UInt>(fvect_id);
  } catch (debug::Exception & e) {
    AKANTU_SILENT_EXCEPTION("The material " << name << "(" << getID()
                                            << ") does not contain a vector "
                                            << vect_id << "(" << fvect_id
                                            << ") [" << e << "]");
  }
}

/* -------------------------------------------------------------------------- */
template <typename T>
const InternalField<T> &
Material::getInternal([[gnu::unused]] const ID & int_id) const {
  AKANTU_TO_IMPLEMENT();
  return NULL;
}

/* -------------------------------------------------------------------------- */
template <typename T>
InternalField<T> & Material::getInternal([[gnu::unused]] const ID & int_id) {
  AKANTU_TO_IMPLEMENT();
  return NULL;
}

/* -------------------------------------------------------------------------- */
template <>
const InternalField<Real> & Material::getInternal(const ID & int_id) const {
  auto it = internal_vectors_real.find(getID() + ":" + int_id);
  if (it == internal_vectors_real.end()) {
    AKANTU_SILENT_EXCEPTION("The material " << name << "(" << getID()
                                            << ") does not contain an internal "
                                            << int_id << " ("
                                            << (getID() + ":" + int_id) << ")");
  }
  return *it->second;
}

/* -------------------------------------------------------------------------- */
template <> InternalField<Real> & Material::getInternal(const ID & int_id) {
  auto it = internal_vectors_real.find(getID() + ":" + int_id);
  if (it == internal_vectors_real.end()) {
    AKANTU_SILENT_EXCEPTION("The material " << name << "(" << getID()
                                            << ") does not contain an internal "
                                            << int_id << " ("
                                            << (getID() + ":" + int_id) << ")");
  }
  return *it->second;
}

/* -------------------------------------------------------------------------- */
template <>
const InternalField<UInt> & Material::getInternal(const ID & int_id) const {
  auto it = internal_vectors_uint.find(getID() + ":" + int_id);
  if (it == internal_vectors_uint.end()) {
    AKANTU_SILENT_EXCEPTION("The material " << name << "(" << getID()
                                            << ") does not contain an internal "
                                            << int_id << " ("
                                            << (getID() + ":" + int_id) << ")");
  }
  return *it->second;
}

/* -------------------------------------------------------------------------- */
template <> InternalField<UInt> & Material::getInternal(const ID & int_id) {
  auto it = internal_vectors_uint.find(getID() + ":" + int_id);
  if (it == internal_vectors_uint.end()) {
    AKANTU_SILENT_EXCEPTION("The material " << name << "(" << getID()
                                            << ") does not contain an internal "
                                            << int_id << " ("
                                            << (getID() + ":" + int_id) << ")");
  }
  return *it->second;
}

/* -------------------------------------------------------------------------- */
void Material::addElements(const Array<Element> & elements_to_add) {
  AKANTU_DEBUG_IN();
  UInt mat_id = model.getInternalIndexFromID(getID());
  Array<Element>::const_iterator<Element> el_begin = elements_to_add.begin();
  Array<Element>::const_iterator<Element> el_end = elements_to_add.end();
  for (; el_begin != el_end; ++el_begin) {
    const Element & element = *el_begin;
    Array<UInt> & mat_indexes =
        model.getMaterialByElement(element.type, element.ghost_type);
    Array<UInt> & mat_loc_num =
        model.getMaterialLocalNumbering(element.type, element.ghost_type);

    UInt index =
        this->addElement(element.type, element.element, element.ghost_type);
    mat_indexes(element.element) = mat_id;
    mat_loc_num(element.element) = index;
  }

  this->resizeInternals();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::removeElements(const Array<Element> & elements_to_remove) {
  AKANTU_DEBUG_IN();

  Array<Element>::const_iterator<Element> el_begin = elements_to_remove.begin();
  Array<Element>::const_iterator<Element> el_end = elements_to_remove.end();

  if (el_begin == el_end)
    return;

  ElementTypeMapArray<UInt> material_local_new_numbering(
      "remove mat filter elem", getID(), getMemoryID());

  Element element;
  for (auto ghost_type : ghost_types) {
    element.ghost_type = ghost_type;

    for (auto & type : element_filter.elementTypes(_ghost_type = ghost_type,
                       _element_kind = _ek_not_defined)) {
      element.type = type;

      Array<UInt> & elem_filter = this->element_filter(type, ghost_type);
      Array<UInt> & mat_loc_num =
          this->model.getMaterialLocalNumbering(type, ghost_type);

      if (!material_local_new_numbering.exists(type, ghost_type))
        material_local_new_numbering.alloc(elem_filter.size(), 1, type,
                                           ghost_type);
      Array<UInt> & mat_renumbering =
          material_local_new_numbering(type, ghost_type);

      UInt nb_element = elem_filter.size();
      Array<UInt> elem_filter_tmp;

      UInt new_id = 0;
      for (UInt el = 0; el < nb_element; ++el) {
        element.element = elem_filter(el);

        if (std::find(el_begin, el_end, element) == el_end) {
          elem_filter_tmp.push_back(element.element);

          mat_renumbering(el) = new_id;
          mat_loc_num(element.element) = new_id;
          ++new_id;
        } else {
          mat_renumbering(el) = UInt(-1);
        }
      }

      elem_filter.resize(elem_filter_tmp.size());
      elem_filter.copy(elem_filter_tmp);
    }
  }

  for (auto it = internal_vectors_real.begin();
       it != internal_vectors_real.end(); ++it)
    it->second->removeIntegrationPoints(material_local_new_numbering);

  for (auto it = internal_vectors_uint.begin();
       it != internal_vectors_uint.end(); ++it)
    it->second->removeIntegrationPoints(material_local_new_numbering);

  for (auto it = internal_vectors_bool.begin();
       it != internal_vectors_bool.end(); ++it)
    it->second->removeIntegrationPoints(material_local_new_numbering);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::resizeInternals() {
  AKANTU_DEBUG_IN();
  for (auto it = internal_vectors_real.begin();
       it != internal_vectors_real.end(); ++it)
    it->second->resize();

  for (auto it = internal_vectors_uint.begin();
       it != internal_vectors_uint.end(); ++it)
    it->second->resize();

  for (auto it = internal_vectors_bool.begin();
       it != internal_vectors_bool.end(); ++it)
    it->second->resize();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::onElementsAdded(const Array<Element> &,
                               const NewElementsEvent &) {
  this->resizeInternals();
}

/* -------------------------------------------------------------------------- */
void Material::onElementsRemoved(
    const Array<Element> & element_list,
    const ElementTypeMapArray<UInt> & new_numbering,
    [[gnu::unused]] const RemovedElementsEvent & event) {
  UInt my_num = model.getInternalIndexFromID(getID());

  ElementTypeMapArray<UInt> material_local_new_numbering(
      "remove mat filter elem", getID(), getMemoryID());

  auto el_begin = element_list.begin();
  auto el_end = element_list.end();

  for (auto && gt : ghost_types) {
    for (auto && type :
         new_numbering.elementTypes(_all_dimensions, gt, _ek_not_defined)) {

      if (not element_filter.exists(type, gt) ||
          element_filter(type, gt).size() == 0)
        continue;

      auto & elem_filter = element_filter(type, gt);
      auto & mat_indexes = this->model.getMaterialByElement(type, gt);
      auto & mat_loc_num = this->model.getMaterialLocalNumbering(type, gt);
      auto nb_element = this->model.getMesh().getNbElement(type, gt);

      // all materials will resize of the same size...
      mat_indexes.resize(nb_element);
      mat_loc_num.resize(nb_element);

      if (!material_local_new_numbering.exists(type, gt))
        material_local_new_numbering.alloc(elem_filter.size(), 1, type, gt);

      auto & mat_renumbering = material_local_new_numbering(type, gt);
      const auto & renumbering = new_numbering(type, gt);
      Array<UInt> elem_filter_tmp;
      UInt ni = 0;
      Element el{type, 0, gt};

      for (UInt i = 0; i < elem_filter.size(); ++i) {
        el.element = elem_filter(i);
        if (std::find(el_begin, el_end, el) == el_end) {
          UInt new_el = renumbering(el.element);
          AKANTU_DEBUG_ASSERT(new_el != UInt(-1),
                              "A not removed element as been badly renumbered");
          elem_filter_tmp.push_back(new_el);
          mat_renumbering(i) = ni;

          mat_indexes(new_el) = my_num;
          mat_loc_num(new_el) = ni;
          ++ni;
        } else {
          mat_renumbering(i) = UInt(-1);
        }
      }

      elem_filter.resize(elem_filter_tmp.size());
      elem_filter.copy(elem_filter_tmp);
    }
  }

  for (auto it = internal_vectors_real.begin();
       it != internal_vectors_real.end(); ++it)
    it->second->removeIntegrationPoints(material_local_new_numbering);

  for (auto it = internal_vectors_uint.begin();
       it != internal_vectors_uint.end(); ++it)
    it->second->removeIntegrationPoints(material_local_new_numbering);

  for (auto it = internal_vectors_bool.begin();
       it != internal_vectors_bool.end(); ++it)
    it->second->removeIntegrationPoints(material_local_new_numbering);
}

/* -------------------------------------------------------------------------- */
void Material::beforeSolveStep() { this->savePreviousState(); }

/* -------------------------------------------------------------------------- */
void Material::afterSolveStep(bool converged) {
  if (not converged) {
    this->restorePreviousState();
    return;
  }

  for (auto & type : element_filter.elementTypes(_all_dimensions, _not_ghost,
                                                 _ek_not_defined)) {
    this->updateEnergies(type);
  }
}
/* -------------------------------------------------------------------------- */
void Material::onDamageIteration() { this->savePreviousState(); }

/* -------------------------------------------------------------------------- */
void Material::onDamageUpdate() {
  for (auto & type : element_filter.elementTypes(_all_dimensions, _not_ghost,
                                                 _ek_not_defined)) {
    this->updateEnergiesAfterDamage(type);
  }
}

/* -------------------------------------------------------------------------- */
void Material::onDump() {
  if (this->isFiniteDeformation())
    this->computeAllCauchyStresses(_not_ghost);
}

/* -------------------------------------------------------------------------- */
void Material::printself(std::ostream & stream, int indent) const {
  std::string space(indent, AKANTU_INDENT);
  std::string type = getID().substr(getID().find_last_of(':') + 1);

  stream << space << "Material " << type << " [" << std::endl;
  Parsable::printself(stream, indent);
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
/// extrapolate internal values
void Material::extrapolateInternal(const ID & id, const Element & element,
                                   [[gnu::unused]] const Matrix<Real> & point,
                                   Matrix<Real> & extrapolated) {
  if (this->isInternal<Real>(id, element.kind())) {
    UInt nb_element =
        this->element_filter(element.type, element.ghost_type).size();
    const ID name = this->getID() + ":" + id;
    UInt nb_quads =
        this->internal_vectors_real[name]->getFEEngine().getNbIntegrationPoints(
            element.type, element.ghost_type);
    const Array<Real> & internal =
        this->getArray<Real>(id, element.type, element.ghost_type);
    UInt nb_component = internal.getNbComponent();
    Array<Real>::const_matrix_iterator internal_it =
        internal.begin_reinterpret(nb_component, nb_quads, nb_element);
    Element local_element = this->convertToLocalElement(element);

    /// instead of really extrapolating, here the value of the first GP
    /// is copied into the result vector. This works only for linear
    /// elements
    /// @todo extrapolate!!!!
    AKANTU_DEBUG_WARNING("This is a fix, values are not truly extrapolated");

    const Matrix<Real> & values = internal_it[local_element.element];
    UInt index = 0;
    Vector<Real> tmp(nb_component);
    for (UInt j = 0; j < values.cols(); ++j) {
      tmp = values(j);
      if (tmp.norm() > 0) {
        index = j;
        break;
      }
    }

    for (UInt i = 0; i < extrapolated.size(); ++i) {
      extrapolated(i) = values(index);
    }
  } else {
    Matrix<Real> default_values(extrapolated.rows(), extrapolated.cols(), 0.);
    extrapolated = default_values;
  }
}

/* -------------------------------------------------------------------------- */
void Material::applyEigenGradU(const Matrix<Real> & prescribed_eigen_grad_u,
                               const GhostType ghost_type) {

  for (auto && type : element_filter.elementTypes(_all_dimensions, _not_ghost,
                                                  _ek_not_defined)) {
    if (!element_filter(type, ghost_type).size())
      continue;
    auto eigen_it = this->eigengradu(type, ghost_type)
                        .begin(spatial_dimension, spatial_dimension);
    auto eigen_end = this->eigengradu(type, ghost_type)
                         .end(spatial_dimension, spatial_dimension);
    for (; eigen_it != eigen_end; ++eigen_it) {
      auto & current_eigengradu = *eigen_it;
      current_eigengradu = prescribed_eigen_grad_u;
    }
  }
}

/* -------------------------------------------------------------------------- */
MaterialFactory & Material::getFactory() {
  return MaterialFactory::getInstance();
}

} // namespace akantu
