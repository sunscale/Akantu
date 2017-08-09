/**
 * @file   material.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Tue Jul 27 2010
 * @date last modification: Tue Nov 24 2015
 *
 * @brief  Implementation of the common part of the material class
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "material.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
Material::Material(SolidMechanicsModel & model, const ID & id)
    : Memory(id, model.getMemoryID()), Parsable(_st_material, id),
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
      interpolation_points_matrices("interpolation points matrices", *this) {
  AKANTU_DEBUG_IN();

  /// for each connectivity types allocate the element filer array of the
  /// material
  element_filter.initialize(model.getMesh(),
                            _spatial_dimension = spatial_dimension);
  // model.getMesh().initElementTypeMapArray(element_filter, 1,
  // spatial_dimension,
  //                                         false, _ek_regular);

  this->initialize();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Material::Material(SolidMechanicsModel & model, UInt dim, const Mesh & mesh,
                   FEEngine & fe_engine, const ID & id)
    : Memory(id, model.getMemoryID()), Parsable(_st_material, id),
      is_init(false), fem(model.getFEEngine()), finite_deformation(false),
      name(""), model(model), spatial_dimension(dim),
      element_filter("element_filter", id, this->memory_id),
      stress("stress", *this, dim, fe_engine, this->element_filter),
      eigengradu("eigen_grad_u", *this, dim, fe_engine, this->element_filter),
      gradu("gradu", *this, dim, fe_engine, this->element_filter),
      green_strain("green_strain", *this, dim, fe_engine, this->element_filter),
      piola_kirchhoff_2("poila_kirchhoff_2", *this, dim, fe_engine,
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
  element_filter.initialize(mesh, _spatial_dimension = spatial_dimension);
  // mesh.initElementTypeMapArray(element_filter, 1, spatial_dimension, false,
  //                              _ek_regular);

  this->initialize();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Material::~Material() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

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

  for (std::map<ID, InternalField<Real> *>::iterator it =
           internal_vectors_real.begin();
       it != internal_vectors_real.end(); ++it)
    it->second->resize();

  for (std::map<ID, InternalField<UInt> *>::iterator it =
           internal_vectors_uint.begin();
       it != internal_vectors_uint.end(); ++it)
    it->second->resize();

  for (std::map<ID, InternalField<bool> *>::iterator it =
           internal_vectors_bool.begin();
       it != internal_vectors_bool.end(); ++it)
    it->second->resize();

  is_init = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::savePreviousState() {
  AKANTU_DEBUG_IN();

  for (std::map<ID, InternalField<Real> *>::iterator it =
           internal_vectors_real.begin();
       it != internal_vectors_real.end(); ++it) {
    if (it->second->hasHistory())
      it->second->saveCurrentValues();
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Compute  the  residual  by  assembling  @f$\int_{e}  \sigma_e  \frac{\partial
 * \varphi}{\partial X} dX @f$
 *
 * @param[in] displacements nodes displacements
 * @param[in] ghost_type compute the residual for _ghost or _not_ghost element
 */
// void Material::updateResidual(GhostType ghost_type) {
//   AKANTU_DEBUG_IN();

//   computeAllStresses(ghost_type);
//   assembleResidual(ghost_type);

//   AKANTU_DEBUG_OUT();
// }

/* -------------------------------------------------------------------------- */
void Material::assembleInternalForces(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = model.getSpatialDimension();

  if (!finite_deformation) {

    Array<Real> & internal_force =
        const_cast<Array<Real> &>(model.getInternalForce());

    Mesh & mesh = fem.getMesh();
    Mesh::type_iterator it =
        element_filter.firstType(spatial_dimension, ghost_type);
    Mesh::type_iterator last_type =
        element_filter.lastType(spatial_dimension, ghost_type);
    for (; it != last_type; ++it) {
      Array<UInt> & elem_filter = element_filter(*it, ghost_type);
      UInt nb_element = elem_filter.getSize();
      if (nb_element) {
        const Array<Real> & shapes_derivatives =
            fem.getShapesDerivatives(*it, ghost_type);

        UInt size_of_shapes_derivatives = shapes_derivatives.getNbComponent();
        UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(*it);
        UInt nb_quadrature_points = fem.getNbIntegrationPoints(*it, ghost_type);

        /// compute @f$\sigma \frac{\partial \varphi}{\partial X}@f$ by
        /// @f$\mathbf{B}^t \mathbf{\sigma}_q@f$
        Array<Real> * sigma_dphi_dx =
            new Array<Real>(nb_element * nb_quadrature_points,
                            size_of_shapes_derivatives, "sigma_x_dphi_/_dX");

        Array<Real> * shapesd_filtered =
            new Array<Real>(0, size_of_shapes_derivatives, "filtered shapesd");

        FEEngine::filterElementalData(mesh, shapes_derivatives,
                                      *shapesd_filtered, *it, ghost_type,
                                      elem_filter);

        Array<Real> & stress_vect = this->stress(*it, ghost_type);

        Array<Real>::matrix_iterator sigma =
            stress_vect.begin(spatial_dimension, spatial_dimension);
        Array<Real>::matrix_iterator B =
            shapesd_filtered->begin(spatial_dimension, nb_nodes_per_element);
        Array<Real>::matrix_iterator Bt_sigma_it =
            sigma_dphi_dx->begin(spatial_dimension, nb_nodes_per_element);

        for (UInt q = 0; q < nb_element * nb_quadrature_points;
             ++q, ++sigma, ++B, ++Bt_sigma_it)
          Bt_sigma_it->mul<false, false>(*sigma, *B);

        delete shapesd_filtered;

        /**
         * compute @f$\int \sigma  * \frac{\partial \varphi}{\partial X}dX@f$ by
         * @f$ \sum_q \mathbf{B}^t
         * \mathbf{\sigma}_q \overline w_q J_q@f$
         */
        Array<Real> * int_sigma_dphi_dx = new Array<Real>(
            nb_element, nb_nodes_per_element * spatial_dimension,
            "int_sigma_x_dphi_/_dX");

        fem.integrate(*sigma_dphi_dx, *int_sigma_dphi_dx,
                      size_of_shapes_derivatives, *it, ghost_type, elem_filter);
        delete sigma_dphi_dx;

        /// assemble
        model.getDOFManager().assembleElementalArrayLocalArray(
            *int_sigma_dphi_dx, internal_force, *it, ghost_type, -1,
            elem_filter);
        delete int_sigma_dphi_dx;
      }
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
 * @param[in] current_position nodes postition + displacements
 * @param[in] ghost_type compute the residual for _ghost or _not_ghost element
 */
void Material::computeAllStresses(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = model.getSpatialDimension();

  for (const auto & type :
       element_filter.elementTypes(spatial_dimension, ghost_type)) {
    Array<UInt> & elem_filter = element_filter(type, ghost_type);

    if (elem_filter.getSize() == 0)
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
      this->computeCauchyStress<1>(type, ghost_type);
      break;
    case 2:
      this->computeCauchyStress<2>(type, ghost_type);
      break;
    case 3:
      this->computeCauchyStress<3>(type, ghost_type);
      break;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void Material::computeCauchyStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Array<Real>::matrix_iterator gradu_it =
      this->gradu(el_type, ghost_type).begin(dim, dim);

  Array<Real>::matrix_iterator gradu_end =
      this->gradu(el_type, ghost_type).end(dim, dim);

  Array<Real>::matrix_iterator piola_it =
      this->piola_kirchhoff_2(el_type, ghost_type).begin(dim, dim);

  Array<Real>::matrix_iterator stress_it =
      this->stress(el_type, ghost_type).begin(dim, dim);

  Matrix<Real> F_tensor(dim, dim);

  for (; gradu_it != gradu_end; ++gradu_it, ++piola_it, ++stress_it) {
    Matrix<Real> & grad_u = *gradu_it;
    Matrix<Real> & piola = *piola_it;
    Matrix<Real> & sigma = *stress_it;

    gradUToF<dim>(grad_u, F_tensor);
    this->computeCauchyStressOnQuad<dim>(F_tensor, piola, sigma);
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
 * @param[in] current_position nodes postition + displacements
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
  if (elem_filter.getSize()) {

    const Array<Real> & shapes_derivatives =
        fem.getShapesDerivatives(type, ghost_type);

    Array<Real> & gradu_vect = gradu(type, ghost_type);

    UInt nb_element = elem_filter.getSize();
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    UInt nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);

    gradu_vect.resize(nb_quadrature_points * nb_element);

    fem.gradientOnIntegrationPoints(model.getDisplacement(), gradu_vect, dim,
                                    type, ghost_type, elem_filter);

    UInt tangent_size = getTangentStiffnessVoigtSize(dim);

    Array<Real> * tangent_stiffness_matrix = new Array<Real>(
        nb_element * nb_quadrature_points, tangent_size * tangent_size,
        "tangent_stiffness_matrix");

    tangent_stiffness_matrix->clear();

    computeTangentModuli(type, *tangent_stiffness_matrix, ghost_type);

    Array<Real> * shapesd_filtered = new Array<Real>(
        nb_element, dim * nb_nodes_per_element, "filtered shapesd");

    FEEngine::filterElementalData(fem.getMesh(), shapes_derivatives,
                                  *shapesd_filtered, type, ghost_type,
                                  elem_filter);

    /// compute @f$\mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
    UInt bt_d_b_size = dim * nb_nodes_per_element;

    Array<Real> * bt_d_b =
        new Array<Real>(nb_element * nb_quadrature_points,
                        bt_d_b_size * bt_d_b_size, "B^t*D*B");

    Matrix<Real> B(tangent_size, dim * nb_nodes_per_element);
    Matrix<Real> Bt_D(dim * nb_nodes_per_element, tangent_size);

    Array<Real>::matrix_iterator shapes_derivatives_filtered_it =
        shapesd_filtered->begin(dim, nb_nodes_per_element);

    Array<Real>::matrix_iterator Bt_D_B_it =
        bt_d_b->begin(dim * nb_nodes_per_element, dim * nb_nodes_per_element);

    Array<Real>::matrix_iterator D_it =
        tangent_stiffness_matrix->begin(tangent_size, tangent_size);

    Array<Real>::matrix_iterator D_end =
        tangent_stiffness_matrix->end(tangent_size, tangent_size);

    for (; D_it != D_end;
         ++D_it, ++Bt_D_B_it, ++shapes_derivatives_filtered_it) {
      Matrix<Real> & D = *D_it;
      Matrix<Real> & Bt_D_B = *Bt_D_B_it;

      VoigtHelper<dim>::transferBMatrixToSymVoigtBMatrix(
          *shapes_derivatives_filtered_it, B, nb_nodes_per_element);
      Bt_D.mul<true, false>(B, D);
      Bt_D_B.mul<false, false>(Bt_D, B);
    }

    delete tangent_stiffness_matrix;
    delete shapesd_filtered;

    /// compute @f$ k_e = \int_e \mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
    Array<Real> * K_e =
        new Array<Real>(nb_element, bt_d_b_size * bt_d_b_size, "K_e");

    fem.integrate(*bt_d_b, *K_e, bt_d_b_size * bt_d_b_size, type, ghost_type,
                  elem_filter);

    delete bt_d_b;

    model.getDOFManager().assembleElementalMatricesToMatrix(
        "K", "displacement", *K_e, type, ghost_type, _symmetric, elem_filter);
    delete K_e;
  }
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

  UInt nb_element = elem_filter.getSize();
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);

  // gradu_vect.resize(nb_quadrature_points * nb_element);

  // fem.gradientOnIntegrationPoints(model.getIncrement(), gradu_vect,
  //        dim, type, ghost_type, &elem_filter);

  Array<Real> * shapes_derivatives_filtered = new Array<Real>(
      nb_element * nb_quadrature_points, dim * nb_nodes_per_element,
      "shapes derivatives filtered");

  Array<Real>::const_matrix_iterator shapes_derivatives_it =
      shapes_derivatives.begin(spatial_dimension, nb_nodes_per_element);

  Array<Real>::matrix_iterator shapes_derivatives_filtered_it =
      shapes_derivatives_filtered->begin(spatial_dimension,
                                         nb_nodes_per_element);
  UInt * elem_filter_val = elem_filter.storage();
  for (UInt e = 0; e < nb_element; ++e, ++elem_filter_val)
    for (UInt q = 0; q < nb_quadrature_points;
         ++q, ++shapes_derivatives_filtered_it)
      *shapes_derivatives_filtered_it =
          shapes_derivatives_it[*elem_filter_val * nb_quadrature_points + q];

  /// compute @f$\mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  UInt bt_s_b_size = dim * nb_nodes_per_element;

  Array<Real> * bt_s_b = new Array<Real>(nb_element * nb_quadrature_points,
                                         bt_s_b_size * bt_s_b_size, "B^t*D*B");

  UInt piola_matrix_size = getCauchyStressMatrixSize(dim);

  Matrix<Real> B(piola_matrix_size, bt_s_b_size);
  Matrix<Real> Bt_S(bt_s_b_size, piola_matrix_size);
  Matrix<Real> S(piola_matrix_size, piola_matrix_size);

  shapes_derivatives_filtered_it = shapes_derivatives_filtered->begin(
      spatial_dimension, nb_nodes_per_element);

  Array<Real>::matrix_iterator Bt_S_B_it =
      bt_s_b->begin(bt_s_b_size, bt_s_b_size);
  Array<Real>::matrix_iterator Bt_S_B_end =
      bt_s_b->end(bt_s_b_size, bt_s_b_size);

  Array<Real>::matrix_iterator piola_it =
      piola_kirchhoff_2(type, ghost_type).begin(dim, dim);

  for (; Bt_S_B_it != Bt_S_B_end;
       ++Bt_S_B_it, ++shapes_derivatives_filtered_it, ++piola_it) {
    Matrix<Real> & Bt_S_B = *Bt_S_B_it;
    Matrix<Real> & Piola_kirchhoff_matrix = *piola_it;

    setCauchyStressMatrix<dim>(Piola_kirchhoff_matrix, S);
    VoigtHelper<dim>::transferBMatrixToBNL(*shapes_derivatives_filtered_it, B,
                                           nb_nodes_per_element);
    Bt_S.mul<true, false>(B, S);
    Bt_S_B.mul<false, false>(Bt_S, B);
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

  UInt nb_element = elem_filter.getSize();
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

  Array<Real>::const_matrix_iterator shapes_derivatives_it =
      shapes_derivatives.begin(spatial_dimension, nb_nodes_per_element);

  Array<Real>::matrix_iterator shapes_derivatives_filtered_it =
      shapes_derivatives_filtered->begin(spatial_dimension,
                                         nb_nodes_per_element);
  UInt * elem_filter_val = elem_filter.storage();
  for (UInt e = 0; e < nb_element; ++e, ++elem_filter_val)
    for (UInt q = 0; q < nb_quadrature_points;
         ++q, ++shapes_derivatives_filtered_it)
      *shapes_derivatives_filtered_it =
          shapes_derivatives_it[*elem_filter_val * nb_quadrature_points + q];

  /// compute @f$\mathbf{B}^t * \mathbf{D} * \mathbf{B}@f$
  UInt bt_d_b_size = dim * nb_nodes_per_element;

  Array<Real> * bt_d_b = new Array<Real>(nb_element * nb_quadrature_points,
                                         bt_d_b_size * bt_d_b_size, "B^t*D*B");

  Matrix<Real> B(tangent_size, dim * nb_nodes_per_element);
  Matrix<Real> B2(tangent_size, dim * nb_nodes_per_element);
  Matrix<Real> Bt_D(dim * nb_nodes_per_element, tangent_size);

  shapes_derivatives_filtered_it = shapes_derivatives_filtered->begin(
      spatial_dimension, nb_nodes_per_element);

  Array<Real>::matrix_iterator Bt_D_B_it =
      bt_d_b->begin(dim * nb_nodes_per_element, dim * nb_nodes_per_element);

  Array<Real>::matrix_iterator grad_u_it = gradu_vect.begin(dim, dim);

  Array<Real>::matrix_iterator D_it =
      tangent_stiffness_matrix->begin(tangent_size, tangent_size);
  Array<Real>::matrix_iterator D_end =
      tangent_stiffness_matrix->end(tangent_size, tangent_size);

  for (; D_it != D_end;
       ++D_it, ++Bt_D_B_it, ++shapes_derivatives_filtered_it, ++grad_u_it) {
    Matrix<Real> & grad_u = *grad_u_it;
    Matrix<Real> & D = *D_it;
    Matrix<Real> & Bt_D_B = *Bt_D_B_it;

    // transferBMatrixToBL1<dim > (*shapes_derivatives_filtered_it, B,
    // nb_nodes_per_element);
    VoigtHelper<dim>::transferBMatrixToSymVoigtBMatrix(
        *shapes_derivatives_filtered_it, B, nb_nodes_per_element);
    VoigtHelper<dim>::transferBMatrixToBL2(*shapes_derivatives_filtered_it,
                                           grad_u, B2, nb_nodes_per_element);
    B += B2;
    Bt_D.mul<true, false>(B, D);
    Bt_D_B.mul<false, false>(Bt_D, B);
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
  for (auto type : element_filter.elementTypes(spatial_dimension, ghost_type)) {
    const Array<Real> & shapes_derivatives =
        fem.getShapesDerivatives(type, ghost_type);

    Array<UInt> & elem_filter = element_filter(type, ghost_type);
    if (elem_filter.getSize() == 0)
      continue;
    UInt size_of_shapes_derivatives = shapes_derivatives.getNbComponent();
    UInt nb_element = elem_filter.getSize();
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    UInt nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);

    Array<Real> * shapesd_filtered = new Array<Real>(
        nb_element, size_of_shapes_derivatives, "filtered shapesd");

    FEEngine::filterElementalData(mesh, shapes_derivatives, *shapesd_filtered,
                                  type, ghost_type, elem_filter);

    Array<Real>::matrix_iterator shapes_derivatives_filtered_it =
        shapesd_filtered->begin(dim, nb_nodes_per_element);

    // Set stress vectors

    UInt stress_size = getTangentStiffnessVoigtSize(dim);

    // Set matrices B and BNL*
    UInt bt_s_size = dim * nb_nodes_per_element;

    Array<Real> * bt_s =
        new Array<Real>(nb_element * nb_quadrature_points, bt_s_size, "B^t*S");

    Array<Real>::matrix_iterator grad_u_it =
        this->gradu(type, ghost_type).begin(dim, dim);

    Array<Real>::matrix_iterator grad_u_end =
        this->gradu(type, ghost_type).end(dim, dim);

    Array<Real>::matrix_iterator stress_it =
        this->piola_kirchhoff_2(type, ghost_type).begin(dim, dim);

    shapes_derivatives_filtered_it =
        shapesd_filtered->begin(dim, nb_nodes_per_element);

    Array<Real>::matrix_iterator bt_s_it = bt_s->begin(bt_s_size, 1);

    Matrix<Real> S_vect(stress_size, 1);
    Matrix<Real> B_tensor(stress_size, bt_s_size);
    Matrix<Real> B2_tensor(stress_size, bt_s_size);

    for (; grad_u_it != grad_u_end; ++grad_u_it, ++stress_it,
                                    ++shapes_derivatives_filtered_it,
                                    ++bt_s_it) {
      Matrix<Real> & grad_u = *grad_u_it;
      Matrix<Real> & r_it = *bt_s_it;
      Matrix<Real> & S_it = *stress_it;

      setCauchyStressArray<dim>(S_it, S_vect);
      VoigtHelper<dim>::transferBMatrixToSymVoigtBMatrix(
          *shapes_derivatives_filtered_it, B_tensor, nb_nodes_per_element);
      VoigtHelper<dim>::transferBMatrixToBL2(*shapes_derivatives_filtered_it,
                                             grad_u, B2_tensor,
                                             nb_nodes_per_element);

      B_tensor += B2_tensor;

      r_it.mul<true, false>(B_tensor, S_vect);
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
void Material::computeAllStressesFromTangentModuli(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = model.getSpatialDimension();

  for (auto type : element_filter.elementTypes(spatial_dimension, ghost_type)) {
    switch (spatial_dimension) {
    case 1: {
      computeAllStressesFromTangentModuli<1>(type, ghost_type);
      break;
    }
    case 2: {
      computeAllStressesFromTangentModuli<2>(type, ghost_type);
      break;
    }
    case 3: {
      computeAllStressesFromTangentModuli<3>(type, ghost_type);
      break;
    }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void Material::computeAllStressesFromTangentModuli(const ElementType & type,
                                                   GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  const Array<Real> & shapes_derivatives =
      fem.getShapesDerivatives(type, ghost_type);
  Array<UInt> & elem_filter = element_filter(type, ghost_type);
  Array<Real> & gradu_vect = gradu(type, ghost_type);

  UInt nb_element = elem_filter.getSize();
  if (nb_element) {
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    UInt nb_quadrature_points = fem.getNbIntegrationPoints(type, ghost_type);

    gradu_vect.resize(nb_quadrature_points * nb_element);

    Array<Real> & disp = model.getDisplacement();

    fem.gradientOnIntegrationPoints(disp, gradu_vect, dim, type, ghost_type,
                                    elem_filter);

    UInt tangent_moduli_size = getTangentStiffnessVoigtSize(dim);

    Array<Real> * tangent_moduli_tensors = new Array<Real>(
        nb_element * nb_quadrature_points,
        tangent_moduli_size * tangent_moduli_size, "tangent_moduli_tensors");

    tangent_moduli_tensors->clear();
    computeTangentModuli(type, *tangent_moduli_tensors, ghost_type);

    Array<Real> * shapesd_filtered = new Array<Real>(
        nb_element, dim * nb_nodes_per_element, "filtered shapesd");

    FEEngine::filterElementalData(fem.getMesh(), shapes_derivatives,
                                  *shapesd_filtered, type, ghost_type,
                                  elem_filter);

    Array<Real> filtered_u(nb_element,
                           nb_nodes_per_element * spatial_dimension);

    FEEngine::extractNodalToElementField(fem.getMesh(), disp, filtered_u, type,
                                         ghost_type, elem_filter);

    /// compute @f$\mathbf{D} \mathbf{B} \mathbf{u}@f$
    Array<Real>::matrix_iterator shapes_derivatives_filtered_it =
        shapesd_filtered->begin(dim, nb_nodes_per_element);

    Array<Real>::matrix_iterator D_it =
        tangent_moduli_tensors->begin(tangent_moduli_size, tangent_moduli_size);
    Array<Real>::matrix_iterator sigma_it =
        stress(type, ghost_type).begin(spatial_dimension, spatial_dimension);
    Array<Real>::vector_iterator u_it =
        filtered_u.begin(spatial_dimension * nb_nodes_per_element);

    Matrix<Real> B(tangent_moduli_size,
                   spatial_dimension * nb_nodes_per_element);
    Vector<Real> Bu(tangent_moduli_size);
    Vector<Real> DBu(tangent_moduli_size);

    for (UInt e = 0; e < nb_element; ++e, ++u_it) {
      for (UInt q = 0; q < nb_quadrature_points;
           ++q, ++D_it, ++shapes_derivatives_filtered_it, ++sigma_it) {
        Vector<Real> & u = *u_it;
        Matrix<Real> & sigma = *sigma_it;
        Matrix<Real> & D = *D_it;

        VoigtHelper<dim>::transferBMatrixToSymVoigtBMatrix(
            *shapes_derivatives_filtered_it, B, nb_nodes_per_element);

        Bu.mul<false>(B, u);
        DBu.mul<false>(D, Bu);

        // Voigt notation to full symmetric tensor
        for (UInt i = 0; i < dim; ++i)
          sigma(i, i) = DBu(i);
        if (dim == 2) {
          sigma(0, 1) = sigma(1, 0) = DBu(2);
        } else if (dim == 3) {
          sigma(1, 2) = sigma(2, 1) = DBu(3);
          sigma(0, 2) = sigma(2, 0) = DBu(4);
          sigma(0, 1) = sigma(1, 0) = DBu(5);
        }
      }
    }
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
void Material::computePotentialEnergy(ElementType, GhostType) {
  AKANTU_DEBUG_IN();

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
Real Material::getEnergy(std::string type) {
  AKANTU_DEBUG_IN();
  if (type == "potential")
    return getPotentialEnergy();
  AKANTU_DEBUG_OUT();
  return 0.;
}

/* -------------------------------------------------------------------------- */
Real Material::getEnergy(std::string energy_id, ElementType type, UInt index) {
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
    UInt nb_element = elem_fil.getSize();
    UInt nb_element_full = this->model.getMesh().getNbElement(type, ghost_type);
    UInt nb_interpolation_points_per_elem =
        by_elem_res.getSize() / nb_element_full;

    const Array<Element> & facet_to_element =
        mesh_facets.getSubelementToElement(type, ghost_type);
    ElementType type_facet = Mesh::getFacetType(type);
    UInt nb_facet_per_elem = facet_to_element.getNbComponent();
    UInt nb_quad_per_facet =
        nb_interpolation_points_per_elem / nb_facet_per_elem;
    Element element_for_comparison(type, 0, ghost_type);
    const Array<std::vector<Element>> * element_to_facet = NULL;
    GhostType current_ghost_type = _casper;
    Array<Real> * result_vec = NULL;

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
const Array<T> &
Material::getArray(__attribute__((unused)) const ID & vect_id,
                   __attribute__((unused)) const ElementType & type,
                   __attribute__((unused)) const GhostType & ghost_type) const {
  AKANTU_DEBUG_TO_IMPLEMENT();
  return NULL;
}

/* -------------------------------------------------------------------------- */
template <typename T>
Array<T> & Material::getArray(__attribute__((unused)) const ID & vect_id,
                              __attribute__((unused)) const ElementType & type,
                              __attribute__((unused))
                              const GhostType & ghost_type) {
  AKANTU_DEBUG_TO_IMPLEMENT();
  return NULL;
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
                                            << vect_id << "(" << fvect_id
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
                                            << vect_id << "(" << fvect_id
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
                                            << vect_id << "(" << fvect_id
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
const InternalField<T> & Material::getInternal(__attribute__((unused))
                                               const ID & int_id) const {
  AKANTU_DEBUG_TO_IMPLEMENT();
  return NULL;
}

/* -------------------------------------------------------------------------- */
template <typename T>
InternalField<T> & Material::getInternal(__attribute__((unused))
                                         const ID & int_id) {
  AKANTU_DEBUG_TO_IMPLEMENT();
  return NULL;
}

/* -------------------------------------------------------------------------- */
template <>
const InternalField<Real> & Material::getInternal(const ID & int_id) const {
  std::map<ID, InternalField<Real> *>::const_iterator it =
      internal_vectors_real.find(getID() + ":" + int_id);
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
  std::map<ID, InternalField<Real> *>::iterator it =
      internal_vectors_real.find(getID() + ":" + int_id);
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
  std::map<ID, InternalField<UInt> *>::const_iterator it =
      internal_vectors_uint.find(getID() + ":" + int_id);
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
  std::map<ID, InternalField<UInt> *>::iterator it =
      internal_vectors_uint.find(getID() + ":" + int_id);
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
  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {
    GhostType ghost_type = *gt;
    element.ghost_type = ghost_type;

    ElementTypeMapArray<UInt>::type_iterator it =
        element_filter.firstType(_all_dimensions, ghost_type, _ek_not_defined);
    ElementTypeMapArray<UInt>::type_iterator end =
        element_filter.lastType(_all_dimensions, ghost_type, _ek_not_defined);

    for (; it != end; ++it) {
      ElementType type = *it;
      element.type = type;

      Array<UInt> & elem_filter = this->element_filter(type, ghost_type);
      Array<UInt> & mat_loc_num =
          this->model.getMaterialLocalNumbering(type, ghost_type);

      if (!material_local_new_numbering.exists(type, ghost_type))
        material_local_new_numbering.alloc(elem_filter.getSize(), 1, type,
                                           ghost_type);
      Array<UInt> & mat_renumbering =
          material_local_new_numbering(type, ghost_type);

      UInt nb_element = elem_filter.getSize();
      element.kind = (*el_begin).kind;
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

      elem_filter.resize(elem_filter_tmp.getSize());
      elem_filter.copy(elem_filter_tmp);
    }
  }

  for (std::map<ID, InternalField<Real> *>::iterator it =
           internal_vectors_real.begin();
       it != internal_vectors_real.end(); ++it)
    it->second->removeIntegrationPoints(material_local_new_numbering);

  for (std::map<ID, InternalField<UInt> *>::iterator it =
           internal_vectors_uint.begin();
       it != internal_vectors_uint.end(); ++it)
    it->second->removeIntegrationPoints(material_local_new_numbering);

  for (std::map<ID, InternalField<bool> *>::iterator it =
           internal_vectors_bool.begin();
       it != internal_vectors_bool.end(); ++it)
    it->second->removeIntegrationPoints(material_local_new_numbering);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Material::resizeInternals() {
  AKANTU_DEBUG_IN();
  for (std::map<ID, InternalField<Real> *>::iterator it =
           internal_vectors_real.begin();
       it != internal_vectors_real.end(); ++it)
    it->second->resize();

  for (std::map<ID, InternalField<UInt> *>::iterator it =
           internal_vectors_uint.begin();
       it != internal_vectors_uint.end(); ++it)
    it->second->resize();

  for (std::map<ID, InternalField<bool> *>::iterator it =
           internal_vectors_bool.begin();
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
    __attribute__((unused)) const RemovedElementsEvent & event) {
  UInt my_num = model.getInternalIndexFromID(getID());

  ElementTypeMapArray<UInt> material_local_new_numbering(
      "remove mat filter elem", getID(), getMemoryID());

  Array<Element>::const_iterator<Element> el_begin = element_list.begin();
  Array<Element>::const_iterator<Element> el_end = element_list.end();

  for (ghost_type_t::iterator g = ghost_type_t::begin();
       g != ghost_type_t::end(); ++g) {
    GhostType gt = *g;

    ElementTypeMapArray<UInt>::type_iterator it =
        new_numbering.firstType(_all_dimensions, gt, _ek_not_defined);
    ElementTypeMapArray<UInt>::type_iterator end =
        new_numbering.lastType(_all_dimensions, gt, _ek_not_defined);
    for (; it != end; ++it) {
      ElementType type = *it;
      if (element_filter.exists(type, gt) &&
          element_filter(type, gt).getSize()) {
        Array<UInt> & elem_filter = element_filter(type, gt);

        Array<UInt> & mat_indexes = this->model.getMaterialByElement(*it, gt);
        Array<UInt> & mat_loc_num =
            this->model.getMaterialLocalNumbering(*it, gt);
        UInt nb_element = this->model.getMesh().getNbElement(type, gt);

        // all materials will resize of the same size...
        mat_indexes.resize(nb_element);
        mat_loc_num.resize(nb_element);

        if (!material_local_new_numbering.exists(type, gt))
          material_local_new_numbering.alloc(elem_filter.getSize(), 1, type,
                                             gt);

        Array<UInt> & mat_renumbering = material_local_new_numbering(type, gt);
        const Array<UInt> & renumbering = new_numbering(type, gt);
        Array<UInt> elem_filter_tmp;
        UInt ni = 0;
        Element el;
        el.type = type;
        el.ghost_type = gt;
        el.kind = Mesh::getKind(type);
        for (UInt i = 0; i < elem_filter.getSize(); ++i) {
          el.element = elem_filter(i);
          if (std::find(el_begin, el_end, el) == el_end) {
            UInt new_el = renumbering(el.element);
            AKANTU_DEBUG_ASSERT(
                new_el != UInt(-1),
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

        elem_filter.resize(elem_filter_tmp.getSize());
        elem_filter.copy(elem_filter_tmp);
      }
    }
  }

  for (std::map<ID, InternalField<Real> *>::iterator it =
           internal_vectors_real.begin();
       it != internal_vectors_real.end(); ++it)
    it->second->removeIntegrationPoints(material_local_new_numbering);

  for (std::map<ID, InternalField<UInt> *>::iterator it =
           internal_vectors_uint.begin();
       it != internal_vectors_uint.end(); ++it)
    it->second->removeIntegrationPoints(material_local_new_numbering);

  for (std::map<ID, InternalField<bool> *>::iterator it =
           internal_vectors_bool.begin();
       it != internal_vectors_bool.end(); ++it)
    it->second->removeIntegrationPoints(material_local_new_numbering);
}

/* -------------------------------------------------------------------------- */
void Material::onBeginningSolveStep(__attribute__((unused))
                                    const AnalysisMethod & method) {
  this->savePreviousState();
}

/* -------------------------------------------------------------------------- */
void Material::onEndSolveStep(__attribute__((unused))
                              const AnalysisMethod & method) {
  ElementTypeMapArray<UInt>::type_iterator it = this->element_filter.firstType(
      _all_dimensions, _not_ghost, _ek_not_defined);
  ElementTypeMapArray<UInt>::type_iterator end =
      element_filter.lastType(_all_dimensions, _not_ghost, _ek_not_defined);

  for (; it != end; ++it) {
    this->updateEnergies(*it, _not_ghost);
  }
}
/* -------------------------------------------------------------------------- */
void Material::onDamageIteration() { this->savePreviousState(); }

/* -------------------------------------------------------------------------- */
void Material::onDamageUpdate() {
  ElementTypeMapArray<UInt>::type_iterator it = this->element_filter.firstType(
      _all_dimensions, _not_ghost, _ek_not_defined);
  ElementTypeMapArray<UInt>::type_iterator end =
      element_filter.lastType(_all_dimensions, _not_ghost, _ek_not_defined);

  for (; it != end; ++it) {
    this->updateEnergiesAfterDamage(*it, _not_ghost);
  }
}

/* -------------------------------------------------------------------------- */
void Material::onDump() {
  if (this->isFiniteDeformation())
    this->computeAllCauchyStresses(_not_ghost);
}

/* -------------------------------------------------------------------------- */
void Material::printself(std::ostream & stream, int indent) const {
  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
    ;

  std::string type = getID().substr(getID().find_last_of(":") + 1);

  stream << space << "Material " << type << " [" << std::endl;
  Parsable::printself(stream, indent);
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
/// extrapolate internal values
void Material::extrapolateInternal(const ID & id, const Element & element,
                                   __attribute__((unused))
                                   const Matrix<Real> & point,
                                   Matrix<Real> & extrapolated) {
  if (this->isInternal<Real>(id, element.kind)) {
    UInt nb_element =
        this->element_filter(element.type, element.ghost_type).getSize();
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

  ElementTypeMapArray<UInt>::type_iterator it = this->element_filter.firstType(
      _all_dimensions, _not_ghost, _ek_not_defined);
  ElementTypeMapArray<UInt>::type_iterator end =
      element_filter.lastType(_all_dimensions, _not_ghost, _ek_not_defined);

  for (; it != end; ++it) {
    ElementType type = *it;
    if (!element_filter(type, ghost_type).getSize())
      continue;
    Array<Real>::matrix_iterator eigen_it =
        this->eigengradu(type, ghost_type)
            .begin(spatial_dimension, spatial_dimension);
    Array<Real>::matrix_iterator eigen_end =
        this->eigengradu(type, ghost_type)
            .end(spatial_dimension, spatial_dimension);
    for (; eigen_it != eigen_end; ++eigen_it) {
      Matrix<Real> & current_eigengradu = *eigen_it;
      current_eigengradu = prescribed_eigen_grad_u;
    }
  }
}

} // namespace akantu
