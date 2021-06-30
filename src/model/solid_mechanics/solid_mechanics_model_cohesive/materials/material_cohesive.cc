/**
 * @file   material_cohesive.cc
 *
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Seyedeh Mohadeseh Taheri Mousavi <mohadeseh.taherimousavi@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Feb 22 2012
 * @date last modification: Mon Feb 19 2018
 *
 * @brief  Specialization of the material class for cohesive elements
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
#include "material_cohesive.hh"
#include "aka_random_generator.hh"
#include "dof_synchronizer.hh"
#include "fe_engine_template.hh"
#include "integrator_gauss.hh"
#include "shape_cohesive.hh"
#include "solid_mechanics_model_cohesive.hh"
#include "sparse_matrix.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
MaterialCohesive::MaterialCohesive(SolidMechanicsModel & model, const ID & id)
    : Material(model, id),
      facet_filter("facet_filter", id),
      fem_cohesive(
          model.getFEEngineClass<MyFEEngineCohesiveType>("CohesiveFEEngine")),
      reversible_energy("reversible_energy", *this),
      total_energy("total_energy", *this), opening("opening", *this),
      tractions("tractions", *this),
      contact_tractions("contact_tractions", *this),
      contact_opening("contact_opening", *this), delta_max("delta max", *this),
      use_previous_delta_max(false), use_previous_opening(false),
      damage("damage", *this), sigma_c("sigma_c", *this),
      normal(0, spatial_dimension, "normal") {

  AKANTU_DEBUG_IN();

  this->model = dynamic_cast<SolidMechanicsModelCohesive *>(&model);

  this->registerParam("sigma_c", sigma_c, _pat_parsable | _pat_readable,
                      "Critical stress");
  this->registerParam("delta_c", delta_c, Real(0.),
                      _pat_parsable | _pat_readable, "Critical displacement");

  this->element_filter.initialize(this->model->getMesh(),
                                  _spatial_dimension = spatial_dimension,
                                  _element_kind = _ek_cohesive);
  // this->model->getMesh().initElementTypeMapArray(
  //     this->element_filter, 1, spatial_dimension, false, _ek_cohesive);

  if (this->model->getIsExtrinsic()) {
    this->facet_filter.initialize(this->model->getMeshFacets(),
                                  _spatial_dimension = spatial_dimension - 1,
                                  _element_kind = _ek_regular);
  }
  // this->model->getMeshFacets().initElementTypeMapArray(facet_filter, 1,
  //                                                      spatial_dimension -
  //                                                      1);

  this->reversible_energy.initialize(1);
  this->total_energy.initialize(1);

  this->tractions.initialize(spatial_dimension);
  this->tractions.initializeHistory();

  this->contact_tractions.initialize(spatial_dimension);
  this->contact_opening.initialize(spatial_dimension);

  this->opening.initialize(spatial_dimension);
  this->opening.initializeHistory();

  this->delta_max.initialize(1);
  this->damage.initialize(1);

  if (this->model->getIsExtrinsic()) {
    this->sigma_c.initialize(1);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
MaterialCohesive::~MaterialCohesive() = default;

/* -------------------------------------------------------------------------- */
void MaterialCohesive::initMaterial() {
  AKANTU_DEBUG_IN();
  Material::initMaterial();
  if (this->use_previous_delta_max) {
    this->delta_max.initializeHistory();
  }
  if (this->use_previous_opening) {
    this->opening.initializeHistory();
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialCohesive::assembleInternalForces(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

#if defined(AKANTU_DEBUG_TOOLS)
  debug::element_manager.printData(debug::_dm_material_cohesive,
                                   "Cohesive Tractions", tractions);
#endif

  auto & internal_force = const_cast<Array<Real> &>(model->getInternalForce());

  for (auto type : element_filter.elementTypes(spatial_dimension, ghost_type,
                                               _ek_cohesive)) {
    auto & elem_filter = element_filter(type, ghost_type);
    UInt nb_element = elem_filter.size();
    if (nb_element == 0) {
      continue;
    }

    const auto & shapes = fem_cohesive.getShapes(type, ghost_type);
    auto & traction = tractions(type, ghost_type);
    auto & contact_traction = contact_tractions(type, ghost_type);

    UInt size_of_shapes = shapes.getNbComponent();
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
    UInt nb_quadrature_points =
        fem_cohesive.getNbIntegrationPoints(type, ghost_type);

    /// compute @f$t_i N_a@f$

    auto * traction_cpy = new Array<Real>(nb_element * nb_quadrature_points,
                                          spatial_dimension * size_of_shapes);

    auto traction_it = traction.begin(spatial_dimension, 1);
    auto contact_traction_it = contact_traction.begin(spatial_dimension, 1);
    auto shapes_filtered_begin = shapes.begin(1, size_of_shapes);
    auto traction_cpy_it =
        traction_cpy->begin(spatial_dimension, size_of_shapes);

    Matrix<Real> traction_tmp(spatial_dimension, 1);

    for (UInt el = 0; el < nb_element; ++el) {

      UInt current_quad = elem_filter(el) * nb_quadrature_points;

      for (UInt q = 0; q < nb_quadrature_points; ++q, ++traction_it,
                ++contact_traction_it, ++current_quad, ++traction_cpy_it) {

        const Matrix<Real> & shapes_filtered =
            shapes_filtered_begin[current_quad];

        traction_tmp.copy(*traction_it);
        traction_tmp += *contact_traction_it;

        traction_cpy_it->mul<false, false>(traction_tmp, shapes_filtered);
      }
    }

    /**
     * compute @f$\int t \cdot N\, dS@f$ by  @f$ \sum_q \mathbf{N}^t
     * \mathbf{t}_q \overline w_q J_q@f$
     */
    auto * partial_int_t_N = new Array<Real>(
        nb_element, spatial_dimension * size_of_shapes, "int_t_N");

    fem_cohesive.integrate(*traction_cpy, *partial_int_t_N,
                           spatial_dimension * size_of_shapes, type, ghost_type,
                           elem_filter);

    delete traction_cpy;

    auto * int_t_N = new Array<Real>(
        nb_element, 2 * spatial_dimension * size_of_shapes, "int_t_N");

    Real * int_t_N_val = int_t_N->storage();
    Real * partial_int_t_N_val = partial_int_t_N->storage();
    for (UInt el = 0; el < nb_element; ++el) {
      std::copy_n(partial_int_t_N_val, size_of_shapes * spatial_dimension,
                  int_t_N_val);
      std::copy_n(partial_int_t_N_val, size_of_shapes * spatial_dimension,
                  int_t_N_val + size_of_shapes * spatial_dimension);

      for (UInt n = 0; n < size_of_shapes * spatial_dimension; ++n) {
        int_t_N_val[n] *= -1.;
      }

      int_t_N_val += nb_nodes_per_element * spatial_dimension;
      partial_int_t_N_val += size_of_shapes * spatial_dimension;
    }

    delete partial_int_t_N;

    /// assemble
    model->getDOFManager().assembleElementalArrayLocalArray(
        *int_t_N, internal_force, type, ghost_type, 1, elem_filter);

    delete int_t_N;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialCohesive::assembleStiffnessMatrix(GhostType ghost_type) {

  AKANTU_DEBUG_IN();

  for (auto type : element_filter.elementTypes(spatial_dimension, ghost_type,
                                               _ek_cohesive)) {
    UInt nb_quadrature_points =
        fem_cohesive.getNbIntegrationPoints(type, ghost_type);
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

    const Array<Real> & shapes = fem_cohesive.getShapes(type, ghost_type);
    Array<UInt> & elem_filter = element_filter(type, ghost_type);

    UInt nb_element = elem_filter.size();

    if (nb_element == 0U) {
      continue;
    }

    UInt size_of_shapes = shapes.getNbComponent();

    auto * shapes_filtered = new Array<Real>(nb_element * nb_quadrature_points,
                                             size_of_shapes, "filtered shapes");

    Real * shapes_filtered_val = shapes_filtered->storage();
    UInt * elem_filter_val = elem_filter.storage();

    for (UInt el = 0; el < nb_element; ++el) {
      auto * shapes_val = shapes.storage() + elem_filter_val[el] *
                                                 size_of_shapes *
                                                 nb_quadrature_points;
      memcpy(shapes_filtered_val, shapes_val,
             size_of_shapes * nb_quadrature_points * sizeof(Real));
      shapes_filtered_val += size_of_shapes * nb_quadrature_points;
    }

    Matrix<Real> A(spatial_dimension * size_of_shapes,
                   spatial_dimension * nb_nodes_per_element);

    for (UInt i = 0; i < spatial_dimension * size_of_shapes; ++i) {
      A(i, i) = 1;
      A(i, i + spatial_dimension * size_of_shapes) = -1;
    }

    /// get the tangent matrix @f$\frac{\partial{(t/\delta)}}{\partial{\delta}}
    /// @f$
    auto * tangent_stiffness_matrix = new Array<Real>(
        nb_element * nb_quadrature_points,
        spatial_dimension * spatial_dimension, "tangent_stiffness_matrix");

    //    Array<Real> * normal = new Array<Real>(nb_element *
    //    nb_quadrature_points, spatial_dimension, "normal");
    normal.resize(nb_quadrature_points);

    computeNormal(model->getCurrentPosition(), normal, type, ghost_type);

    /// compute openings @f$\mathbf{\delta}@f$
    // computeOpening(model->getDisplacement(), opening(type, ghost_type), type,
    // ghost_type);

    tangent_stiffness_matrix->zero();

    computeTangentTraction(type, *tangent_stiffness_matrix, normal, ghost_type);

    // delete normal;

    UInt size_at_nt_d_n_a = spatial_dimension * nb_nodes_per_element *
                            spatial_dimension * nb_nodes_per_element;
    auto * at_nt_d_n_a = new Array<Real>(nb_element * nb_quadrature_points,
                                         size_at_nt_d_n_a, "A^t*N^t*D*N*A");

    Array<Real>::iterator<Vector<Real>> shapes_filt_it =
        shapes_filtered->begin(size_of_shapes);

    Array<Real>::matrix_iterator D_it =
        tangent_stiffness_matrix->begin(spatial_dimension, spatial_dimension);

    Array<Real>::matrix_iterator At_Nt_D_N_A_it =
        at_nt_d_n_a->begin(spatial_dimension * nb_nodes_per_element,
                           spatial_dimension * nb_nodes_per_element);

    Array<Real>::matrix_iterator At_Nt_D_N_A_end =
        at_nt_d_n_a->end(spatial_dimension * nb_nodes_per_element,
                         spatial_dimension * nb_nodes_per_element);

    Matrix<Real> N(spatial_dimension, spatial_dimension * size_of_shapes);
    Matrix<Real> N_A(spatial_dimension,
                     spatial_dimension * nb_nodes_per_element);
    Matrix<Real> D_N_A(spatial_dimension,
                       spatial_dimension * nb_nodes_per_element);

    for (; At_Nt_D_N_A_it != At_Nt_D_N_A_end;
         ++At_Nt_D_N_A_it, ++D_it, ++shapes_filt_it) {
      N.zero();
      /**
       * store  the   shapes  in  voigt   notations  matrix  @f$\mathbf{N}  =
       * \begin{array}{cccccc} N_0(\xi) & 0 & N_1(\xi)  &0 & N_2(\xi) & 0 \\
       * 0 & * N_0(\xi)& 0 &N_1(\xi)& 0 & N_2(\xi) \end{array} @f$
       **/
      for (UInt i = 0; i < spatial_dimension; ++i) {
        for (UInt n = 0; n < size_of_shapes; ++n) {
          N(i, i + spatial_dimension * n) = (*shapes_filt_it)(n);
        }
      }

      /**
       * compute stiffness matrix  @f$   \mathbf{K}    =    \delta \mathbf{U}^T
       * \int_{\Gamma_c}    {\mathbf{P}^t \frac{\partial{\mathbf{t}}}
       *{\partial{\delta}}
       * \mathbf{P} d\Gamma \Delta \mathbf{U}}  @f$
       **/
      N_A.mul<false, false>(N, A);
      D_N_A.mul<false, false>(*D_it, N_A);
      (*At_Nt_D_N_A_it).mul<true, false>(D_N_A, N_A);
    }

    delete tangent_stiffness_matrix;
    delete shapes_filtered;

    auto * K_e = new Array<Real>(nb_element, size_at_nt_d_n_a, "K_e");

    fem_cohesive.integrate(*at_nt_d_n_a, *K_e, size_at_nt_d_n_a, type,
                           ghost_type, elem_filter);

    delete at_nt_d_n_a;

    model->getDOFManager().assembleElementalMatricesToMatrix(
        "K", "displacement", *K_e, type, ghost_type, _unsymmetric, elem_filter);

    delete K_e;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- *
 * Compute traction from displacements
 *
 * @param[in] ghost_type compute the residual for _ghost or _not_ghost element
 */
void MaterialCohesive::computeTraction(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

#if defined(AKANTU_DEBUG_TOOLS)
  debug::element_manager.printData(debug::_dm_material_cohesive,
                                   "Cohesive Openings", opening);
#endif

  for (const auto & type : element_filter.elementTypes(
           spatial_dimension, ghost_type, _ek_cohesive)) {
    Array<UInt> & elem_filter = element_filter(type, ghost_type);

    UInt nb_element = elem_filter.size();
    if (nb_element == 0) {
      continue;
    }

    UInt nb_quadrature_points =
        nb_element * fem_cohesive.getNbIntegrationPoints(type, ghost_type);

    normal.resize(nb_quadrature_points);

    /// compute normals @f$\mathbf{n}@f$
    computeNormal(model->getCurrentPosition(), normal, type, ghost_type);

    /// compute openings @f$\mathbf{\delta}@f$
    computeOpening(model->getDisplacement(), opening(type, ghost_type), type,
                   ghost_type);

    /// compute traction @f$\mathbf{t}@f$
    computeTraction(normal, type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialCohesive::computeNormal(const Array<Real> & position,
                                     Array<Real> & normal, ElementType type,
                                     GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto & fem_cohesive =
      this->model->getFEEngineClass<MyFEEngineCohesiveType>("CohesiveFEEngine");

  normal.zero();

#define COMPUTE_NORMAL(type)                                                   \
  fem_cohesive.getShapeFunctions()                                             \
      .computeNormalsOnIntegrationPoints<type, CohesiveReduceFunctionMean>(    \
          position, normal, ghost_type, element_filter(type, ghost_type));

  AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(COMPUTE_NORMAL);
#undef COMPUTE_NORMAL

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialCohesive::computeOpening(const Array<Real> & displacement,
                                      Array<Real> & opening, ElementType type,
                                      GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto & fem_cohesive =
      this->model->getFEEngineClass<MyFEEngineCohesiveType>("CohesiveFEEngine");

#define COMPUTE_OPENING(type)                                                  \
  fem_cohesive.getShapeFunctions()                                             \
      .interpolateOnIntegrationPoints<type, CohesiveReduceFunctionOpening>(    \
          displacement, opening, spatial_dimension, ghost_type,                \
          element_filter(type, ghost_type));

  AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(COMPUTE_OPENING);
#undef COMPUTE_OPENING

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialCohesive::updateEnergies(ElementType type) {
  AKANTU_DEBUG_IN();

  if (Mesh::getKind(type) != _ek_cohesive) {
    return;
  }

  Vector<Real> b(spatial_dimension);
  Vector<Real> h(spatial_dimension);
  auto erev = reversible_energy(type).begin();
  auto etot = total_energy(type).begin();
  auto traction_it = tractions(type).begin(spatial_dimension);
  auto traction_old_it = tractions.previous(type).begin(spatial_dimension);
  auto opening_it = opening(type).begin(spatial_dimension);
  auto opening_old_it = opening.previous(type).begin(spatial_dimension);

  auto traction_end = tractions(type).end(spatial_dimension);

  /// loop on each quadrature point
  for (; traction_it != traction_end; ++traction_it, ++traction_old_it,
                                      ++opening_it, ++opening_old_it, ++erev,
                                      ++etot) {
    /// trapezoidal integration
    b = *opening_it;
    b -= *opening_old_it;

    h = *traction_old_it;
    h += *traction_it;

    *etot += .5 * b.dot(h);
    *erev = .5 * traction_it->dot(*opening_it);
  }

  /// update old values

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Real MaterialCohesive::getReversibleEnergy() {
  AKANTU_DEBUG_IN();
  Real erev = 0.;

  /// integrate reversible energy for each type of elements
  for (const auto & type : element_filter.elementTypes(
           spatial_dimension, _not_ghost, _ek_cohesive)) {
    erev +=
        fem_cohesive.integrate(reversible_energy(type, _not_ghost), type,
                               _not_ghost, element_filter(type, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return erev;
}

/* -------------------------------------------------------------------------- */
Real MaterialCohesive::getDissipatedEnergy() {
  AKANTU_DEBUG_IN();
  Real edis = 0.;

  /// integrate dissipated energy for each type of elements
  for (const auto & type : element_filter.elementTypes(
           spatial_dimension, _not_ghost, _ek_cohesive)) {
    Array<Real> dissipated_energy(total_energy(type, _not_ghost));
    dissipated_energy -= reversible_energy(type, _not_ghost);
    edis += fem_cohesive.integrate(dissipated_energy, type, _not_ghost,
                                   element_filter(type, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return edis;
}

/* -------------------------------------------------------------------------- */
Real MaterialCohesive::getContactEnergy() {
  AKANTU_DEBUG_IN();
  Real econ = 0.;

  /// integrate contact energy for each type of elements
  for (const auto & type : element_filter.elementTypes(
           spatial_dimension, _not_ghost, _ek_cohesive)) {

    auto & el_filter = element_filter(type, _not_ghost);
    UInt nb_quad_per_el = fem_cohesive.getNbIntegrationPoints(type, _not_ghost);
    UInt nb_quad_points = el_filter.size() * nb_quad_per_el;
    Array<Real> contact_energy(nb_quad_points);

    auto contact_traction_it =
        contact_tractions(type, _not_ghost).begin(spatial_dimension);
    auto contact_opening_it =
        contact_opening(type, _not_ghost).begin(spatial_dimension);

    /// loop on each quadrature point
    for (UInt q = 0; q < nb_quad_points;
         ++contact_traction_it, ++contact_opening_it, ++q) {

      contact_energy(q) = .5 * contact_traction_it->dot(*contact_opening_it);
    }

    econ += fem_cohesive.integrate(contact_energy, type, _not_ghost, el_filter);
  }

  AKANTU_DEBUG_OUT();
  return econ;
}

/* -------------------------------------------------------------------------- */
Real MaterialCohesive::getEnergy(const std::string & type) {
  if (type == "reversible") {
    return getReversibleEnergy();
  }
  if (type == "dissipated") {
    return getDissipatedEnergy();
  }
  if (type == "cohesive contact") {
    return getContactEnergy();
  }

  return 0.;
}

/* -------------------------------------------------------------------------- */
} // namespace akantu
