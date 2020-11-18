/**
 * @file   material_reinforcement_tmpl.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Wed Mar 25 2015
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Reinforcement material
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#include "aka_common.hh"
#include "aka_voigthelper.hh"
#include "material_reinforcement.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class Mat, UInt dim>
MaterialReinforcement<Mat, dim>::MaterialReinforcement(
    EmbeddedInterfaceModel & model, const ID & id)
    : Mat(model, 1, model.getInterfaceMesh(),
          model.getFEEngine("EmbeddedInterfaceFEEngine"), id),
      emodel(model),
      gradu_embedded("gradu_embedded", *this, 1,
                     model.getFEEngine("EmbeddedInterfaceFEEngine"),
                     this->element_filter),
      directing_cosines("directing_cosines", *this, 1,
                        model.getFEEngine("EmbeddedInterfaceFEEngine"),
                        this->element_filter),
      pre_stress("pre_stress", *this, 1,
                 model.getFEEngine("EmbeddedInterfaceFEEngine"),
                 this->element_filter),
      area(1.0) {
  AKANTU_DEBUG_IN();
  this->initialize();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Mat, UInt dim>
void MaterialReinforcement<Mat, dim>::initialize() {
  AKANTU_DEBUG_IN();

  this->registerParam("area", area, _pat_parsable | _pat_modifiable,
                      "Reinforcement cross-sectional area");
  this->registerParam("pre_stress", pre_stress, _pat_parsable | _pat_modifiable,
                      "Uniform pre-stress");

  // this->unregisterInternal(this->stress);

  // Fool the AvgHomogenizingFunctor
  // stress.initialize(dim * dim);

  // Reallocate the element filter
  this->element_filter.initialize(this->emodel.getInterfaceMesh(),
                                  _spatial_dimension = 1);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template <class Mat, UInt dim>
MaterialReinforcement<Mat, dim>::~MaterialReinforcement() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template <class Mat, UInt dim>
void MaterialReinforcement<Mat, dim>::initMaterial() {
  Mat::initMaterial();

  gradu_embedded.initialize(dim * dim);
  pre_stress.initialize(1);

  /// We initialise the stuff that is not going to change during the simulation
  this->initFilters();
  this->allocBackgroundShapeDerivatives();
  this->initBackgroundShapeDerivatives();
  this->initDirectingCosines();
}

/* -------------------------------------------------------------------------- */
/// Initialize the filter for background elements
template <class Mat, UInt dim>
void MaterialReinforcement<Mat, dim>::initFilters() {
  for (auto gt : ghost_types) {
    for (auto && type : emodel.getInterfaceMesh().elementTypes(1, gt)) {
      std::string shaped_id = "filter";
      if (gt == _ghost) {
        shaped_id += ":ghost";
      }
      auto & background =
          background_filter(std::make_unique<ElementTypeMapArray<UInt>>(
                                "bg_" + shaped_id, this->name),
                            type, gt);
      auto & foreground = foreground_filter(
          std::make_unique<ElementTypeMapArray<UInt>>(shaped_id, this->name),
          type, gt);
      foreground->initialize(emodel.getMesh(), _nb_component = 1,
                             _ghost_type = gt);
      background->initialize(emodel.getMesh(), _nb_component = 1,
                             _ghost_type = gt);

      // Computing filters
      for (auto && bg_type : background->elementTypes(dim, gt)) {
        filterInterfaceBackgroundElements(
            (*foreground)(bg_type), (*background)(bg_type), bg_type, type, gt);
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
/// Construct a filter for a (interface_type, background_type) pair
template <class Mat, UInt dim>
void MaterialReinforcement<Mat, dim>::filterInterfaceBackgroundElements(
    Array<UInt> & foreground, Array<UInt> & background, ElementType type,
    ElementType interface_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  foreground.resize(0);
  background.resize(0);

  Array<Element> & elements =
      emodel.getInterfaceAssociatedElements(interface_type, ghost_type);
  Array<UInt> & elem_filter = this->element_filter(interface_type, ghost_type);

  for (auto & elem_id : elem_filter) {
    Element & elem = elements(elem_id);
    if (elem.type == type) {
      background.push_back(elem.element);
      foreground.push_back(elem_id);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

namespace detail {
  class BackgroundShapeDInitializer : public ElementTypeMapArrayInitializer {
  public:
    BackgroundShapeDInitializer(UInt spatial_dimension, FEEngine & engine,
                                ElementType foreground_type,
                                const ElementTypeMapArray<UInt> & filter,
                                GhostType ghost_type)
        : ElementTypeMapArrayInitializer(
              [](ElementType bgtype, GhostType /*unused*/) {
                return ShapeFunctions::getShapeDerivativesSize(bgtype);
              },
              spatial_dimension, ghost_type, _ek_regular) {
      auto nb_quad = engine.getNbIntegrationPoints(foreground_type);
      // Counting how many background elements are affected by elements of
      // interface_type
      for (auto type : filter.elementTypes(this->spatial_dimension)) {
        // Inserting size
        array_size_per_bg_type(filter(type).size() * nb_quad, type,
                               this->ghost_type);
      }
    }

    auto elementTypes() const -> decltype(auto) {
      return array_size_per_bg_type.elementTypes();
    }

    UInt size(ElementType bgtype) const {
      return array_size_per_bg_type(bgtype, this->ghost_type);
    }

  protected:
    ElementTypeMap<UInt> array_size_per_bg_type;
  };
} // namespace detail

/**
 * Background shape derivatives need to be stored per background element
 * types but also per embedded element type, which is why they are stored
 * in an ElementTypeMap<ElementTypeMapArray<Real> *>. The outer ElementTypeMap
 * refers to the embedded types, and the inner refers to the background types.
 */
template <class Mat, UInt dim>
void MaterialReinforcement<Mat, dim>::allocBackgroundShapeDerivatives() {
  AKANTU_DEBUG_IN();

  for (auto gt : ghost_types) {
    for (auto && type : emodel.getInterfaceMesh().elementTypes(1, gt)) {
      std::string shaped_id = "embedded_shape_derivatives";
      if (gt == _ghost) {
        shaped_id += ":ghost";
      }

      auto & shaped_etma = shape_derivatives(
          std::make_unique<ElementTypeMapArray<Real>>(shaped_id, this->name),
          type, gt);
      shaped_etma->initialize(
          detail::BackgroundShapeDInitializer(
              emodel.getSpatialDimension(),
              emodel.getFEEngine("EmbeddedInterfaceFEEngine"), type,
              *background_filter(type, gt), gt),
          0, true);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template <class Mat, UInt dim>
void MaterialReinforcement<Mat, dim>::initBackgroundShapeDerivatives() {
  AKANTU_DEBUG_IN();

  for (auto interface_type :
       this->element_filter.elementTypes(this->spatial_dimension)) {
    for (auto type : background_filter(interface_type)->elementTypes(dim)) {
      computeBackgroundShapeDerivatives(interface_type, type, _not_ghost,
                                        this->element_filter(interface_type));
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Mat, UInt dim>
void MaterialReinforcement<Mat, dim>::computeBackgroundShapeDerivatives(
    ElementType interface_type, ElementType bg_type, GhostType ghost_type,
    const Array<UInt> & filter) {
  auto & interface_engine = emodel.getFEEngine("EmbeddedInterfaceFEEngine");
  auto & engine = emodel.getFEEngine();
  auto & interface_mesh = emodel.getInterfaceMesh();

  const auto nb_nodes_elem_bg = Mesh::getNbNodesPerElement(bg_type);
  // const auto nb_strss = VoigtHelper<dim>::size;
  const auto nb_quads_per_elem =
      interface_engine.getNbIntegrationPoints(interface_type);

  Array<Real> quad_pos(0, dim, "interface_quad_pos");
  interface_engine.interpolateOnIntegrationPoints(interface_mesh.getNodes(),
                                                  quad_pos, dim, interface_type,
                                                  ghost_type, filter);
  auto & background_shapesd =
      (*shape_derivatives(interface_type, ghost_type))(bg_type, ghost_type);
  auto & background_elements =
      (*background_filter(interface_type, ghost_type))(bg_type, ghost_type);
  auto & foreground_elements =
      (*foreground_filter(interface_type, ghost_type))(bg_type, ghost_type);

  auto shapesd_begin =
      background_shapesd.begin(dim, nb_nodes_elem_bg, nb_quads_per_elem);
  auto quad_begin = quad_pos.begin(dim, nb_quads_per_elem);

  for (auto && tuple : zip(background_elements, foreground_elements)) {
    auto bg = std::get<0>(tuple);
    auto fg = std::get<1>(tuple);
    for (UInt i = 0; i < nb_quads_per_elem; ++i) {
      Matrix<Real> shapesd = Tensor3<Real>(shapesd_begin[fg])(i);
      Vector<Real> quads = Matrix<Real>(quad_begin[fg])(i);

      engine.computeShapeDerivatives(quads, bg, bg_type, shapesd, ghost_type);
    }
  }
}

/* -------------------------------------------------------------------------- */

template <class Mat, UInt dim>
void MaterialReinforcement<Mat, dim>::initDirectingCosines() {
  AKANTU_DEBUG_IN();

  Mesh & mesh = emodel.getInterfaceMesh();

  const UInt voigt_size = VoigtHelper<dim>::size;
  directing_cosines.initialize(voigt_size);

  for (auto && type : mesh.elementTypes(1, _not_ghost)) {
    computeDirectingCosines(type, _not_ghost);
    // computeDirectingCosines(*type_it, _ghost);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template <class Mat, UInt dim>
void MaterialReinforcement<Mat, dim>::assembleStiffnessMatrix(
    GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Mesh & interface_mesh = emodel.getInterfaceMesh();

  for (auto && type : interface_mesh.elementTypes(1, _not_ghost)) {
    assembleStiffnessMatrix(type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template <class Mat, UInt dim>
void MaterialReinforcement<Mat, dim>::assembleInternalForces(
    GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Mesh & interface_mesh = emodel.getInterfaceMesh();

  for (auto && type : interface_mesh.elementTypes(1, _not_ghost)) {
    this->assembleInternalForces(type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Mat, UInt dim>
void MaterialReinforcement<Mat, dim>::computeAllStresses(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Mesh & interface_mesh = emodel.getInterfaceMesh();
  for (auto && type : interface_mesh.elementTypes(_ghost_type = ghost_type)) {
    computeGradU(type, ghost_type);
    this->computeStress(type, ghost_type);
    addPrestress(type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Mat, UInt dim>
void MaterialReinforcement<Mat, dim>::addPrestress(ElementType type,
                                                   GhostType ghost_type) {
  auto & stress = this->stress(type, ghost_type);
  auto & sigma_p = this->pre_stress(type, ghost_type);

  for (auto && tuple : zip(stress, sigma_p)) {
    std::get<0>(tuple) += std::get<1>(tuple);
  }
}

/* -------------------------------------------------------------------------- */
template <class Mat, UInt dim>
void MaterialReinforcement<Mat, dim>::assembleInternalForces(
    ElementType type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Mesh & mesh = emodel.getMesh();

  for (auto && mesh_type : mesh.elementTypes(dim, ghost_type)) {
    assembleInternalForcesInterface(type, mesh_type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

/**
 * Computes and assemble the residual. Residual in reinforcement is computed as:
 *
 * \f[
 * \vec{r} = A_s \int_S{\mathbf{B}^T\mathbf{C}^T \vec{\sigma_s}\,\mathrm{d}s}
 * \f]
 */
template <class Mat, UInt dim>
void MaterialReinforcement<Mat, dim>::assembleInternalForcesInterface(
    ElementType interface_type, ElementType background_type,
    GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  UInt voigt_size = VoigtHelper<dim>::size;

  FEEngine & interface_engine = emodel.getFEEngine("EmbeddedInterfaceFEEngine");

  Array<UInt> & elem_filter = this->element_filter(interface_type, ghost_type);

  UInt nodes_per_background_e = Mesh::getNbNodesPerElement(background_type);
  UInt nb_quadrature_points =
      interface_engine.getNbIntegrationPoints(interface_type, ghost_type);
  UInt nb_element = elem_filter.size();

  UInt back_dof = dim * nodes_per_background_e;

  Array<Real> & shapesd = (*shape_derivatives(interface_type, ghost_type))(
      background_type, ghost_type);

  Array<Real> integrant(nb_quadrature_points * nb_element, back_dof,
                        "integrant");

  auto integrant_it = integrant.begin(back_dof);
  auto integrant_end = integrant.end(back_dof);

  Array<Real>::matrix_iterator B_it =
      shapesd.begin(dim, nodes_per_background_e);
  auto C_it = directing_cosines(interface_type, ghost_type).begin(voigt_size);

  auto sigma_it = this->stress(interface_type, ghost_type).begin();

  Matrix<Real> Bvoigt(voigt_size, back_dof);

  for (; integrant_it != integrant_end;
       ++integrant_it, ++B_it, ++C_it, ++sigma_it) {
    VoigtHelper<dim>::transferBMatrixToSymVoigtBMatrix(*B_it, Bvoigt,
                                                       nodes_per_background_e);
    Vector<Real> & C = *C_it;
    Vector<Real> & BtCt_sigma = *integrant_it;

    BtCt_sigma.mul<true>(Bvoigt, C);
    BtCt_sigma *= *sigma_it * area;
  }

  Array<Real> residual_interface(nb_element, back_dof, "residual_interface");
  interface_engine.integrate(integrant, residual_interface, back_dof,
                             interface_type, ghost_type, elem_filter);
  integrant.resize(0);

  Array<UInt> background_filter(nb_element, 1, "background_filter");

  auto & filter =
      getBackgroundFilter(interface_type, background_type, ghost_type);

  emodel.getDOFManager().assembleElementalArrayLocalArray(
      residual_interface, emodel.getInternalForce(), background_type,
      ghost_type, -1., filter);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template <class Mat, UInt dim>
void MaterialReinforcement<Mat, dim>::computeDirectingCosines(
    ElementType type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Mesh & interface_mesh = emodel.getInterfaceMesh();

  const UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  const UInt steel_dof = dim * nb_nodes_per_element;
  const UInt voigt_size = VoigtHelper<dim>::size;
  const UInt nb_quad_points = emodel.getFEEngine("EmbeddedInterfaceFEEngine")
                                  .getNbIntegrationPoints(type, ghost_type);

  Array<Real> node_coordinates(this->element_filter(type, ghost_type).size(),
                               steel_dof);

  this->emodel.getFEEngine().template extractNodalToElementField<Real>(
      interface_mesh, interface_mesh.getNodes(), node_coordinates, type,
      ghost_type, this->element_filter(type, ghost_type));

  Array<Real>::matrix_iterator directing_cosines_it =
      directing_cosines(type, ghost_type).begin(1, voigt_size);

  Array<Real>::matrix_iterator node_coordinates_it =
      node_coordinates.begin(dim, nb_nodes_per_element);
  Array<Real>::matrix_iterator node_coordinates_end =
      node_coordinates.end(dim, nb_nodes_per_element);

  for (; node_coordinates_it != node_coordinates_end; ++node_coordinates_it) {
    for (UInt i = 0; i < nb_quad_points; i++, ++directing_cosines_it) {
      Matrix<Real> & nodes = *node_coordinates_it;
      Matrix<Real> & cosines = *directing_cosines_it;

      computeDirectingCosinesOnQuad(nodes, cosines);
    }
  }

  // Mauro: the directing_cosines internal is defined on the quadrature points
  // of each element

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template <class Mat, UInt dim>
void MaterialReinforcement<Mat, dim>::assembleStiffnessMatrix(
    ElementType type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Mesh & mesh = emodel.getMesh();

  for (auto && mesh_type : mesh.elementTypes(dim, ghost_type)) {
    assembleStiffnessMatrixInterface(type, mesh_type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Computes the reinforcement stiffness matrix (Gomes & Awruch, 2001)
 * \f[
 * \mathbf{K}_e = \sum_{i=1}^R{A_i\int_{S_i}{\mathbf{B}^T
 * \mathbf{C}_i^T \mathbf{D}_{s, i} \mathbf{C}_i \mathbf{B}\,\mathrm{d}s}}
 * \f]
 */
template <class Mat, UInt dim>
void MaterialReinforcement<Mat, dim>::assembleStiffnessMatrixInterface(
    ElementType interface_type, ElementType background_type,
    GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  UInt voigt_size = VoigtHelper<dim>::size;

  FEEngine & interface_engine = emodel.getFEEngine("EmbeddedInterfaceFEEngine");

  Array<UInt> & elem_filter = this->element_filter(interface_type, ghost_type);
  Array<Real> & grad_u = gradu_embedded(interface_type, ghost_type);

  UInt nb_element = elem_filter.size();
  UInt nodes_per_background_e = Mesh::getNbNodesPerElement(background_type);
  UInt nb_quadrature_points =
      interface_engine.getNbIntegrationPoints(interface_type, ghost_type);

  UInt back_dof = dim * nodes_per_background_e;

  UInt integrant_size = back_dof;

  grad_u.resize(nb_quadrature_points * nb_element);

  Array<Real> tangent_moduli(nb_element * nb_quadrature_points, 1,
                             "interface_tangent_moduli");
  this->computeTangentModuli(interface_type, tangent_moduli, ghost_type);

  Array<Real> & shapesd = (*shape_derivatives(interface_type, ghost_type))(
      background_type, ghost_type);

  Array<Real> integrant(nb_element * nb_quadrature_points,
                        integrant_size * integrant_size, "B^t*C^t*D*C*B");

  /// Temporary matrices for integrant product
  Matrix<Real> Bvoigt(voigt_size, back_dof);
  Matrix<Real> DCB(1, back_dof);
  Matrix<Real> CtDCB(voigt_size, back_dof);

  Array<Real>::scalar_iterator D_it = tangent_moduli.begin();
  Array<Real>::scalar_iterator D_end = tangent_moduli.end();

  Array<Real>::matrix_iterator C_it =
      directing_cosines(interface_type, ghost_type).begin(1, voigt_size);
  Array<Real>::matrix_iterator B_it =
      shapesd.begin(dim, nodes_per_background_e);
  Array<Real>::matrix_iterator integrant_it =
      integrant.begin(integrant_size, integrant_size);

  for (; D_it != D_end; ++D_it, ++C_it, ++B_it, ++integrant_it) {
    Real & D = *D_it;
    Matrix<Real> & C = *C_it;
    Matrix<Real> & B = *B_it;
    Matrix<Real> & BtCtDCB = *integrant_it;

    VoigtHelper<dim>::transferBMatrixToSymVoigtBMatrix(B, Bvoigt,
                                                       nodes_per_background_e);

    DCB.mul<false, false>(C, Bvoigt);
    DCB *= D * area;
    CtDCB.mul<true, false>(C, DCB);
    BtCtDCB.mul<true, false>(Bvoigt, CtDCB);
  }

  tangent_moduli.resize(0);

  Array<Real> K_interface(nb_element, integrant_size * integrant_size,
                          "K_interface");
  interface_engine.integrate(integrant, K_interface,
                             integrant_size * integrant_size, interface_type,
                             ghost_type, elem_filter);

  integrant.resize(0);

  // Mauro: Here K_interface contains the local stiffness matrices,
  // directing_cosines contains the information about the orientation
  // of the reinforcements, any rotation of the local stiffness matrix
  // can be done here

  auto & filter =
      getBackgroundFilter(interface_type, background_type, ghost_type);

  emodel.getDOFManager().assembleElementalMatricesToMatrix(
      "K", "displacement", K_interface, background_type, ghost_type, _symmetric,
      filter);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Mat, UInt dim>
Real MaterialReinforcement<Mat, dim>::getEnergy(const std::string & id) {
  AKANTU_DEBUG_IN();
  if (id == "potential") {
    Real epot = 0.;

    this->computePotentialEnergyByElements();

    for (auto && type :
         this->element_filter.elementTypes(this->spatial_dimension)) {
      FEEngine & interface_engine =
          emodel.getFEEngine("EmbeddedInterfaceFEEngine");
      epot += interface_engine.integrate(
          this->potential_energy(type, _not_ghost), type, _not_ghost,
          this->element_filter(type, _not_ghost));
      epot *= area;
    }

    return epot;
  }

  AKANTU_DEBUG_OUT();
  return 0;
}

/* -------------------------------------------------------------------------- */

template <class Mat, UInt dim>
void MaterialReinforcement<Mat, dim>::computeGradU(ElementType interface_type,
                                                   GhostType ghost_type) {
  // Looping over background types
  for (auto && bg_type :
       background_filter(interface_type, ghost_type)->elementTypes(dim)) {

    const UInt nodes_per_background_e = Mesh::getNbNodesPerElement(bg_type);
    const UInt voigt_size = VoigtHelper<dim>::size;

    auto & bg_shapesd =
        (*shape_derivatives(interface_type, ghost_type))(bg_type, ghost_type);

    auto & filter = getBackgroundFilter(interface_type, bg_type, ghost_type);

    Array<Real> disp_per_element(0, dim * nodes_per_background_e, "disp_elem");

    FEEngine::extractNodalToElementField(
        emodel.getMesh(), emodel.getDisplacement(), disp_per_element, bg_type,
        ghost_type, filter);

    Matrix<Real> concrete_du(dim, dim);
    Matrix<Real> epsilon(dim, dim);
    Vector<Real> evoigt(voigt_size);

    for (auto && tuple :
         zip(make_view(disp_per_element, dim, nodes_per_background_e),
             make_view(bg_shapesd, dim, nodes_per_background_e),
             this->gradu(interface_type, ghost_type),
             make_view(this->directing_cosines(interface_type, ghost_type),
                       voigt_size))) {
      auto & u = std::get<0>(tuple);
      auto & B = std::get<1>(tuple);
      auto & du = std::get<2>(tuple);
      auto & C = std::get<3>(tuple);

      concrete_du.mul<false, true>(u, B);
      auto epsilon = 0.5 * (concrete_du + concrete_du.transpose());
      strainTensorToVoigtVector(epsilon, evoigt);
      du = C.dot(evoigt);
    }
  }
}

/* -------------------------------------------------------------------------- */

/**
 * The structure of the directing cosines matrix is :
 * \f{eqnarray*}{
 *  C_{1,\cdot} & = & (l^2, m^2, n^2, mn, ln, lm) \\
 *  C_{i,j} & = & 0
 * \f}
 *
 * with :
 * \f[
 * (l, m, n) = \frac{1}{\|\frac{\mathrm{d}\vec{r}(s)}{\mathrm{d}s}\|} \cdot
 * \frac{\mathrm{d}\vec{r}(s)}{\mathrm{d}s}
 * \f]
 */
template <class Mat, UInt dim>
inline void MaterialReinforcement<Mat, dim>::computeDirectingCosinesOnQuad(
    const Matrix<Real> & nodes, Matrix<Real> & cosines) {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(nodes.cols() == 2,
                      "Higher order reinforcement elements not implemented");

  const Vector<Real> a = nodes(0);
  const Vector<Real> b = nodes(1);
  Vector<Real> delta = b - a;

  Real sq_length = 0.;
  for (UInt i = 0; i < dim; i++) {
    sq_length += delta(i) * delta(i);
  }

  if (dim == 2) {
    cosines(0, 0) = delta(0) * delta(0); // l^2
    cosines(0, 1) = delta(1) * delta(1); // m^2
    cosines(0, 2) = delta(0) * delta(1); // lm
  } else if (dim == 3) {
    cosines(0, 0) = delta(0) * delta(0); // l^2
    cosines(0, 1) = delta(1) * delta(1); // m^2
    cosines(0, 2) = delta(2) * delta(2); // n^2

    cosines(0, 3) = delta(1) * delta(2); // mn
    cosines(0, 4) = delta(0) * delta(2); // ln
    cosines(0, 5) = delta(0) * delta(1); // lm
  }

  cosines /= sq_length;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template <class Mat, UInt dim>
inline void MaterialReinforcement<Mat, dim>::stressTensorToVoigtVector(
    const Matrix<Real> & tensor, Vector<Real> & vector) {
  AKANTU_DEBUG_IN();

  for (UInt i = 0; i < dim; i++) {
    vector(i) = tensor(i, i);
  }

  if (dim == 2) {
    vector(2) = tensor(0, 1);
  } else if (dim == 3) {
    vector(3) = tensor(1, 2);
    vector(4) = tensor(0, 2);
    vector(5) = tensor(0, 1);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template <class Mat, UInt dim>
inline void MaterialReinforcement<Mat, dim>::strainTensorToVoigtVector(
    const Matrix<Real> & tensor, Vector<Real> & vector) {
  AKANTU_DEBUG_IN();

  for (UInt i = 0; i < dim; i++) {
    vector(i) = tensor(i, i);
  }

  if (dim == 2) {
    vector(2) = 2 * tensor(0, 1);
  } else if (dim == 3) {
    vector(3) = 2 * tensor(1, 2);
    vector(4) = 2 * tensor(0, 2);
    vector(5) = 2 * tensor(0, 1);
  }

  AKANTU_DEBUG_OUT();
}

} // namespace akantu
