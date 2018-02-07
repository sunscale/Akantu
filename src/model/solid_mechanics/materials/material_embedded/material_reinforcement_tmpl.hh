/**
 * @file   material_reinforcement_tmpl.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Thu Feb 1 2018
 *
 * @brief  Reinforcement material
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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

#include "aka_common.hh"
#include "aka_voigthelper.hh"
#include "material_reinforcement.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class Mat, UInt dim>
MaterialReinforcement<Mat, dim>::MaterialReinforcement(
    EmbeddedInterfaceModel & model, const ID & id)
    : Mat(model, 1,
          model.getInterfaceMesh(),
          model.getFEEngine("EmbeddedInterfaceFEEngine"), id),
      emodel(model),
      stress_embedded("stress_embedded", *this, 1,
                      model.getFEEngine("EmbeddedInterfaceFEEngine"),
                      this->element_filter),
      gradu_embedded("gradu_embedded", *this, 1,
                     model.getFEEngine("EmbeddedInterfaceFEEngine"),
                     this->element_filter),
      directing_cosines("directing_cosines", *this, 1,
                        model.getFEEngine("EmbeddedInterfaceFEEngine"),
                        this->element_filter),
      pre_stress("pre_stress", *this, 1,
                 model.getFEEngine("EmbeddedInterfaceFEEngine"),
                 this->element_filter),
      area(1.0), shape_derivatives() {
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

  this->unregisterInternal(this->stress);

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

  stress_embedded.initialize(dim * dim);
  gradu_embedded.initialize(dim * dim);
  pre_stress.initialize(1);

  /// We initialise the stuff that is not going to change during the simulation
  this->initFilters();
  this->allocBackgroundShapeDerivatives();
  this->initBackgroundShapeDerivatives();
  this->initDirectingCosines();
}

/* -------------------------------------------------------------------------- */

namespace detail {
  class FilterInitializer : public MeshElementTypeMapArrayInializer {
  public:
    FilterInitializer(EmbeddedInterfaceModel & emodel,
                      const GhostType & ghost_type)
        : MeshElementTypeMapArrayInializer(emodel.getMesh(),
                                           1, emodel.getSpatialDimension(),
                                           ghost_type, _ek_regular) {}

    UInt size(const ElementType & /*bgtype*/) const override { return 0; }
  };
}

/* -------------------------------------------------------------------------- */
/// Initialize the filter for background elements
template <class Mat, UInt dim>
void MaterialReinforcement<Mat, dim>::initFilters() {
  for (auto gt : ghost_types) {
    for (auto && type : emodel.getInterfaceMesh().elementTypes(1, gt)) {
      std::string shaped_id = "filter";
      if (gt == _ghost)
        shaped_id += ":ghost";

      auto & background =
          background_filter(std::make_unique<ElementTypeMapArray<UInt>>(
                                "bg_" + shaped_id, this->name),
                            type, gt);
      auto & foreground = foreground_filter(
          std::make_unique<ElementTypeMapArray<UInt>>(shaped_id, this->name),
          type, gt);
      foreground->initialize(detail::FilterInitializer(emodel, gt), 0, true);
      background->initialize(detail::FilterInitializer(emodel, gt), 0, true);

      // Computing filters
      for (auto && bg_type : background->elementTypes(dim, gt)) {
	std::cout << "Type " << bg_type << std::endl;
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
    Array<UInt> & foreground, Array<UInt> & background,
    const ElementType & type, const ElementType & interface_type,
    GhostType ghost_type) {
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
  class BackgroundShapeDInitializer : public ElementTypeMapArrayInializer {
  public:
    BackgroundShapeDInitializer(UInt spatial_dimension,
				FEEngine & engine,
				const ElementType & foreground_type,
                                ElementTypeMapArray<UInt> & filter,
                                const GhostType & ghost_type)
        : ElementTypeMapArrayInializer(spatial_dimension, 0, ghost_type,
                                       _ek_regular) {
      auto nb_quad = engine.getNbIntegrationPoints(foreground_type);
      // Counting how many background elements are affected by elements of
      // interface_type
      for (auto type : filter.elementTypes(this->spatial_dimension)) {
        // Inserting size
        array_size_per_bg_type(filter(type).size() * nb_quad, type,
                               this->ghost_type);
      }

      std::cout << array_size_per_bg_type << std::endl;
      std::cout << this->spatial_dimension << std::endl;
    }

    auto elementTypes() const -> decltype(auto) {
      return array_size_per_bg_type.elementTypes();
    }

    UInt size(const ElementType & bgtype) const {
      return array_size_per_bg_type(bgtype, this->ghost_type);
    }

    UInt nbComponent(const ElementType & bgtype) const {
      return ShapeFunctions::getShapeDerivativesSize(bgtype);
    }

  protected:
    ElementTypeMap<UInt> array_size_per_bg_type;
  };
}

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
      if (gt == _ghost)
        shaped_id += ":ghost";

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
    const ElementType & interface_type, const ElementType & bg_type,
    GhostType ghost_type, const Array<UInt> & filter) {
  auto & interface_engine = emodel.getFEEngine("EmbeddedInterfaceFEEngine");
  auto & engine = emodel.getFEEngine();
  auto & interface_mesh = emodel.getInterfaceMesh();

  const auto ndof = dim * Mesh::getNbNodesPerElement(bg_type);
  const auto nb_stress =
      ShapeFunctions::getShapeDerivativesSize(bg_type) / ndof;
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
      background_shapesd.begin(nb_stress, ndof, nb_quads_per_elem);
  auto quad_begin = quad_pos.begin(dim, Mesh::getNbNodesPerElement(bg_type));

  for (auto && tuple : zip(background_elements, foreground_elements)) {
    UInt bg = std::get<0>(tuple), fg = std::get<1>(tuple);
    for (UInt i = 0; i < nb_quads_per_elem; ++i) {
      Matrix<Real> shapesd = Tensor3<Real>(shapesd_begin[fg])(i);
      Vector<Real> quads = Matrix<Real>(quad_begin[fg])(i);

      engine.computeShapeDerivatives(quads, bg, bg_type, shapesd,
                                     ghost_type);
    }
  }
}

/* -------------------------------------------------------------------------- */

template <class Mat, UInt dim>
void MaterialReinforcement<Mat, dim>::initDirectingCosines() {
  AKANTU_DEBUG_IN();

  Mesh & mesh = emodel.getInterfaceMesh();

  Mesh::type_iterator type_it = mesh.firstType(1, _not_ghost);
  Mesh::type_iterator type_end = mesh.lastType(1, _not_ghost);

  const UInt voigt_size = VoigtHelper<dim>::size;
  directing_cosines.initialize(voigt_size * voigt_size);

  for (; type_it != type_end; ++type_it) {
    computeDirectingCosines(*type_it, _not_ghost);
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

  Mesh::type_iterator type_it = interface_mesh.firstType(1, _not_ghost);
  Mesh::type_iterator type_end = interface_mesh.lastType(1, _not_ghost);

  for (; type_it != type_end; ++type_it) {
    assembleStiffnessMatrix(*type_it, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template <class Mat, UInt dim>
void MaterialReinforcement<Mat, dim>::assembleInternalForces(
    GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Mesh & interface_mesh = emodel.getInterfaceMesh();

  Mesh::type_iterator type_it = interface_mesh.firstType(1, _not_ghost);
  Mesh::type_iterator type_end = interface_mesh.lastType(1, _not_ghost);

  for (; type_it != type_end; ++type_it) {
    this->assembleInternalForces(*type_it, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Mat, UInt dim>
void MaterialReinforcement<Mat, dim>::computeGradU(const ElementType & type,
                                                   GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  UInt nb_quad_points = emodel.getFEEngine("EmbeddedInterfaceFEEngine")
                            .getNbIntegrationPoints(type);

  Array<Real> & gradu_vec = gradu_embedded(type, ghost_type);

  Mesh::type_iterator back_it = emodel.getMesh().firstType(dim, ghost_type);
  Mesh::type_iterator back_end = emodel.getMesh().lastType(dim, ghost_type);

  for (; back_it != back_end; ++back_it) {
    UInt nodes_per_background_e = Mesh::getNbNodesPerElement(*back_it);

    Array<Real> & shapesd =
        shape_derivatives(type, ghost_type)->operator()(*back_it, ghost_type);

    auto & filter = getBackgroundFilter(type, *back_it, ghost_type);

    Array<Real> disp_per_element(0, dim * nodes_per_background_e, "disp_elem");

    FEEngine::extractNodalToElementField(
        emodel.getMesh(), emodel.getDisplacement(), disp_per_element, *back_it,
        ghost_type, filter);

    Array<Real>::matrix_iterator disp_it =
        disp_per_element.begin(dim, nodes_per_background_e);
    Array<Real>::matrix_iterator disp_end =
        disp_per_element.end(dim, nodes_per_background_e);

    Array<Real>::matrix_iterator shapes_it =
        shapesd.begin(dim, nodes_per_background_e);
    Array<Real>::matrix_iterator grad_u_it = gradu_vec.begin(dim, dim);

    for (; disp_it != disp_end; ++disp_it) {
      for (UInt i = 0; i < nb_quad_points; i++, ++shapes_it, ++grad_u_it) {
        Matrix<Real> & B = *shapes_it;
        Matrix<Real> & du = *grad_u_it;
        Matrix<Real> & u = *disp_it;

        du.mul<false, true>(u, B);
      }
    }

  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Mat, UInt dim>
void MaterialReinforcement<Mat, dim>::computeAllStresses(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Mesh::type_iterator it = emodel.getInterfaceMesh().firstType();
  Mesh::type_iterator last_type = emodel.getInterfaceMesh().lastType();

  for (; it != last_type; ++it) {
    computeGradU(*it, ghost_type);
    this->computeStress(*it, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Mat, UInt dim>
void MaterialReinforcement<Mat, dim>::assembleInternalForces(
    const ElementType & type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Mesh & mesh = emodel.getMesh();

  Mesh::type_iterator type_it = mesh.firstType(dim, ghost_type);
  Mesh::type_iterator type_end = mesh.lastType(dim, ghost_type);

  for (; type_it != type_end; ++type_it) {
    assembleInternalForcesInterface(type, *type_it, ghost_type);
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
    const ElementType & interface_type, const ElementType & background_type,
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

  Array<Real>::vector_iterator integrant_it = integrant.begin(back_dof);
  Array<Real>::vector_iterator integrant_end = integrant.end(back_dof);

  Array<Real>::matrix_iterator B_it =
      shapesd.begin(dim, nodes_per_background_e);
  Array<Real>::matrix_iterator C_it =
      directing_cosines(interface_type, ghost_type)
          .begin(voigt_size, voigt_size);
  Array<Real>::matrix_iterator sigma_it =
      stress_embedded(interface_type, ghost_type).begin(dim, dim);

  Vector<Real> sigma(voigt_size);
  Matrix<Real> Bvoigt(voigt_size, back_dof);
  Vector<Real> Ct_sigma(voigt_size);

  for (; integrant_it != integrant_end;
       ++integrant_it, ++B_it, ++C_it, ++sigma_it) {
    VoigtHelper<dim>::transferBMatrixToSymVoigtBMatrix(*B_it, Bvoigt,
                                                       nodes_per_background_e);
    Matrix<Real> & C = *C_it;
    Vector<Real> & BtCt_sigma = *integrant_it;

    stressTensorToVoigtVector(*sigma_it, sigma);

    Ct_sigma.mul<true>(C, sigma);
    BtCt_sigma.mul<true>(Bvoigt, Ct_sigma);
    BtCt_sigma *= area;
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
    const ElementType & type, GhostType ghost_type) {
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
      directing_cosines(type, ghost_type).begin(voigt_size, voigt_size);

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
    const ElementType & type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Mesh & mesh = emodel.getMesh();

  Mesh::type_iterator type_it = mesh.firstType(dim, ghost_type);
  Mesh::type_iterator type_end = mesh.lastType(dim, ghost_type);

  for (; type_it != type_end; ++type_it) {
    assembleStiffnessMatrixInterface(type, *type_it, ghost_type);
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
    const ElementType & interface_type, const ElementType & background_type,
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
  Matrix<Real> DC(voigt_size, voigt_size);
  Matrix<Real> DCB(voigt_size, back_dof);
  Matrix<Real> CtDCB(voigt_size, back_dof);

  Array<Real>::scalar_iterator D_it = tangent_moduli.begin();
  Array<Real>::scalar_iterator D_end = tangent_moduli.end();

  Array<Real>::matrix_iterator C_it =
      directing_cosines(interface_type, ghost_type)
          .begin(voigt_size, voigt_size);
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

    DC.clear();
    DC(0, 0) = D * area;
    DC *= C;
    DCB.mul<false, false>(DC, Bvoigt);
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

    Mesh::type_iterator it = this->element_filter.firstType(
                            this->spatial_dimension),
                        end = this->element_filter.lastType(
                            this->spatial_dimension);

    for (; it != end; ++it) {
      FEEngine & interface_engine =
          emodel.getFEEngine("EmbeddedInterfaceFEEngine");
      epot += interface_engine.integrate(
          this->potential_energy(*it, _not_ghost), *it, _not_ghost,
          this->element_filter(*it, _not_ghost));
      epot *= area;
    }

    return epot;
  }

  AKANTU_DEBUG_OUT();
  return 0;
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

  const Vector<Real> a = nodes(0), b = nodes(1);
  Vector<Real> delta = b - a;

  cosines.clear();

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

} // akantu
