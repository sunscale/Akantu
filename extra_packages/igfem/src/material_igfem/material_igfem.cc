/**
 * @file   element_class_igfem.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 *
 * @brief  Implementation parent material for IGFEM
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

#include "material_igfem.hh"
#include "aka_math.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
MaterialIGFEM::MaterialIGFEM(SolidMechanicsModel & model, const ID & id)
    : Material(model, id), nb_sub_materials(2),
      sub_material("sub_material", *this), name_sub_mat_1(""),
      name_sub_mat_2("") {
  AKANTU_DEBUG_IN();

  this->model = dynamic_cast<SolidMechanicsModelIGFEM *>(&model);
  this->fem = &(model.getFEEngineClass<MyFEEngineIGFEMType>("IGFEMFEEngine"));

  this->model->getMesh().initElementTypeMapArray(
      element_filter, 1, spatial_dimension, false, _ek_igfem);

  this->initialize();

  AKANTU_DEBUG_OUT();
};

/* -------------------------------------------------------------------------- */

MaterialIGFEM::~MaterialIGFEM() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialIGFEM::initialize() {

  this->gradu.setElementKind(_ek_igfem);
  this->stress.setElementKind(_ek_igfem);
  this->eigengradu.setElementKind(_ek_igfem);

  this->gradu.setFEEngine(*fem);
  this->stress.setFEEngine(*fem);
  this->eigengradu.setFEEngine(*fem);

  registerParam("name_sub_mat_1", name_sub_mat_1, std::string(),
                _pat_parsable | _pat_readable);
  registerParam("name_sub_mat_2", name_sub_mat_2, std::string(),
                _pat_parsable | _pat_readable);

  this->sub_material.initialize(1);
}

/* -------------------------------------------------------------------------- */
void MaterialIGFEM::computeQuadraturePointsCoordinates(
    ElementTypeMapArray<Real> & quadrature_points_coordinates,
    const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  /// compute quadrature points position in undeformed configuration
  Array<Real> & nodes_coordinates = this->fem->getMesh().getNodes();
  Mesh::type_iterator it =
      this->element_filter.firstType(spatial_dimension, ghost_type, _ek_igfem);
  Mesh::type_iterator last_type =
      this->element_filter.lastType(spatial_dimension, ghost_type, _ek_igfem);
  for (; it != last_type; ++it) {
    const Array<UInt> & elem_filter = this->element_filter(*it, ghost_type);
    UInt nb_element = elem_filter.getSize();
    if (nb_element) {
      UInt nb_tot_quad =
          this->fem->getNbIntegrationPoints(*it, ghost_type) * nb_element;

      Array<Real> & quads = quadrature_points_coordinates(*it, ghost_type);
      quads.resize(nb_tot_quad);

      this->model->getFEEngine("IGFEMFEEngine")
          .interpolateOnIntegrationPoints(nodes_coordinates, quads,
                                          spatial_dimension, *it, ghost_type,
                                          elem_filter);
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
void MaterialIGFEM::computeAllStresses(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = model->getSpatialDimension();

  Mesh::type_iterator it =
      this->fem->getMesh().firstType(spatial_dimension, ghost_type, _ek_igfem);
  Mesh::type_iterator last_type =
      this->fem->getMesh().lastType(spatial_dimension, ghost_type, _ek_igfem);

  for (; it != last_type; ++it) {
    Array<UInt> & elem_filter = element_filter(*it, ghost_type);
    if (elem_filter.getSize()) {
      Array<Real> & gradu_vect = gradu(*it, ghost_type);

      /// compute @f$\nabla u@f$
      this->fem->gradientOnIntegrationPoints(model->getDisplacement(),
                                             gradu_vect, spatial_dimension, *it,
                                             ghost_type, elem_filter);

      gradu_vect -= eigengradu(*it, ghost_type);

      /// compute @f$\mathbf{\sigma}_q@f$ from @f$\nabla u@f$
      computeStress(*it, ghost_type);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/// extrapolate internal values
void MaterialIGFEM::extrapolateInternal(const ID & id, const Element & element,
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
      if (tmp.norm() > Math::getTolerance()) {
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
template <ElementType type>
void MaterialIGFEM::setSubMaterial(const Array<Element> & element_list,
                                   GhostType ghost_type) {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template <>
void MaterialIGFEM::setSubMaterial<_igfem_triangle_5>(
    const Array<Element> & element_list, GhostType ghost_type) {
  SolidMechanicsModelIGFEM * igfem_model =
      static_cast<SolidMechanicsModelIGFEM *>(this->model);
  Vector<UInt> sub_material_index(this->nb_sub_materials);
  Array<Element>::const_iterator<Element> el_begin = element_list.begin();
  Array<Element>::const_iterator<Element> el_end = element_list.end();
  const Mesh & mesh = this->model->getMesh();
  Array<Real> nodes_coordinates(mesh.getNodes(), true);
  Array<Real>::const_vector_iterator nodes_it =
      nodes_coordinates.begin(spatial_dimension);
  Element el;
  el.kind = _ek_igfem;
  el.type = _igfem_triangle_5;
  el.ghost_type = ghost_type;
  UInt nb_nodes_per_el = mesh.getNbNodesPerElement(el.type);
  UInt nb_parent_nodes = IGFEMHelper::getNbParentNodes(el.type);
  Vector<bool> is_inside(nb_parent_nodes);
  const Array<UInt> & connectivity = mesh.getConnectivity(el.type, ghost_type);
  Array<UInt>::const_vector_iterator connec_it =
      connectivity.begin(nb_nodes_per_el);

  /// get the number of quadrature points for the two sub-elements
  UInt quads_1 = IGFEMHelper::getNbQuadraturePoints(el.type, 0);
  UInt quads_2 = IGFEMHelper::getNbQuadraturePoints(el.type, 1);
  UInt nb_total_quads = quads_1 + quads_2;

  UInt * sub_mat_ptr = this->sub_material(el.type, ghost_type).storage();

  /// loop all elements for the given type
  const Array<UInt> & filter = this->element_filter(el.type, ghost_type);
  UInt nb_elements = filter.getSize();
  for (UInt e = 0; e < nb_elements; ++e, ++connec_it) {
    el.element = filter(e);
    if (std::find(el_begin, el_end, el) == el_end) {
      sub_mat_ptr += nb_total_quads;
      continue;
    }

    for (UInt i = 0; i < nb_parent_nodes; ++i) {
      Vector<Real> node = nodes_it[(*connec_it)(i)];
      is_inside(i) = igfem_model->isInside(node, this->name_sub_mat_1);
    }

    UInt orientation = IGFEMHelper::getElementOrientation(el.type, is_inside);

    switch (orientation) {
    case 0: {
      sub_material_index(0) = 0;
      sub_material_index(1) = 1;
      break;
    }

    case 1: {
      sub_material_index(0) = 1;
      sub_material_index(1) = 0;
      break;
    }

    case 2: {
      sub_material_index(0) = 0;
      sub_material_index(1) = 0;
      break;
    }
    case 3: {
      sub_material_index(0) = 1;
      sub_material_index(0) = 1;
      break;
    }
    }

    for (UInt q = 0; q < quads_1; ++q, ++sub_mat_ptr) {
      UInt index = sub_material_index(0);
      *sub_mat_ptr = index;
    }
    for (UInt q = 0; q < quads_2; ++q, ++sub_mat_ptr) {
      UInt index = sub_material_index(1);
      *sub_mat_ptr = index;
    }
  }
}

/* -------------------------------------------------------------------------- */
template <>
void MaterialIGFEM::setSubMaterial<_igfem_triangle_4>(
    const Array<Element> & element_list, GhostType ghost_type) {
  SolidMechanicsModelIGFEM * igfem_model =
      static_cast<SolidMechanicsModelIGFEM *>(this->model);
  Vector<UInt> sub_material_index(this->nb_sub_materials);
  Array<Element>::const_iterator<Element> el_begin = element_list.begin();
  Array<Element>::const_iterator<Element> el_end = element_list.end();
  const Mesh & mesh = this->model->getMesh();
  Element el;
  el.kind = _ek_igfem;
  el.ghost_type = ghost_type;
  el.type = _igfem_triangle_4;
  UInt nb_nodes_per_el = mesh.getNbNodesPerElement(el.type);
  Vector<Real> barycenter(spatial_dimension);
  const Array<UInt> & connectivity = mesh.getConnectivity(el.type, ghost_type);
  Array<UInt>::const_vector_iterator connec_it =
      connectivity.begin(nb_nodes_per_el);

  /// get the number of quadrature points for the two sub-elements
  UInt quads_1 = IGFEMHelper::getNbQuadraturePoints(el.type, 0);
  UInt quads_2 = IGFEMHelper::getNbQuadraturePoints(el.type, 1);
  UInt nb_total_quads = quads_1 + quads_2;

  UInt * sub_mat_ptr = this->sub_material(el.type, ghost_type).storage();

  /// loop all elements for the given type
  const Array<UInt> & filter = this->element_filter(el.type, ghost_type);
  UInt nb_elements = filter.getSize();
  for (UInt e = 0; e < nb_elements; ++e, ++connec_it) {
    el.element = filter(e);
    if (std::find(el_begin, el_end, el) == el_end) {
      sub_mat_ptr += nb_total_quads;
      continue;
    }

    for (UInt s = 0; s < this->nb_sub_materials; ++s) {
      igfem_model->getSubElementBarycenter(el.element, s, el.type, barycenter,
                                           ghost_type);
      sub_material_index(s) =
          1 - igfem_model->isInside(barycenter, this->name_sub_mat_1);
    }

    for (UInt q = 0; q < quads_1; ++q, ++sub_mat_ptr) {
      UInt index = sub_material_index(0);
      *sub_mat_ptr = index;
    }

    for (UInt q = 0; q < quads_2; ++q, ++sub_mat_ptr) {
      UInt index = sub_material_index(1);
      *sub_mat_ptr = index;
    }
  }
}

/* -------------------------------------------------------------------------- */
void MaterialIGFEM::applyEigenGradU(
    const Matrix<Real> & prescribed_eigen_grad_u, const ID & id,
    const GhostType ghost_type) {

  std::map<UInt, ID>::const_iterator sub_mat_it =
      this->sub_material_names.begin();
  for (; sub_mat_it != sub_material_names.end(); ++sub_mat_it) {
    if (sub_mat_it->second == id) {
      UInt sub_element_index = sub_mat_it->first;

      ElementTypeMapArray<UInt>::type_iterator it =
          this->element_filter.firstType(_all_dimensions, ghost_type,
                                         _ek_not_defined);
      ElementTypeMapArray<UInt>::type_iterator end =
          element_filter.lastType(_all_dimensions, ghost_type, _ek_not_defined);

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
        UInt * sub_mat_ptr = this->sub_material(type, ghost_type).storage();

        for (; eigen_it != eigen_end; ++eigen_it, ++sub_mat_ptr) {
          if (*sub_mat_ptr == sub_element_index) {
            Matrix<Real> & current_eigengradu = *eigen_it;
            current_eigengradu = prescribed_eigen_grad_u;
          }
        }
      }
    }
  }
}

} // namespace akantu
