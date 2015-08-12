/**
 * @file   element_class_igfem.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 *
 * @brief  Implementation parent material for IGFEM
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

#include "material_igfem.hh"
#include "aka_math.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
MaterialIGFEM::MaterialIGFEM(SolidMechanicsModel & model, const ID & id) :
  Material(model, id),
  sub_material("sub_material", *this),
  name_sub_mat_1(""),
  name_sub_mat_2("") {
  AKANTU_DEBUG_IN();

  this->model = dynamic_cast<SolidMechanicsModelIGFEM*>(&model);
  this->fem = &(model.getFEEngineClass<MyFEEngineIGFEMType>("IGFEMFEEngine"));
 
  this->model->getMesh().initElementTypeMapArray(element_filter,
						 1,
						 spatial_dimension,
						 false,
						 _ek_igfem);

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

  registerParam("name_sub_mat_1"                 , name_sub_mat_1                 , std::string(), _pat_parsable | _pat_readable);
  registerParam("name_sub_mat_2"                 , name_sub_mat_2                 , std::string(), _pat_parsable | _pat_readable);

  this->sub_material.initialize(1);

}

/* -------------------------------------------------------------------------- */
void MaterialIGFEM::computeQuadraturePointsCoordinates(ElementTypeMapArray<Real> & quadrature_points_coordinates,
						       const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  /// compute quadrature points position in undeformed configuration
  Array<Real> nodes_coordinates(model->getIGFEMNodes(), true);

  Mesh::type_iterator it = this->element_filter.firstType(spatial_dimension, ghost_type, _ek_igfem);
  Mesh::type_iterator last_type = this->element_filter.lastType(spatial_dimension, ghost_type, _ek_igfem);
  for(; it != last_type; ++it) {
    const Array<UInt> & elem_filter = this->element_filter(*it, ghost_type);
      UInt nb_element  = elem_filter.getSize();
      if(nb_element) {
      UInt nb_tot_quad = this->fem->getNbQuadraturePoints(*it, ghost_type) * nb_element;

      Array<Real> & quads = quadrature_points_coordinates(*it, ghost_type);
      quads.resize(nb_tot_quad);

      this->model->getFEEngine("IGFEMFEEngine").interpolateOnQuadraturePoints(nodes_coordinates,
							       quads, spatial_dimension,
							       *it, ghost_type, elem_filter);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline ElementTypeMap<UInt> MaterialIGFEM::getInternalDataPerElem(const ID & id,
								  const ElementKind & element_kind,
								  const ID & fe_engine_id) const {
  if (element_kind == _ek_igfem) {
    return Material::getInternalDataPerElem(id, element_kind, "IGFEMFEEngine");
  } else {
    return Material::getInternalDataPerElem(id, element_kind, fe_engine_id);
  }
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

  Mesh::type_iterator it = this->fem->getMesh().firstType(spatial_dimension, ghost_type, _ek_igfem);
  Mesh::type_iterator last_type = this->fem->getMesh().lastType(spatial_dimension, ghost_type, _ek_igfem);

  for(; it != last_type; ++it) {
    Array<UInt> & elem_filter = element_filter(*it, ghost_type);
    if (elem_filter.getSize()) {
      Array<Real> & gradu_vect = gradu(*it, ghost_type);

      /// compute @f$\nabla u@f$
      this->fem->gradientOnQuadraturePoints(model->getDisplacement(), gradu_vect,
						  spatial_dimension,
						  *it, ghost_type, elem_filter);

      gradu_vect -= eigengradu(*it, ghost_type);

      /// compute @f$\mathbf{\sigma}_q@f$ from @f$\nabla u@f$
      computeStress(*it, ghost_type);
    }
  }

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */



__END_AKANTU__
