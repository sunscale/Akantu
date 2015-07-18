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
  sub_material("sub-material", *this),
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
  this->eigenstrain.setElementKind(_ek_igfem);

  this->gradu.setFEEngine(*fem);
  this->stress.setFEEngine(*fem);
  this->eigenstrain.setFEEngine(*fem);

  registerParam("name_sub_mat_1"                 , name_sub_mat_1                 , std::string(), _pat_parsable | _pat_readable);
  registerParam("name_sub_mat_2"                 , name_sub_mat_2                 , std::string(), _pat_parsable | _pat_readable);

  this->sub_material.initialize(1);

}

/* -------------------------------------------------------------------------- */
void transferInternals(const Array<Element> & new_elements, 
		       const Array<Element> & old_elements,
		       const std::string & field_id,
		       ElementTypeMapArray<Real> & internal_flat,
		       const GhostType ghost_type) {

  // typedef ElementTypeMapArray<UInt>::type_iterator iterator;
  // iterator tit = element_filter.firstType(this->spatial_dimension,
  // 						ghost_type, _ek_igfem);
  // iterator end = element_filter.lastType(this->spatial_dimension,
  // 					       ghost_type, _ek_igfem);
  
  // Array<Element>::const_iterator<Element> el_new_begin = new_elements.begin();
  // Array<Element>::const_iterator<Element> el_new_end   = old_elements.end();
  // Array<Element>::const_iterator<Element> el_old_begin = old_elements.begin();


  // for (; tit != end; ++tit) {
  //   ElementType type = *it;
  //   Element el;
  //   el.type = type;
  //   el.ghost_type = ghost_type;

  //   const Array<UInt> & filter   = this->element_filter(type,ghost_type);
  //   /// get the array
  //   const Array<Real> & dst_vect = this->getArray(field_id, type, ghost_type);
  //   ///@todo ask name of internal flat
  //   const Array<Real> & src_vect = this->internal_flat(field_id,type,ghost_type);

  //   // total number of elements in the mesh for a given type
  //   UInt nb_element_src = this->model->mesh.getNbElement(type,ghost_type);
  //   // number of elements in this material
  //   UInt nb_element = filter.getSize();
  //   // number of quadrature points per elem
  //   UInt nb_quad_per_elem = (dst_vect.getSize()/nb_element);
  //   // number of data per quadrature point
  //   UInt nb_data_per_quad = dst_vect.getNbComponent();

  //   Array<Real>::const_vector_iterator it_src =
  //     src_vect.begin_reinterpret(nb_data,nb_element_src);
  //   Array<Real>::vector_iterator it_dst =
  //     dst_vect.begin_reinterpret(nb_data,nb_element);

  //   /// loop over all the elements in the filter
  //   for (UInt i = 0; i < elem_filter.getSize(); ++i; ++it_dst) {
  //     /// get current element
  //     el.element = elem_filter(i);
  //     /// check if current element is in elements added list
  //     Array<Element>::const_iterator<Element> el_added_it =std::find(el_new_begin, el_new_end, el);
  //     if (el_added_it != el_new_end) {
  // 	Vector<Real> & to_interpolate = *it_src[i];
  // 	Vector<Real> & interpolated = *it_dst;
  // 	UInt sub_element = !this->model->getIntersector.isInside(el.element, 0);
  // 	this->interpolateInternal(type, to_interpolate, interpolated, sub_element);
  //     }
      
  //   }

  // }
  
}

/* -------------------------------------------------------------------------- */
void MaterialIGFEM::computeQuadraturePointsCoordinates(ElementTypeMapArray<Real> & quadrature_points_coordinates,
						       const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  Array<Real> nodes_coordinates(model->getIGFEMNodes(), true);
  nodes_coordinates += this->model->getDisplacement();

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

      gradu_vect -= eigenstrain(*it, ghost_type);

      /// compute @f$\mathbf{\sigma}_q@f$ from @f$\nabla u@f$
      computeStress(*it, ghost_type);
    }
  }

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */



__END_AKANTU__
