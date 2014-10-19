/**
 * @file   solid_mechanics_model_material.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Nov 26 2010
 * @date last modification: Tue Jun 24 2014
 *
 * @brief  instatiation of materials
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
#include "solid_mechanics_model.hh"
#include "material_list.hh"
#include "aka_math.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
#define AKANTU_INTANTIATE_MATERIAL_BY_DIM_NO_TMPL(dim, elem)		\
  registerNewMaterial< BOOST_PP_ARRAY_ELEM(1, elem)< dim > >(section)

#define AKANTU_INTANTIATE_MATERIAL_BY_DIM_TMPL_EACH(r, data, i, elem)	\
  BOOST_PP_EXPR_IF(BOOST_PP_NOT_EQUAL(0, i), else )			\
  if(opt_param == BOOST_PP_STRINGIZE(BOOST_PP_TUPLE_ELEM(2, 0, elem))) { \
    registerNewMaterial< BOOST_PP_ARRAY_ELEM(1, data)< BOOST_PP_ARRAY_ELEM(0, data), \
						       BOOST_PP_SEQ_ENUM(BOOST_PP_TUPLE_ELEM(2, 1, elem)) > >(section); \
    }


#define AKANTU_INTANTIATE_MATERIAL_BY_DIM_TMPL(dim, elem)		\
  BOOST_PP_SEQ_FOR_EACH_I(AKANTU_INTANTIATE_MATERIAL_BY_DIM_TMPL_EACH,	\
			  (2, (dim, BOOST_PP_ARRAY_ELEM(1, elem))),	\
			  BOOST_PP_ARRAY_ELEM(2, elem))			\
  else {								\
    AKANTU_INTANTIATE_MATERIAL_BY_DIM_NO_TMPL(dim, elem);		\
  }

#define AKANTU_INTANTIATE_MATERIAL_BY_DIM(dim, elem)			\
  BOOST_PP_IF(BOOST_PP_EQUAL(3, BOOST_PP_ARRAY_SIZE(elem) ),		\
	      AKANTU_INTANTIATE_MATERIAL_BY_DIM_TMPL,			\
	      AKANTU_INTANTIATE_MATERIAL_BY_DIM_NO_TMPL)(dim, elem)


#define AKANTU_INTANTIATE_MATERIAL(elem)				\
  switch(spatial_dimension) {						\
  case 1: { AKANTU_INTANTIATE_MATERIAL_BY_DIM(1, elem); break; }	\
  case 2: { AKANTU_INTANTIATE_MATERIAL_BY_DIM(2, elem); break; }	\
  case 3: { AKANTU_INTANTIATE_MATERIAL_BY_DIM(3, elem); break; }	\
  }


#define AKANTU_INTANTIATE_MATERIAL_IF(elem)				\
  if (mat_type == BOOST_PP_STRINGIZE(BOOST_PP_ARRAY_ELEM(0, elem))) {	\
    AKANTU_INTANTIATE_MATERIAL(elem);					\
  }

#define AKANTU_INTANTIATE_OTHER_MATERIAL(r, data, elem)			\
  else AKANTU_INTANTIATE_MATERIAL_IF(elem)

#define AKANTU_INSTANTIATE_MATERIALS()					\
  do {									\
    AKANTU_INTANTIATE_MATERIAL_IF(BOOST_PP_SEQ_HEAD(AKANTU_MATERIAL_LIST)) \
      BOOST_PP_SEQ_FOR_EACH(AKANTU_INTANTIATE_OTHER_MATERIAL,		\
			    _,						\
			    BOOST_PP_SEQ_TAIL(AKANTU_MATERIAL_LIST))	\
    else {								\
      if(getStaticParser().isPermissive())				\
	AKANTU_DEBUG_INFO("Malformed material file " <<			\
			  ": unknown material type '"			\
			  << mat_type << "'");				\
      else								\
	AKANTU_DEBUG_WARNING("Malformed material file "			\
			     <<": unknown material type " << mat_type	\
			     << ". This is perhaps a user"		\
			     << " defined material ?");			\
    }									\
  } while(0)

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::instantiateMaterials() {
  std::pair<Parser::const_section_iterator, Parser::const_section_iterator>
    sub_sect = this->parser->getSubSections(_st_material);

  Parser::const_section_iterator it = sub_sect.first;
  for (; it != sub_sect.second; ++it) {
    const ParserSection & section = *it;
    std::string mat_type  = section.getName();
    std::string opt_param = section.getOption();
    AKANTU_INSTANTIATE_MATERIALS();
  }

  are_materials_instantiated = true;
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::initMaterials() {
  AKANTU_DEBUG_ASSERT(materials.size() != 0, "No material to initialize !");

  if(!are_materials_instantiated) instantiateMaterials();

  Material ** mat_val = &(materials.at(0));

  Element element;
  element.ghost_type = _not_ghost;
  Mesh::type_iterator it  = mesh.firstType(spatial_dimension, _not_ghost, _ek_not_defined);
  Mesh::type_iterator end = mesh.lastType(spatial_dimension, _not_ghost, _ek_not_defined);

  // Fill the element material array from the material selector
  for(; it != end; ++it) {
    UInt nb_element = mesh.getNbElement(*it, _not_ghost);
    element.type = *it;
    element.kind = mesh.getKind(element.type);
    Array<UInt> & el_id_by_mat = element_index_by_material(*it, _not_ghost);
    for (UInt el = 0; el < nb_element; ++el) {
      element.element = el;
      UInt mat_index = (*material_selector)(element);
      AKANTU_DEBUG_ASSERT(mat_index < materials.size(),
			  "The material selector returned an index that does not exists");
      el_id_by_mat(el, 0) = mat_index;
    }
  }

  // synchronize the element material arrays
  synch_registry->synchronize(_gst_material_id);


  /// fill the element filters of the materials using the element_material arrays
  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;
    it  = mesh.firstType(spatial_dimension, gt, _ek_not_defined);
    end = mesh.lastType(spatial_dimension, gt, _ek_not_defined);

    for(; it != end; ++it) {
      UInt nb_element = mesh.getNbElement(*it, gt);
      Array<UInt> & el_id_by_mat = element_index_by_material(*it, gt);
      for (UInt el = 0; el < nb_element; ++el) {
	UInt mat_index = el_id_by_mat(el, 0);
	UInt index = mat_val[mat_index]->addElement(*it, el, gt);
	el_id_by_mat(el, 1) = index;
      }
    }
  }

  std::vector<Material *>::iterator mat_it;
  for(mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    /// init internals properties
    (*mat_it)->initMaterial();
  }

  synch_registry->synchronize(_gst_smm_init_mat);

  // initialize mass
  switch(method) {
    case _explicit_lumped_mass:
      assembleMassLumped();
      break;
    case _explicit_consistent_mass:
    case _implicit_dynamic:
      assembleMass();
      break;
    case _static:
      break;
    default:
      AKANTU_EXCEPTION("analysis method not recognised by SolidMechanicsModel");
      break;
  }

  // initialize the previous displacement array if at least on material needs it
  for (mat_it = materials.begin(); mat_it != materials.end(); ++mat_it) {
    Material & mat = **mat_it;
    if (mat.isFiniteDeformation() || mat.isInelasticDeformation()) {
      initArraysPreviousDisplacment();
      break;
    }
  }
}

/* -------------------------------------------------------------------------- */
Int SolidMechanicsModel::getInternalIndexFromID(const ID & id) const {
  AKANTU_DEBUG_IN();

  std::vector<Material *>::const_iterator first = materials.begin();
  std::vector<Material *>::const_iterator last  = materials.end();

  for (; first != last; ++first)
    if ((*first)->getID() == id) {
      AKANTU_DEBUG_OUT();
      return (first - materials.begin());
    }

  AKANTU_DEBUG_OUT();
  return -1;
}

__END_AKANTU__
