/**
 * @file   solid_mechanics_model_tmpl.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Dana Christen <dana.christen@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Fri Nov 13 2015
 *
 * @brief  template part of solid mechanics model
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
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SOLID_MECHANICS_MODEL_TMPL_HH__
#define __AKANTU_SOLID_MECHANICS_MODEL_TMPL_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
// template <typename M>
// Material &
// SolidMechanicsModel::registerNewMaterial(const ParserSection & section) {
//   std::string mat_name;
//   std::string mat_type = section.getName();
//   try {
//     std::string tmp = section.getParameter("name");
//     mat_name = tmp; /** this can seam weird, but there is an ambiguous operator
//              * overload that i couldn't solve. @todo remove the
//              * weirdness of this code
//              */
//   } catch (debug::Exception &) {
//     AKANTU_DEBUG_ERROR(
//         "A material of type \'"
//         << mat_type << "\' in the input file has been defined without a name!");
//   }
//   AKANTU_DEBUG_ASSERT(materials_names_to_id.find(mat_name) ==
//                           materials_names_to_id.end(),
//                       "A material with this name '"
//                           << mat_name << "' has already been registered. "
//                           << "Please use unique names for materials");

//   UInt mat_count = materials.size();
//   materials_names_to_id[mat_name] = mat_count;

//   std::stringstream sstr_mat;
//   sstr_mat << this->id << ":" << mat_count << ":" << mat_type;
//   ID mat_id = sstr_mat.str();

//   std::unique_ptr<Material> material(new M(*this, mat_id));
//   material->parseSection(section);

//   materials.push_back(std::move(material));

//   return *(materials.back());
// }

/* -------------------------------------------------------------------------- */
// template <typename M>
// Material &
// SolidMechanicsModel::registerNewEmptyMaterial(const std::string & mat_name) {

//   AKANTU_DEBUG_ASSERT(materials_names_to_id.find(mat_name) ==
//                           materials_names_to_id.end(),
//                       "A material with this name '"
//                           << mat_name << "' has already been registered. "
//                           << "Please use unique names for materials");

//   UInt mat_count = materials.size();
//   materials_names_to_id[mat_name] = mat_count;

//   std::stringstream sstr_mat;
//   sstr_mat << id << ":" << mat_count << ":" << mat_name;
//   ID mat_id = sstr_mat.str();

//   std::unique_ptr<Material> material(new M(*this, mat_id));
//   materials.push_back(std::move(material));
//   return *(materials.back());
// }

/* -------------------------------------------------------------------------- */
// template <typename M>
// void SolidMechanicsModel::registerNewCustomMaterials(const ID & mat_type) {
//   std::pair<Parser::const_section_iterator, Parser::const_section_iterator>
//       sub_sect = getStaticParser().getSubSections(_st_material);

//   Parser::const_section_iterator it = sub_sect.first;
//   for (; it != sub_sect.second; ++it) {
//     if (it->getName() == mat_type) {
//       registerNewMaterial<M>(*it);
//     }
//   }
// }

/* -------------------------------------------------------------------------- */
} // akantu

#endif /* __AKANTU_SOLID_MECHANICS_MODEL_TMPL_HH__ */
