/**
 * @file   material_selector_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Nov 13 2013
 * @date last modification: Thu Jun 05 2014
 *
 * @brief  Implementation of the template MaterialSelector
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef __AKANTU_MATERIAL_SELECTOR_TMPL_HH__
#define __AKANTU_MATERIAL_SELECTOR_TMPL_HH__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<>
class MeshDataMaterialSelector<std::string> : public MaterialSelector {
public:
  MeshDataMaterialSelector(const std::string & name, const SolidMechanicsModel & model) :
    names(model.getMesh().getData<std::string>(name)), model(model) { }
  UInt operator() (const Element & element) {
    try {
      DebugLevel dbl = debug::getDebugLevel();
      debug::setDebugLevel(dblError);
      std::string material_name = names(element.type, element.ghost_type)(element.element);
      debug::setDebugLevel(dbl);

      return model.getMaterialIndex(material_name);
    } catch (...) {
      return MaterialSelector::operator()(element);
    }
  }

private:
  const ElementTypeMapArray<std::string> & names;
protected:
  const SolidMechanicsModel & model;
};

/* -------------------------------------------------------------------------- */
template<>
class MeshDataMaterialSelector<UInt> : public MaterialSelector {
public:
  MeshDataMaterialSelector(const std::string & name, const SolidMechanicsModel & model, UInt first_index = 1) :
    indexes(model.getMesh().getData<UInt>(name)), model(model), first_index(first_index) { }
  UInt operator() (const Element & element) {
    try {
      DebugLevel dbl = debug::getDebugLevel();
      debug::setDebugLevel(dblError);
      UInt mat = indexes(element.type, element.ghost_type)(element.element) - first_index;
      debug::setDebugLevel(dbl);
      return mat;
    } catch (...) {
      return MaterialSelector::operator()(element);
    }
  }
protected:
  const ElementTypeMapArray<UInt> indexes;
  const SolidMechanicsModel & model;
  UInt first_index;
};


__END_AKANTU__

#endif /* __AKANTU_MATERIAL_SELECTOR_TMPL_HH__ */
