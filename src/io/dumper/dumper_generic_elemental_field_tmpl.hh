/**
 * @file   dumper_generic_elemental_field_tmpl.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Sep 02 2014
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  Implementation of the template functions of the ElementalField
 *
 * @section LICENSE
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

template <class types, template <class> class iterator>
void GenericElementalField<types, iterator>::checkHomogeneity() {
  using field_type_iterator = typename field_type::type_iterator;
  using array_type = typename field_type::array_type;

  field_type_iterator tit;
  field_type_iterator end;

  tit = field.firstType(spatial_dimension, ghost_type, element_kind);
  end = field.lastType(spatial_dimension, ghost_type, element_kind);

  this->nb_total_element = 0;
  UInt nb_comp = 0;

  bool homogen = true;
  if (tit != end) {
    nb_comp = this->field(*tit, ghost_type).getNbComponent();
    for (; tit != end; ++tit) {
      const array_type & vect = this->field(*tit, ghost_type);
      UInt nb_element = vect.size();
      UInt nb_comp_cur = vect.getNbComponent();
      if (homogen && nb_comp != nb_comp_cur)
        homogen = false;
      this->nb_total_element += nb_element;

      //      this->nb_data_per_elem(*tit,this->ghost_type) = nb_comp_cur;
    }

    if (!homogen)
      nb_comp = 0;
  }

  this->homogeneous = homogen;
}

/* -------------------------------------------------------------------------- */
