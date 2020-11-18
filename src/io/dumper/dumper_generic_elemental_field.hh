/**
 * @file   dumper_generic_elemental_field.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Sep 02 2014
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  Generic interface for elemental fields
 *
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

#ifndef AKANTU_DUMPER_GENERIC_ELEMENTAL_FIELD_HH_
#define AKANTU_DUMPER_GENERIC_ELEMENTAL_FIELD_HH_
/* -------------------------------------------------------------------------- */
#include "dumper_element_iterator.hh"
#include "dumper_field.hh"
#include "dumper_homogenizing_field.hh"
#include "element_type_map_filter.hh"
/* -------------------------------------------------------------------------- */
namespace akantu {
namespace dumpers {
/* -------------------------------------------------------------------------- */

template <class _types, template <class> class iterator_type>
class GenericElementalField : public Field {
  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  // check dumper_type_traits.hh for additional information over these types
  using types = _types;
  using data_type = typename types::data_type;
  using it_type = typename types::it_type;
  using field_type = typename types::field_type;
  using array_type = typename types::array_type;
  using array_iterator = typename types::array_iterator;
  using field_type_iterator = typename field_type::type_iterator;
  using iterator = iterator_type<types>;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  GenericElementalField(const field_type & field,
                        UInt spatial_dimension = _all_dimensions,
                        GhostType ghost_type = _not_ghost,
                        ElementKind element_kind = _ek_not_defined)
      : field(field), spatial_dimension(spatial_dimension),
        ghost_type(ghost_type), element_kind(element_kind) {
    this->checkHomogeneity();
  }

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// get the number of components of the hosted field
  ElementTypeMap<UInt>
  getNbComponents(UInt dim = _all_dimensions, GhostType ghost_type = _not_ghost,
                  ElementKind kind = _ek_not_defined) override {
    return this->field.getNbComponents(dim, ghost_type, kind);
  };

  /// return the size of the contained data: i.e. the number of elements ?
  virtual UInt size() {
    checkHomogeneity();
    return this->nb_total_element;
  }

  /// return the iohelper datatype to be dumped
  iohelper::DataType getDataType() {
    return iohelper::getDataType<data_type>();
  }

protected:
  /// return the number of entries per element
  UInt getNbDataPerElem(ElementType type,
                        GhostType ghost_type = _not_ghost) const {
    if (!nb_data_per_elem.exists(type, ghost_type)) {
      return field(type, ghost_type).getNbComponent();
    }

    return nb_data_per_elem(type, this->ghost_type);
  }

  /// check if the same quantity of data for all element types
  void checkHomogeneity() override;

public:
  void registerToDumper(const std::string & id,
                        iohelper::Dumper & dumper) override {
    dumper.addElemDataField(id, *this);
  };

  /// for connection to a FieldCompute
  inline std::shared_ptr<Field> connect(FieldComputeProxy & proxy) override {
    return proxy.connectToField(this);
  }

  /// for connection to a Homogenizer
  inline std::unique_ptr<ComputeFunctorInterface>
  connect(HomogenizerProxy & proxy) override {
    return proxy.connectToField(this);
  };

  virtual iterator begin() {
    /// type iterators on the elemental field
    auto types = this->field.elementTypes(this->spatial_dimension,
                                          this->ghost_type, this->element_kind);
    auto tit = types.begin();
    auto end = types.end();

    /// skip all types without data
    for (; tit != end and this->field(*tit, this->ghost_type).empty();
         ++tit) {
    }

    auto type = *tit;

    if (tit == end) {
      return this->end();
    }

    /// getting information for the field of the given type
    const auto & vect = this->field(type, this->ghost_type);
    UInt nb_data_per_elem = this->getNbDataPerElem(type);

    /// define element-wise iterator
    auto view = make_view(vect, nb_data_per_elem);
    auto it = view.begin();
    auto it_end = view.end();
    /// define data iterator
    iterator rit =
        iterator(this->field, tit, end, it, it_end, this->ghost_type);
    rit.setNbDataPerElem(this->nb_data_per_elem);
    return rit;
  }

  virtual iterator end() {
    auto types = this->field.elementTypes(this->spatial_dimension,
                                          this->ghost_type, this->element_kind);
    auto tit = types.begin();
    auto end = types.end();

    auto type = *tit;
    for (; tit != end; ++tit) {
      type = *tit;
    }

    const array_type & vect = this->field(type, this->ghost_type);
    UInt nb_data = this->getNbDataPerElem(type);
    auto it = make_view(vect, nb_data).end();
    auto rit = iterator(this->field, end, end, it, it, this->ghost_type);
    rit.setNbDataPerElem(this->nb_data_per_elem);

    return rit;
  }

  virtual UInt getDim() {
    if (this->homogeneous) {
      auto tit = this->field
                     .elementTypes(this->spatial_dimension, this->ghost_type,
                                   this->element_kind)
                     .begin();
      return this->getNbDataPerElem(*tit);
    }

    throw;
    return 0;
  }

  void setNbDataPerElem(const ElementTypeMap<UInt> & nb_data) override {
    nb_data_per_elem = nb_data;
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// the ElementTypeMapArray embedded in the field
  const field_type & field;
  /// total number of elements
  UInt nb_total_element;
  /// the spatial dimension of the problem
  UInt spatial_dimension;
  /// whether this is a ghost field or not (for type selection)
  GhostType ghost_type;
  /// The element kind to operate on
  ElementKind element_kind;
  /// The number of data per element type
  ElementTypeMap<UInt> nb_data_per_elem;
};

} // namespace dumpers
} // namespace akantu

/* -------------------------------------------------------------------------- */
#include "dumper_generic_elemental_field_tmpl.hh"
/* -------------------------------------------------------------------------- */

#endif /* AKANTU_DUMPER_GENERIC_ELEMENTAL_FIELD_HH_ */
