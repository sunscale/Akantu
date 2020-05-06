/**
 * @file   dumper_nodal_field.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Oct 26 2012
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  Description of nodal fields
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef __AKANTU_DUMPER_NODAL_FIELD_HH__
#define __AKANTU_DUMPER_NODAL_FIELD_HH__

#include "dumper_field.hh"
#include <io_helper.hh>
/* -------------------------------------------------------------------------- */

namespace akantu {
namespace dumpers {

// This represents a iohelper compatible field
template <typename T, bool filtered = false, class Container = Array<T>,
          class Filter = Array<UInt>>
class NodalField;

/* -------------------------------------------------------------------------- */
template <typename T, class Container, class Filter>
class NodalField<T, false, Container, Filter> : public dumpers::Field {
public:
  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */

  /// associated iterator with any nodal field (non filetered)
  class iterator : public iohelper::iterator<T, iterator, Vector<T>> {
  public:
    iterator(T * vect, UInt offset, UInt n, UInt stride,
             __attribute__((unused)) const UInt * filter = nullptr)
        :

          internal_it(vect), offset(offset), n(n), stride(stride) {}

    bool operator!=(const iterator & it) const override {
      return internal_it != it.internal_it;
    }
    iterator & operator++() override {
      internal_it += offset;
      return *this;
    };
    Vector<T> operator*() override {
      return Vector<T>(internal_it + stride, n);
    };

  private:
    T * internal_it;
    UInt offset, n, stride;
  };

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
  NodalField(const Container & field, UInt n = 0, UInt stride = 0,
             [[gnu::unused]] const Filter * filter = nullptr)
      : field(field), n(n), stride(stride), padding(0) {
    AKANTU_DEBUG_ASSERT(filter == nullptr,
                        "Filter passed to unfiltered NodalField!");
    if (n == 0) {
      this->n = field.getNbComponent() - stride;
    }
  }

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
  void registerToDumper(const std::string & id,
                        iohelper::Dumper & dumper) override {
    dumper.addNodeDataField(id, *this);
  }

  inline iterator begin() {
    return iterator(field.storage(), field.getNbComponent(), n, stride);
  }

  inline iterator end() {
    return iterator(field.storage() + field.getNbComponent() * field.size(),
                    field.getNbComponent(), n, stride);
  }

  bool isHomogeneous() override { return true; }
  void checkHomogeneity() override { this->homogeneous = true; }

  virtual UInt getDim() {
    if (this->padding)
      return this->padding;
    else
      return n;
  }

  void setPadding(UInt padding) { this->padding = padding; }

  UInt size() { return field.size(); }

  iohelper::DataType getDataType() { return iohelper::getDataType<T>(); }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

private:
  const Container & field;
  UInt n, stride;
  UInt padding;
};

/* -------------------------------------------------------------------------- */
template <typename T, class Container, class Filter>
class NodalField<T, true, Container, Filter> : public dumpers::Field {

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */

public:
  class iterator : public iohelper::iterator<T, iterator, Vector<T>> {

  public:
    iterator(T * const vect, UInt _offset, UInt _n, UInt _stride,
             const UInt * filter)
        :

          internal_it(vect), offset(_offset), n(_n), stride(_stride),
          filter(filter) {}

    bool operator!=(const iterator & it) const override {
      return filter != it.filter;
    }

    iterator & operator++() override {
      ++filter;
      return *this;
    }

    Vector<T> operator*() override {
      return Vector<T>(internal_it + *(filter)*offset + stride, n);
    }

  private:
    T * const internal_it;
    UInt offset, n, stride;
    const UInt * filter;
  };

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

  NodalField(const Container & _field, UInt _n = 0, UInt _stride = 0,
             const Filter * filter = NULL)
      : field(_field), n(_n), stride(_stride), filter(filter), padding(0) {
    AKANTU_DEBUG_ASSERT(this->filter != nullptr,
                        "No filter passed to filtered NodalField!");

    AKANTU_DEBUG_ASSERT(this->filter->getNbComponent() == 1,
                        "Multi-component filter given to NodalField ("
                            << this->filter->getNbComponent()
                            << " components detected, sould be 1");
    if (n == 0) {
      this->n = field.getNbComponent() - stride;
    }
  }

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

  void registerToDumper(const std::string & id,
                        iohelper::Dumper & dumper) override {
    dumper.addNodeDataField(id, *this);
  }

  inline iterator begin() {
    return iterator(field.storage(), field.getNbComponent(), n, stride,
                    filter->storage());
  }

  inline iterator end() {
    return iterator(field.storage(), field.getNbComponent(), n, stride,
                    filter->storage() + filter->size());
  }

  bool isHomogeneous() override { return true; }
  void checkHomogeneity() override { this->homogeneous = true; }

  virtual UInt getDim() {
    if (this->padding)
      return this->padding;
    else
      return n;
  }

  void setPadding(UInt padding) { this->padding = padding; }

  UInt size() { return filter->size(); }

  iohelper::DataType getDataType() { return iohelper::getDataType<T>(); }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

private:
  const Container & field;
  UInt n, stride;
  const Filter * filter;

  UInt padding;
};

} // namespace dumpers
} // namespace akantu
/* -------------------------------------------------------------------------- */
#endif /* __AKANTU_DUMPER_NODAL_FIELD_HH__ */
