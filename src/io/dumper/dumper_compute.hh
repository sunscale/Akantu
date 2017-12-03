/**
 * @file   dumper_compute.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Tue Sep 02 2014
 * @date last modification: Tue Jan 19 2016
 *
 * @brief  Field that map a function to another field
 *
 * @section LICENSE
 *
 * Copyright  (©)  2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de Lausanne)
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

#ifndef __AKANTU_DUMPER_COMPUTE_HH__
#define __AKANTU_DUMPER_COMPUTE_HH__
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "dumper_iohelper.hh"
#include "dumper_type_traits.hh"
#include "dumper_field.hh"
#include <io_helper.hh>

/* -------------------------------------------------------------------------- */

namespace akantu {
__BEGIN_AKANTU_DUMPER__

class ComputeFunctorInterface {
public:
  virtual ~ComputeFunctorInterface() = default;

  virtual UInt getDim() = 0;
  virtual UInt getNbComponent(UInt old_nb_comp) = 0;
};

/* -------------------------------------------------------------------------- */

template <typename return_type>
class ComputeFunctorOutput : public ComputeFunctorInterface {
public:
  ComputeFunctorOutput() = default;
  ~ComputeFunctorOutput() override = default;
};

/* -------------------------------------------------------------------------- */
template <typename input_type, typename return_type>
class ComputeFunctor : public ComputeFunctorOutput<return_type> {
public:
  ComputeFunctor() = default;
  ~ComputeFunctor() override = default;

  virtual return_type func(const input_type & d, Element global_index) = 0;
};

/* -------------------------------------------------------------------------- */
template <typename SubFieldCompute, typename _return_type>
class FieldCompute : public Field {
  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using sub_iterator = typename SubFieldCompute::iterator;
  using sub_types = typename SubFieldCompute::types;
  using sub_return_type = typename sub_types::return_type;
  using return_type = _return_type;
  using data_type = typename sub_types::data_type;

  using types = TypeTraits<data_type, return_type, ElementTypeMapArray<data_type> >;

  class iterator {
  public:
    iterator(const sub_iterator & it,
             ComputeFunctor<sub_return_type, return_type> & func)
        : it(it), func(func) {}

    bool operator!=(const iterator & it) const { return it.it != this->it; }
    iterator operator++() {
      ++this->it;
      return *this;
    }

    UInt currentGlobalIndex() { return this->it.currentGlobalIndex(); }

    return_type operator*() { return func.func(*it, it.getCurrentElement()); }

    Element getCurrentElement() { return this->it.getCurrentElement(); }

    UInt element_type() { return this->it.element_type(); }

  protected:
    sub_iterator it;
    ComputeFunctor<sub_return_type, return_type> & func;
  };

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  FieldCompute(SubFieldCompute & cont, ComputeFunctorInterface & func)
      : sub_field(cont),
        func(dynamic_cast<ComputeFunctor<sub_return_type, return_type> &>(
            func)) {
    this->checkHomogeneity();
  };

  ~FieldCompute() override {
    delete &(this->sub_field);
    delete &(this->func);
  }

  void registerToDumper(const std::string & id,
                        iohelper::Dumper & dumper) override {
    dumper.addElemDataField(id, *this);
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
public:
  iterator begin() { return iterator(sub_field.begin(), func); }
  iterator end() { return iterator(sub_field.end(), func); }

  UInt getDim() { return func.getDim(); }

  UInt size() {
    throw;
    // return Functor::size();
    return 0;
  }

  void checkHomogeneity() override { this->homogeneous = true; };

  iohelper::DataType getDataType() {
    return iohelper::getDataType<data_type>();
  }

  /// get the number of components of the hosted field
  ElementTypeMap<UInt>
  getNbComponents(UInt dim = _all_dimensions, GhostType ghost_type = _not_ghost,
                  ElementKind kind = _ek_not_defined) override {
    ElementTypeMap<UInt> nb_components;
    const ElementTypeMap<UInt> & old_nb_components =
        this->sub_field.getNbComponents(dim, ghost_type, kind);

    ElementTypeMap<UInt>::type_iterator tit =
        old_nb_components.firstType(dim, ghost_type, kind);
    ElementTypeMap<UInt>::type_iterator end =
        old_nb_components.lastType(dim, ghost_type, kind);

    while (tit != end) {
      UInt nb_comp = old_nb_components(*tit, ghost_type);
      nb_components(*tit, ghost_type) = func.getNbComponent(nb_comp);
      ++tit;
    }
    return nb_components;
  };

  /// for connection to a FieldCompute
  inline Field * connect(FieldComputeProxy & proxy) override;

  /// for connection to a FieldCompute
  ComputeFunctorInterface * connect(HomogenizerProxy & proxy) override;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
public:
  SubFieldCompute & sub_field;
  ComputeFunctor<sub_return_type, return_type> & func;
};

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */

class FieldComputeProxy {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  FieldComputeProxy(ComputeFunctorInterface & func) : func(func){};

  inline static Field * createFieldCompute(Field * field,
                                           ComputeFunctorInterface & func) {
    /// that looks fishy an object passed as a ref and destroyed at their end of
    /// the function
    FieldComputeProxy compute_proxy(func);
    return field->connect(compute_proxy);
  }

  template <typename T> Field * connectToField(T * ptr) {
    if (dynamic_cast<ComputeFunctorOutput<Vector<Real> > *>(&func)) {
      return this->connectToFunctor<Vector<Real> >(ptr);
    } else if (dynamic_cast<ComputeFunctorOutput<Vector<UInt> > *>(&func)) {
      return this->connectToFunctor<Vector<UInt> >(ptr);
    }

    else if (dynamic_cast<ComputeFunctorOutput<Matrix<UInt> > *>(&func)) {
      return this->connectToFunctor<Matrix<UInt> >(ptr);
    }

    else if (dynamic_cast<ComputeFunctorOutput<Matrix<Real> > *>(&func)) {
      return this->connectToFunctor<Matrix<Real> >(ptr);
    }

    else
      throw;
  }

  template <typename output, typename T> Field * connectToFunctor(T * ptr) {
    auto * functor_ptr = new FieldCompute<T, output>(*ptr, func);
    return functor_ptr;
  }

  template <typename output, typename SubFieldCompute, typename return_type1,
            typename return_type2>
  Field * connectToFunctor(__attribute__((unused)) FieldCompute<
      FieldCompute<SubFieldCompute, return_type1>, return_type2> * ptr) {
    throw; //    return new FieldCompute<T,output>(*ptr,func);
    return nullptr;
  }

  template <typename output, typename SubFieldCompute, typename return_type1,
            typename return_type2, typename return_type3, typename return_type4>
  Field * connectToFunctor(__attribute__((unused)) FieldCompute<
      FieldCompute<FieldCompute<FieldCompute<SubFieldCompute, return_type1>,
                                return_type2>,
                   return_type3>,
      return_type4> * ptr) {
    throw; //    return new FieldCompute<T,output>(*ptr,func);
    return nullptr;
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
public:
  ComputeFunctorInterface & func;
};

/* -------------------------------------------------------------------------- */
/// for connection to a FieldCompute
template <typename SubFieldCompute, typename return_type>
inline Field *
FieldCompute<SubFieldCompute, return_type>::connect(FieldComputeProxy & proxy) {
  return proxy.connectToField(this);
}

/* -------------------------------------------------------------------------- */
__END_AKANTU_DUMPER__
} // akantu

#endif /* __AKANTU_DUMPER_COMPUTE_HH__ */
