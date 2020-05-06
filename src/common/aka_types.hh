/**
 * @file   aka_types.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Feb 17 2011
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  description of the "simple" types
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

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_math.hh"
/* -------------------------------------------------------------------------- */
#include <initializer_list>
#include <iomanip>
#include <type_traits>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_AKA_TYPES_HH__
#define __AKANTU_AKA_TYPES_HH__

namespace akantu {

enum NormType { L_1 = 1, L_2 = 2, L_inf = UInt(-1) };

/**
 * DimHelper is a class to generalize the setup of a dim array from 3
 * values. This gives a common interface in the TensorStorage class
 * independently of its derived inheritance (Vector, Matrix, Tensor3)
 * @tparam dim
 */
template <UInt dim> struct DimHelper {
  static inline void setDims(UInt m, UInt n, UInt p, UInt dims[dim]);
};

/* -------------------------------------------------------------------------- */
template <> struct DimHelper<1> {
  static inline void setDims(UInt m, __attribute__((unused)) UInt n,
                             __attribute__((unused)) UInt p, UInt dims[1]) {
    dims[0] = m;
  }
};

/* -------------------------------------------------------------------------- */
template <> struct DimHelper<2> {
  static inline void setDims(UInt m, UInt n, __attribute__((unused)) UInt p,
                             UInt dims[2]) {
    dims[0] = m;
    dims[1] = n;
  }
};

/* -------------------------------------------------------------------------- */
template <> struct DimHelper<3> {
  static inline void setDims(UInt m, UInt n, UInt p, UInt dims[3]) {
    dims[0] = m;
    dims[1] = n;
    dims[2] = p;
  }
};

/* -------------------------------------------------------------------------- */
template <typename T, UInt ndim, class RetType> class TensorStorage;

/* -------------------------------------------------------------------------- */
/* Proxy classes                                                              */
/* -------------------------------------------------------------------------- */
namespace tensors {
  template <class A, class B> struct is_copyable {
    enum : bool { value = false };
  };

  template <class A> struct is_copyable<A, A> {
    enum : bool { value = true };
  };

  template <class A> struct is_copyable<A, typename A::RetType> {
    enum : bool { value = true };
  };

  template <class A> struct is_copyable<A, typename A::RetType::proxy> {
    enum : bool { value = true };
  };

} // namespace tensors

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
namespace types {
  namespace details {
    template <typename reference_> class vector_iterator {
    public:
      using difference_type = std::ptrdiff_t;
      using value_type = std::decay_t<reference_>;
      using pointer = value_type *;
      using reference = reference_;
      using iterator_category = std::input_iterator_tag;

      vector_iterator(pointer ptr) : ptr(ptr) {}

      // input iterator ++it
      vector_iterator & operator++() {
        ++ptr;
        return *this;
      }
      // input iterator it++
      vector_iterator operator++(int) {
        auto cpy = *this;
        ++ptr;
        return cpy;
      }
      vector_iterator & operator+=(int n) {
        ptr += n;
        return *this;
      }

      vector_iterator operator+(int n) {
        vector_iterator cpy(*this);
        cpy += n;
        return cpy;
      }

      // input iterator it != other_it
      bool operator!=(const vector_iterator & other) const {
        return ptr != other.ptr;
      }
      bool operator==(const vector_iterator & other) const {
        return ptr == other.ptr;
      }

      difference_type operator-(const vector_iterator & other) const {
        return this->ptr - other.ptr;
      }

      // input iterator dereference *it
      reference operator*() { return *ptr; }
      pointer operator->() { return ptr; }

    private:
      pointer ptr;
    };
  } // namespace details
} // namespace types

/**
 * @class TensorProxy aka_types.hh
 * The TensorProxy class is a proxy class to the
 * TensorStorage it handles the wrapped case. That is to say if an accessor
 * should give access to a Tensor wrapped on some data, like the
 * Array<T>::iterator they can return a TensorProxy that will be automatically
 * transformed as a TensorStorage wrapped on the same data
 * @tparam T stored type
 * @tparam ndim order of the tensor
 * @tparam _RetType real derived type
 */
template <typename T, UInt ndim, class _RetType>
class TensorProxy : public TensorProxyTrait {
protected:
  using RetTypeProxy = typename _RetType::proxy;

  constexpr TensorProxy(T * data, UInt m, UInt n, UInt p) {
    DimHelper<ndim>::setDims(m, n, p, this->n);
    this->values = data;
  }

  template <class Other, typename = std::enable_if_t<
                             tensors::is_copyable<TensorProxy, Other>::value>>
  explicit TensorProxy(const Other & other) {
    this->values = other.storage();
    for (UInt i = 0; i < ndim; ++i)
      this->n[i] = other.size(i);
  }

public:
  using RetType = _RetType;

  UInt size(UInt i) const {
    AKANTU_DEBUG_ASSERT(i < ndim, "This tensor has only " << ndim
                                                          << " dimensions, not "
                                                          << (i + 1));
    return n[i];
  }

  inline UInt size() const {
    UInt _size = 1;
    for (UInt d = 0; d < ndim; ++d)
      _size *= this->n[d];
    return _size;
  }

  T * storage() const { return values; }

  template <class Other, typename = std::enable_if_t<
                             tensors::is_copyable<TensorProxy, Other>::value>>
  inline TensorProxy & operator=(const Other & other) {
    AKANTU_DEBUG_ASSERT(
        other.size() == this->size(),
        "You are trying to copy two tensors with different sizes");
    memcpy(this->values, other.storage(), this->size() * sizeof(T));
    return *this;
  }
  // template <class Other, typename = std::enable_if_t<
  //                          tensors::is_copyable<TensorProxy, Other>::value>>
  // inline TensorProxy & operator=(const Other && other) {
  //   AKANTU_DEBUG_ASSERT(
  //       other.size() == this->size(),
  //       "You are trying to copy two tensors with different sizes");
  //   memcpy(this->values, other.storage(), this->size() * sizeof(T));
  //   return *this;
  // }

  template <typename O> inline RetTypeProxy & operator*=(const O & o) {
    RetType(*this) *= o;
    return static_cast<RetTypeProxy &>(*this);
  }

  template <typename O> inline RetTypeProxy & operator/=(const O & o) {
    RetType(*this) /= o;
    return static_cast<RetTypeProxy &>(*this);
  }

protected:
  T * values;
  UInt n[ndim];
};

/* -------------------------------------------------------------------------- */
template <typename T> class VectorProxy : public TensorProxy<T, 1, Vector<T>> {
  using parent = TensorProxy<T, 1, Vector<T>>;
  using type = Vector<T>;

public:
  constexpr VectorProxy(T * data, UInt n) : parent(data, n, 0, 0) {}
  template <class Other> explicit VectorProxy(Other & src) : parent(src) {}

  /* ---------------------------------------------------------------------- */
  template <class Other>
  inline VectorProxy<T> & operator=(const Other & other) {
    parent::operator=(other);
    return *this;
  }

  // inline VectorProxy<T> & operator=(const VectorProxy && other) {
  //   parent::operator=(other);
  //   return *this;
  // }
  using iterator = types::details::vector_iterator<T &>;
  using const_iterator = types::details::vector_iterator<const T &>;

  iterator begin() { return iterator(this->storage()); }
  iterator end() { return iterator(this->storage() + this->size()); }

  const_iterator begin() const { return const_iterator(this->storage()); }
  const_iterator end() const {
    return const_iterator(this->storage() + this->size());
  }

  /* ------------------------------------------------------------------------ */
  T & operator()(UInt index) { return this->values[index]; };
  const T & operator()(UInt index) const { return this->values[index]; };
};

template <typename T> class MatrixProxy : public TensorProxy<T, 2, Matrix<T>> {
  using parent = TensorProxy<T, 2, Matrix<T>>;
  using type = Matrix<T>;

public:
  MatrixProxy(T * data, UInt m, UInt n) : parent(data, m, n, 0) {}
  template <class Other> explicit MatrixProxy(Other & src) : parent(src) {}

  /* ---------------------------------------------------------------------- */
  template <class Other>
  inline MatrixProxy<T> & operator=(const Other & other) {
    parent::operator=(other);
    return *this;
  }
};

template <typename T>
class Tensor3Proxy : public TensorProxy<T, 3, Tensor3<T>> {
  using parent = TensorProxy<T, 3, Tensor3<T>>;
  using type = Tensor3<T>;

public:
  Tensor3Proxy(const T * data, UInt m, UInt n, UInt k)
      : parent(data, m, n, k) {}
  Tensor3Proxy(const Tensor3Proxy & src) : parent(src) {}
  Tensor3Proxy(const Tensor3<T> & src) : parent(src) {}

  /* ---------------------------------------------------------------------- */
  template <class Other>
  inline Tensor3Proxy<T> & operator=(const Other & other) {
    parent::operator=(other);
    return *this;
  }
};

/* -------------------------------------------------------------------------- */
/* Tensor base class                                                          */
/* -------------------------------------------------------------------------- */
template <typename T, UInt ndim, class RetType>
class TensorStorage : public TensorTrait {
public:
  using value_type = T;

  friend class Array<T>;

protected:
  template <class TensorType> void copySize(const TensorType & src) {
    for (UInt d = 0; d < ndim; ++d)
      this->n[d] = src.size(d);
    this->_size = src.size();
  }

  TensorStorage() : values(nullptr) {
    for (UInt d = 0; d < ndim; ++d)
      this->n[d] = 0;
    _size = 0;
  }

  TensorStorage(const TensorProxy<T, ndim, RetType> & proxy) {
    this->copySize(proxy);
    this->values = proxy.storage();
    this->wrapped = true;
  }

public:
  TensorStorage(const TensorStorage & src) = delete;

  TensorStorage(const TensorStorage & src, bool deep_copy) : values(nullptr) {
    if (deep_copy)
      this->deepCopy(src);
    else
      this->shallowCopy(src);
  }

protected:
  TensorStorage(UInt m, UInt n, UInt p, const T & def) {
    static_assert(std::is_trivially_constructible<T>{},
                  "Cannot create a tensor on non trivial types");
    DimHelper<ndim>::setDims(m, n, p, this->n);

    this->computeSize();
    this->values = new T[this->_size];
    this->set(def);
    this->wrapped = false;
  }

  TensorStorage(T * data, UInt m, UInt n, UInt p) {
    DimHelper<ndim>::setDims(m, n, p, this->n);

    this->computeSize();
    this->values = data;
    this->wrapped = true;
  }

public:
  /* ------------------------------------------------------------------------ */
  template <class TensorType> inline void shallowCopy(const TensorType & src) {
    this->copySize(src);
    if (!this->wrapped)
      delete[] this->values;
    this->values = src.storage();
    this->wrapped = true;
  }

  /* ------------------------------------------------------------------------ */
  template <class TensorType> inline void deepCopy(const TensorType & src) {
    this->copySize(src);

    if (!this->wrapped)
      delete[] this->values;

    static_assert(std::is_trivially_constructible<T>{},
                  "Cannot create a tensor on non trivial types");
    this->values = new T[this->_size];

    static_assert(std::is_trivially_copyable<T>{},
                  "Cannot copy a tensor on non trivial types");
    memcpy((void *)this->values, (void *)src.storage(),
           this->_size * sizeof(T));

    this->wrapped = false;
  }

  virtual ~TensorStorage() {
    if (!this->wrapped)
      delete[] this->values;
  }

  /* ------------------------------------------------------------------------ */
  inline TensorStorage & operator=(const TensorStorage & src) {
    return this->operator=(aka::as_type<RetType>(src));
  }

  /* ------------------------------------------------------------------------ */
  inline TensorStorage & operator=(const RetType & src) {
    if (this != &src) {
      if (this->wrapped) {
        static_assert(std::is_trivially_copyable<T>{},
                      "Cannot copy a tensor on non trivial types");
        // this test is not sufficient for Tensor of order higher than 1
        AKANTU_DEBUG_ASSERT(this->_size == src.size(),
                            "Tensors of different size ("
                                << this->_size << " != " << src.size() << ")");
        memcpy((void *)this->values, (void *)src.storage(),
               this->_size * sizeof(T));
      } else {
        deepCopy(src);
      }
    }
    return *this;
  }

  /* ------------------------------------------------------------------------ */
  template <class R>
  inline RetType & operator+=(const TensorStorage<T, ndim, R> & other) {
    T * a = this->storage();
    T * b = other.storage();
    AKANTU_DEBUG_ASSERT(
        _size == other.size(),
        "The two tensors do not have the same size, they cannot be subtracted");
    for (UInt i = 0; i < _size; ++i)
      *(a++) += *(b++);
    return *(static_cast<RetType *>(this));
  }

  /* ------------------------------------------------------------------------ */
  template <class R>
  inline RetType & operator-=(const TensorStorage<T, ndim, R> & other) {
    T * a = this->storage();
    T * b = other.storage();
    AKANTU_DEBUG_ASSERT(
        _size == other.size(),
        "The two tensors do not have the same size, they cannot be subtracted");
    for (UInt i = 0; i < _size; ++i)
      *(a++) -= *(b++);
    return *(static_cast<RetType *>(this));
  }

  /* ------------------------------------------------------------------------ */
  inline RetType & operator+=(const T & x) {
    T * a = this->values;
    for (UInt i = 0; i < _size; ++i)
      *(a++) += x;
    return *(static_cast<RetType *>(this));
  }

  /* ------------------------------------------------------------------------ */
  inline RetType & operator-=(const T & x) {
    T * a = this->values;
    for (UInt i = 0; i < _size; ++i)
      *(a++) -= x;
    return *(static_cast<RetType *>(this));
  }

  /* ------------------------------------------------------------------------ */
  inline RetType & operator*=(const T & x) {
    T * a = this->storage();
    for (UInt i = 0; i < _size; ++i)
      *(a++) *= x;
    return *(static_cast<RetType *>(this));
  }

  /* ---------------------------------------------------------------------- */
  inline RetType & operator/=(const T & x) {
    T * a = this->values;
    for (UInt i = 0; i < _size; ++i)
      *(a++) /= x;
    return *(static_cast<RetType *>(this));
  }

  /// \f[Y = \alpha X + Y\f]
  inline RetType & aXplusY(const TensorStorage & other, const T & alpha = 1.) {
    AKANTU_DEBUG_ASSERT(
        _size == other.size(),
        "The two tensors do not have the same size, they cannot be subtracted");

    Math::aXplusY(this->_size, alpha, other.storage(), this->storage());
    return *(static_cast<RetType *>(this));
  }

  /* ------------------------------------------------------------------------ */
  T * storage() const { return values; }
  UInt size() const { return _size; }
  UInt size(UInt i) const {
    AKANTU_DEBUG_ASSERT(i < ndim, "This tensor has only " << ndim
                                                          << " dimensions, not "
                                                          << (i + 1));
    return n[i];
  };
  /* ------------------------------------------------------------------------ */
  inline void clear() { memset(values, 0, _size * sizeof(T)); };
  inline void set(const T & t) { std::fill_n(values, _size, t); };

  template <class TensorType> inline void copy(const TensorType & other) {
    AKANTU_DEBUG_ASSERT(
        _size == other.size(),
        "The two tensors do not have the same size, they cannot be copied");
    memcpy(values, other.storage(), _size * sizeof(T));
  }

  bool isWrapped() const { return this->wrapped; }

protected:
  inline void computeSize() {
    _size = 1;
    for (UInt d = 0; d < ndim; ++d)
      _size *= this->n[d];
  }

protected:
  template <typename R, NormType norm_type> struct NormHelper {
    template <class Ten> static R norm(const Ten & ten) {
      R _norm = 0.;
      R * it = ten.storage();
      R * end = ten.storage() + ten.size();
      for (; it < end; ++it)
        _norm += std::pow(std::abs(*it), norm_type);
      return std::pow(_norm, 1. / norm_type);
    }
  };

  template <typename R> struct NormHelper<R, L_1> {
    template <class Ten> static R norm(const Ten & ten) {
      R _norm = 0.;
      R * it = ten.storage();
      R * end = ten.storage() + ten.size();
      for (; it < end; ++it)
        _norm += std::abs(*it);
      return _norm;
    }
  };

  template <typename R> struct NormHelper<R, L_2> {
    template <class Ten> static R norm(const Ten & ten) {
      R _norm = 0.;
      R * it = ten.storage();
      R * end = ten.storage() + ten.size();
      for (; it < end; ++it)
        _norm += *it * *it;
      return sqrt(_norm);
    }
  };

  template <typename R> struct NormHelper<R, L_inf> {
    template <class Ten> static R norm(const Ten & ten) {
      R _norm = 0.;
      R * it = ten.storage();
      R * end = ten.storage() + ten.size();
      for (; it < end; ++it)
        _norm = std::max(std::abs(*it), _norm);
      return _norm;
    }
  };

public:
  /*----------------------------------------------------------------------- */
  /// "Entrywise" norm norm<L_p> @f[ \|\boldsymbol{T}\|_p = \left(
  /// \sum_i^{n[0]}\sum_j^{n[1]}\sum_k^{n[2]} |T_{ijk}|^p \right)^{\frac{1}{p}}
  /// @f]
  template <NormType norm_type> inline T norm() const {
    return NormHelper<T, norm_type>::norm(*this);
  }

protected:
  UInt n[ndim];
  UInt _size;
  T * values;
  bool wrapped{false};
};

/* -------------------------------------------------------------------------- */
/* Vector                                                                     */
/* -------------------------------------------------------------------------- */
template <typename T> class Vector : public TensorStorage<T, 1, Vector<T>> {
  using parent = TensorStorage<T, 1, Vector<T>>;

public:
  using value_type = typename parent::value_type;
  using proxy = VectorProxy<T>;

public:
  Vector() : parent() {}
  explicit Vector(UInt n, const T & def = T()) : parent(n, 0, 0, def) {}
  Vector(T * data, UInt n) : parent(data, n, 0, 0) {}
  Vector(const Vector & src, bool deep_copy = true) : parent(src, deep_copy) {}
  Vector(const TensorProxy<T, 1, Vector> & src) : parent(src) {}

  Vector(std::initializer_list<T> list) : parent(list.size(), 0, 0, T()) {
    UInt i = 0;
    for (auto val : list) {
      operator()(i++) = val;
    }
  }

public:
  using iterator = types::details::vector_iterator<T &>;
  using const_iterator = types::details::vector_iterator<const T &>;

  iterator begin() { return iterator(this->storage()); }
  iterator end() { return iterator(this->storage() + this->size()); }

  const_iterator begin() const { return const_iterator(this->storage()); }
  const_iterator end() const {
    return const_iterator(this->storage() + this->size());
  }

public:
  ~Vector() override = default;

  /* ------------------------------------------------------------------------ */
  inline Vector & operator=(const Vector & src) {
    parent::operator=(src);
    return *this;
  }

  /* ------------------------------------------------------------------------ */
  inline T & operator()(UInt i) {
    AKANTU_DEBUG_ASSERT((i < this->n[0]),
                        "Access out of the vector! "
                            << "Index (" << i
                            << ") is out of the vector of size (" << this->n[0]
                            << ")");
    return *(this->values + i);
  }

  inline const T & operator()(UInt i) const {
    AKANTU_DEBUG_ASSERT((i < this->n[0]),
                        "Access out of the vector! "
                            << "Index (" << i
                            << ") is out of the vector of size (" << this->n[0]
                            << ")");
    return *(this->values + i);
  }

  inline T & operator[](UInt i) { return this->operator()(i); }
  inline const T & operator[](UInt i) const { return this->operator()(i); }

  /* ------------------------------------------------------------------------ */
  inline Vector<T> & operator*=(Real x) { return parent::operator*=(x); }
  inline Vector<T> & operator/=(Real x) { return parent::operator/=(x); }
  /* ------------------------------------------------------------------------ */
  inline Vector<T> & operator*=(const Vector<T> & vect) {
    AKANTU_DEBUG_ASSERT(this->_size == vect._size,
                        "The vectors have non matching sizes");
    T * a = this->storage();
    T * b = vect.storage();
    for (UInt i = 0; i < this->_size; ++i)
      *(a++) *= *(b++);
    return *this;
  }

  /* ------------------------------------------------------------------------ */
  inline Real dot(const Vector<T> & vect) const {
    return Math::vectorDot(this->values, vect.storage(), this->_size);
  }

  /* ------------------------------------------------------------------------ */
  inline Real mean() const {
    Real mean = 0;
    T * a = this->storage();
    for (UInt i = 0; i < this->_size; ++i)
      mean += *(a++);
    return mean / this->_size;
  }

  /* ------------------------------------------------------------------------ */
  inline Vector & crossProduct(const Vector<T> & v1, const Vector<T> & v2) {
    AKANTU_DEBUG_ASSERT(this->size() == 3,
                        "crossProduct is only defined in 3D (n=" << this->size()
                                                                 << ")");
    AKANTU_DEBUG_ASSERT(
        this->size() == v1.size() && this->size() == v2.size(),
        "crossProduct is not a valid operation non matching size vectors");
    Math::vectorProduct3(v1.storage(), v2.storage(), this->values);
    return *this;
  }

  inline Vector crossProduct(const Vector<T> & v) {
    Vector<T> tmp(this->size());
    tmp.crossProduct(*this, v);
    return tmp;
  }

  /* ------------------------------------------------------------------------ */
  inline void solve(const Matrix<T> & A, const Vector<T> & b) {
    AKANTU_DEBUG_ASSERT(
        this->size() == A.rows() && this->_size == A.cols(),
        "The size of the solution vector mismatches the size of the matrix");
    AKANTU_DEBUG_ASSERT(
        this->_size == b._size,
        "The rhs vector has a mismatch in size with the matrix");
    Math::solve(this->_size, A.storage(), this->values, b.storage());
  }

  /* ------------------------------------------------------------------------ */
  template <bool tr_A>
  inline void mul(const Matrix<T> & A, const Vector<T> & x, T alpha = T(1));
  /* ------------------------------------------------------------------------ */
  inline Real norm() const { return parent::template norm<L_2>(); }

  template <NormType nt> inline Real norm() const {
    return parent::template norm<nt>();
  }

  /* ------------------------------------------------------------------------ */
  inline Vector<Real> & normalize() {
    Real n = norm();
    operator/=(n);
    return *this;
  }

  /* ------------------------------------------------------------------------ */
  /// norm of (*this - x)
  inline Real distance(const Vector<T> & y) const {
    Real * vx = this->values;
    Real * vy = y.storage();
    Real sum_2 = 0;
    for (UInt i = 0; i < this->_size; ++i, ++vx, ++vy)
      sum_2 += (*vx - *vy) * (*vx - *vy);
    return sqrt(sum_2);
  }

  /* ------------------------------------------------------------------------ */
  inline bool equal(const Vector<T> & v,
                    Real tolerance = Math::getTolerance()) const {
    T * a = this->storage();
    T * b = v.storage();
    UInt i = 0;
    while (i < this->_size && (std::abs(*(a++) - *(b++)) < tolerance))
      ++i;
    return i == this->_size;
  }

  /* ------------------------------------------------------------------------ */
  inline short compare(const Vector<T> & v,
                       Real tolerance = Math::getTolerance()) const {
    T * a = this->storage();
    T * b = v.storage();
    for (UInt i(0); i < this->_size; ++i, ++a, ++b) {
      if (std::abs(*a - *b) > tolerance)
        return (((*a - *b) > tolerance) ? 1 : -1);
    }
    return 0;
  }

  /* ------------------------------------------------------------------------ */
  inline bool operator==(const Vector<T> & v) const { return equal(v); }
  inline bool operator!=(const Vector<T> & v) const { return !operator==(v); }
  inline bool operator<(const Vector<T> & v) const { return compare(v) == -1; }
  inline bool operator>(const Vector<T> & v) const { return compare(v) == 1; }

  template <typename Func, typename Acc>
  decltype(auto) accumulate(const Vector<T> & v, Acc && accumulator,
                            Func && func) const {
    T * a = this->storage();
    T * b = v.storage();
    for (UInt i(0); i < this->_size; ++i, ++a, ++b) {
      accumulator = func(*a, *b, std::forward<Acc>(accumulator));
    }
    return accumulator;
  }

  inline bool operator<=(const Vector<T> & v) const {
    bool res = true;
    return accumulate(v, res, [](auto && a, auto && b, auto && accumulator) {
      return accumulator & (a <= b);
    });
  }

  inline bool operator>=(const Vector<T> & v) const {
    bool res = true;
    return accumulate(v, res, [](auto && a, auto && b, auto && accumulator) {
      return accumulator & (a >= b);
    });
  }

  /* ------------------------------------------------------------------------ */
  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const {
    std::string space;
    for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
      ;

    stream << "[";
    for (UInt i = 0; i < this->_size; ++i) {
      if (i != 0)
        stream << ", ";
      stream << this->values[i];
    }
    stream << "]";
  }

  /* ---------------------------------------------------------------------- */
  static inline Vector<T> zeros(UInt n) {
    Vector<T> tmp(n);
    tmp.set(T());
    return tmp;
  }
};

using RVector = Vector<Real>;

/* ------------------------------------------------------------------------ */
template <>
inline bool Vector<UInt>::equal(const Vector<UInt> & v,
                                __attribute__((unused)) Real tolerance) const {
  UInt * a = this->storage();
  UInt * b = v.storage();
  UInt i = 0;
  while (i < this->_size && (*(a++) == *(b++)))
    ++i;
  return i == this->_size;
}

/* -------------------------------------------------------------------------- */
namespace types {
  namespace details {
    template <typename Mat> class column_iterator {
    public:
      using difference_type = std::ptrdiff_t;
      using value_type = decltype(std::declval<Mat>().operator()(0));
      using pointer = value_type *;
      using reference = value_type &;
      using iterator_category = std::input_iterator_tag;

      column_iterator(Mat & mat, UInt col) : mat(mat), col(col) {}
      decltype(auto) operator*() { return mat(col); }
      decltype(auto) operator++() {
        ++col;
        AKANTU_DEBUG_ASSERT(col <= mat.cols(), "The iterator is out of bound");
        return *this;
      }
      decltype(auto) operator++(int) {
        auto tmp = *this;
        ++col;
        AKANTU_DEBUG_ASSERT(col <= mat.cols(), "The iterator is out of bound");
        return tmp;
      }
      bool operator!=(const column_iterator & other) const {
        return col != other.col;
      }

      bool operator==(const column_iterator & other) const {
        return not operator!=(other);
      }

    private:
      Mat & mat;
      UInt col;
    };
  } // namespace details
} // namespace types

/* ------------------------------------------------------------------------ */
/* Matrix                                                                   */
/* ------------------------------------------------------------------------ */
template <typename T> class Matrix : public TensorStorage<T, 2, Matrix<T>> {
  using parent = TensorStorage<T, 2, Matrix<T>>;

public:
  using value_type = typename parent::value_type;
  using proxy = MatrixProxy<T>;

public:
  Matrix() : parent() {}
  Matrix(UInt m, UInt n, const T & def = T()) : parent(m, n, 0, def) {}
  Matrix(T * data, UInt m, UInt n) : parent(data, m, n, 0) {}
  Matrix(const Matrix & src, bool deep_copy = true) : parent(src, deep_copy) {}
  Matrix(const MatrixProxy<T> & src) : parent(src) {}
  Matrix(std::initializer_list<std::initializer_list<T>> list) {
    static_assert(std::is_trivially_copyable<T>{},
                  "Cannot create a tensor on non trivial types");
    std::size_t n = 0;
    std::size_t m = list.size();
    for (auto row : list) {
      n = std::max(n, row.size());
    }

    DimHelper<2>::setDims(m, n, 0, this->n);
    this->computeSize();
    this->values = new T[this->_size];
    this->set(0);

    UInt i = 0, j = 0;
    for (auto & row : list) {
      for (auto & val : row) {
        at(i, j++) = val;
      }
      ++i;
      j = 0;
    }
  }

  ~Matrix() override = default;
  /* ------------------------------------------------------------------------ */
  inline Matrix & operator=(const Matrix & src) {
    parent::operator=(src);
    return *this;
  }

public:
  /* ---------------------------------------------------------------------- */
  UInt rows() const { return this->n[0]; }
  UInt cols() const { return this->n[1]; }

  /* ---------------------------------------------------------------------- */
  inline T & at(UInt i, UInt j) {
    AKANTU_DEBUG_ASSERT(((i < this->n[0]) && (j < this->n[1])),
                        "Access out of the matrix! "
                            << "Index (" << i << ", " << j
                            << ") is out of the matrix of size (" << this->n[0]
                            << ", " << this->n[1] << ")");
    return *(this->values + i + j * this->n[0]);
  }

  inline const T & at(UInt i, UInt j) const {
    AKANTU_DEBUG_ASSERT(((i < this->n[0]) && (j < this->n[1])),
                        "Access out of the matrix! "
                            << "Index (" << i << ", " << j
                            << ") is out of the matrix of size (" << this->n[0]
                            << ", " << this->n[1] << ")");
    return *(this->values + i + j * this->n[0]);
  }

  /* ------------------------------------------------------------------------ */
  inline T & operator()(UInt i, UInt j) { return this->at(i, j); }
  inline const T & operator()(UInt i, UInt j) const { return this->at(i, j); }

  /// give a line vector wrapped on the column i
  inline VectorProxy<T> operator()(UInt j) {
    AKANTU_DEBUG_ASSERT(j < this->n[1],
                        "Access out of the matrix! "
                            << "You are trying to access the column vector "
                            << j << " in a matrix of size (" << this->n[0]
                            << ", " << this->n[1] << ")");
    return VectorProxy<T>(this->values + j * this->n[0], this->n[0]);
  }

  inline const VectorProxy<T> operator()(UInt j) const {
    AKANTU_DEBUG_ASSERT(j < this->n[1],
                        "Access out of the matrix! "
                            << "You are trying to access the column vector "
                            << j << " in a matrix of size (" << this->n[0]
                            << ", " << this->n[1] << ")");
    return VectorProxy<T>(this->values + j * this->n[0], this->n[0]);
  }

public:
  decltype(auto) begin() {
    return types::details::column_iterator<Matrix<T>>(*this, 0);
  }
  decltype(auto) begin() const {
    return types::details::column_iterator<const Matrix<T>>(*this, 0);
  }

  decltype(auto) end() {
    return types::details::column_iterator<Matrix<T>>(*this, this->cols());
  }
  decltype(auto) end() const {
    return types::details::column_iterator<const Matrix<T>>(*this,
                                                            this->cols());
  }

  /* ------------------------------------------------------------------------ */
  inline void block(const Matrix & block, UInt pos_i, UInt pos_j) {
    AKANTU_DEBUG_ASSERT(pos_i + block.rows() <= rows(),
                        "The block size or position are not correct");
    AKANTU_DEBUG_ASSERT(pos_i + block.cols() <= cols(),
                        "The block size or position are not correct");
    for (UInt i = 0; i < block.rows(); ++i)
      for (UInt j = 0; j < block.cols(); ++j)
        this->at(i + pos_i, j + pos_j) = block(i, j);
  }

  inline Matrix block(UInt pos_i, UInt pos_j, UInt block_rows,
                      UInt block_cols) const {
    AKANTU_DEBUG_ASSERT(pos_i + block_rows <= rows(),
                        "The block size or position are not correct");
    AKANTU_DEBUG_ASSERT(pos_i + block_cols <= cols(),
                        "The block size or position are not correct");
    Matrix block(block_rows, block_cols);
    for (UInt i = 0; i < block_rows; ++i)
      for (UInt j = 0; j < block_cols; ++j)
        block(i, j) = this->at(i + pos_i, j + pos_j);
    return block;
  }

  inline T & operator[](UInt idx) { return *(this->values + idx); };
  inline const T & operator[](UInt idx) const { return *(this->values + idx); };

  /* ---------------------------------------------------------------------- */
  inline Matrix operator*(const Matrix & B) const {
    Matrix C(this->rows(), B.cols());
    C.mul<false, false>(*this, B);
    return C;
  }

  /* ----------------------------------------------------------------------- */
  inline Matrix & operator*=(const T & x) { return parent::operator*=(x); }

  inline Matrix & operator*=(const Matrix & B) {
    Matrix C(*this);
    this->mul<false, false>(C, B);
    return *this;
  }

  /* ---------------------------------------------------------------------- */
  template <bool tr_A, bool tr_B>
  inline void mul(const Matrix & A, const Matrix & B, T alpha = 1.0) {
    UInt k = A.cols();
    if (tr_A)
      k = A.rows();

#ifndef AKANTU_NDEBUG
    if (tr_B) {
      AKANTU_DEBUG_ASSERT(k == B.cols(),
                          "matrices to multiply have no fit dimensions");
      AKANTU_DEBUG_ASSERT(this->cols() == B.rows(),
                          "matrices to multiply have no fit dimensions");
    } else {
      AKANTU_DEBUG_ASSERT(k == B.rows(),
                          "matrices to multiply have no fit dimensions");
      AKANTU_DEBUG_ASSERT(this->cols() == B.cols(),
                          "matrices to multiply have no fit dimensions");
    }
    if (tr_A) {
      AKANTU_DEBUG_ASSERT(this->rows() == A.cols(),
                          "matrices to multiply have no fit dimensions");
    } else {
      AKANTU_DEBUG_ASSERT(this->rows() == A.rows(),
                          "matrices to multiply have no fit dimensions");
    }
#endif // AKANTU_NDEBUG

    Math::matMul<tr_A, tr_B>(this->rows(), this->cols(), k, alpha, A.storage(),
                             B.storage(), 0., this->storage());
  }

  /* ---------------------------------------------------------------------- */
  inline void outerProduct(const Vector<T> & A, const Vector<T> & B) {
    AKANTU_DEBUG_ASSERT(
        A.size() == this->rows() && B.size() == this->cols(),
        "A and B are not compatible with the size of the matrix");
    for (UInt i = 0; i < this->rows(); ++i) {
      for (UInt j = 0; j < this->cols(); ++j) {
        this->values[i + j * this->rows()] += A[i] * B[j];
      }
    }
  }

private:
  class EigenSorter {
  public:
    EigenSorter(const Vector<T> & eigs) : eigs(eigs) {}

    bool operator()(const UInt & a, const UInt & b) const {
      return (eigs(a) > eigs(b));
    }

  private:
    const Vector<T> & eigs;
  };

public:
  /* ---------------------------------------------------------------------- */
  inline void eig(Vector<T> & eigenvalues, Matrix<T> & eigenvectors,
                  bool sort = true) const {
    AKANTU_DEBUG_ASSERT(this->cols() == this->rows(),
                        "eig is not a valid operation on a rectangular matrix");
    AKANTU_DEBUG_ASSERT(eigenvalues.size() == this->cols(),
                        "eigenvalues should be of size " << this->cols()
                                                         << ".");
#ifndef AKANTU_NDEBUG
    if (eigenvectors.storage() != nullptr)
      AKANTU_DEBUG_ASSERT((eigenvectors.rows() == eigenvectors.cols()) &&
                              (eigenvectors.rows() == this->cols()),
                          "Eigenvectors needs to be a square matrix of size "
                              << this->cols() << " x " << this->cols() << ".");
#endif

    Matrix<T> tmp = *this;
    Vector<T> tmp_eigs(eigenvalues.size());
    Matrix<T> tmp_eig_vects(eigenvectors.rows(), eigenvectors.cols());

    if (tmp_eig_vects.rows() == 0 || tmp_eig_vects.cols() == 0)
      Math::matrixEig(tmp.cols(), tmp.storage(), tmp_eigs.storage());
    else
      Math::matrixEig(tmp.cols(), tmp.storage(), tmp_eigs.storage(),
                      tmp_eig_vects.storage());

    if (not sort) {
      eigenvalues = tmp_eigs;
      eigenvectors = tmp_eig_vects;
      return;
    }

    Vector<UInt> perm(eigenvalues.size());
    for (UInt i = 0; i < perm.size(); ++i)
      perm(i) = i;

    std::sort(perm.storage(), perm.storage() + perm.size(),
              EigenSorter(tmp_eigs));

    for (UInt i = 0; i < perm.size(); ++i)
      eigenvalues(i) = tmp_eigs(perm(i));

    if (tmp_eig_vects.rows() != 0 && tmp_eig_vects.cols() != 0)
      for (UInt i = 0; i < perm.size(); ++i) {
        for (UInt j = 0; j < eigenvectors.rows(); ++j) {
          eigenvectors(j, i) = tmp_eig_vects(j, perm(i));
        }
      }
  }

  /* ---------------------------------------------------------------------- */
  inline void eig(Vector<T> & eigenvalues) const {
    Matrix<T> empty;
    eig(eigenvalues, empty);
  }

  /* ---------------------------------------------------------------------- */
  inline void eye(T alpha = 1.) {
    AKANTU_DEBUG_ASSERT(this->cols() == this->rows(),
                        "eye is not a valid operation on a rectangular matrix");
    this->clear();
    for (UInt i = 0; i < this->cols(); ++i) {
      this->values[i + i * this->rows()] = alpha;
    }
  }

  /* ---------------------------------------------------------------------- */
  static inline Matrix<T> eye(UInt m, T alpha = 1.) {
    Matrix<T> tmp(m, m);
    tmp.eye(alpha);
    return tmp;
  }

  /* ---------------------------------------------------------------------- */
  inline T trace() const {
    AKANTU_DEBUG_ASSERT(
        this->cols() == this->rows(),
        "trace is not a valid operation on a rectangular matrix");
    T trace = 0.;
    for (UInt i = 0; i < this->rows(); ++i) {
      trace += this->values[i + i * this->rows()];
    }
    return trace;
  }

  /* ---------------------------------------------------------------------- */
  inline Matrix transpose() const {
    Matrix tmp(this->cols(), this->rows());
    for (UInt i = 0; i < this->rows(); ++i) {
      for (UInt j = 0; j < this->cols(); ++j) {
        tmp(j, i) = operator()(i, j);
      }
    }
    return tmp;
  }

  /* ---------------------------------------------------------------------- */
  inline void inverse(const Matrix & A) {
    AKANTU_DEBUG_ASSERT(A.cols() == A.rows(),
                        "inv is not a valid operation on a rectangular matrix");
    AKANTU_DEBUG_ASSERT(this->cols() == A.cols(),
                        "the matrix should have the same size as its inverse");

    if (this->cols() == 1)
      *this->values = 1. / *A.storage();
    else if (this->cols() == 2)
      Math::inv2(A.storage(), this->values);
    else if (this->cols() == 3)
      Math::inv3(A.storage(), this->values);
    else
      Math::inv(this->cols(), A.storage(), this->values);
  }

  inline Matrix inverse() {
    Matrix inv(this->rows(), this->cols());
    inv.inverse(*this);
    return inv;
  }

  /* --------------------------------------------------------------------- */
  inline T det() const {
    AKANTU_DEBUG_ASSERT(this->cols() == this->rows(),
                        "inv is not a valid operation on a rectangular matrix");
    if (this->cols() == 1)
      return *(this->values);
    else if (this->cols() == 2)
      return Math::det2(this->values);
    else if (this->cols() == 3)
      return Math::det3(this->values);
    else
      return Math::det(this->cols(), this->values);
  }

  /* --------------------------------------------------------------------- */
  inline T doubleDot(const Matrix<T> & other) const {
    AKANTU_DEBUG_ASSERT(
        this->cols() == this->rows(),
        "doubleDot is not a valid operation on a rectangular matrix");
    if (this->cols() == 1)
      return *(this->values) * *(other.storage());
    else if (this->cols() == 2)
      return Math::matrixDoubleDot22(this->values, other.storage());
    else if (this->cols() == 3)
      return Math::matrixDoubleDot33(this->values, other.storage());
    else
      AKANTU_ERROR("doubleDot is not defined for other spatial dimensions"
                   << " than 1, 2 or 3.");
    return T();
  }

  /* ---------------------------------------------------------------------- */
  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const {
    std::string space;
    for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
      ;

    stream << "[";
    for (UInt i = 0; i < this->n[0]; ++i) {
      if (i != 0)
        stream << ", ";
      stream << "[";
      for (UInt j = 0; j < this->n[1]; ++j) {
        if (j != 0)
          stream << ", ";
        stream << operator()(i, j);
      }
      stream << "]";
    }
    stream << "]";
  };
};

/* ------------------------------------------------------------------------ */
template <typename T>
template <bool tr_A>
inline void Vector<T>::mul(const Matrix<T> & A, const Vector<T> & x, T alpha) {
#ifndef AKANTU_NDEBUG
  UInt n = x.size();
  if (tr_A) {
    AKANTU_DEBUG_ASSERT(n == A.rows(),
                        "matrix and vector to multiply have no fit dimensions");
    AKANTU_DEBUG_ASSERT(this->size() == A.cols(),
                        "matrix and vector to multiply have no fit dimensions");
  } else {
    AKANTU_DEBUG_ASSERT(n == A.cols(),
                        "matrix and vector to multiply have no fit dimensions");
    AKANTU_DEBUG_ASSERT(this->size() == A.rows(),
                        "matrix and vector to multiply have no fit dimensions");
  }
#endif
  Math::matVectMul<tr_A>(A.rows(), A.cols(), alpha, A.storage(), x.storage(),
                         0., this->storage());
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline std::ostream & operator<<(std::ostream & stream,
                                 const Matrix<T> & _this) {
  _this.printself(stream);
  return stream;
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline std::ostream & operator<<(std::ostream & stream,
                                 const Vector<T> & _this) {
  _this.printself(stream);
  return stream;
}

/* ------------------------------------------------------------------------ */
/* Tensor3                                                                  */
/* ------------------------------------------------------------------------ */
template <typename T> class Tensor3 : public TensorStorage<T, 3, Tensor3<T>> {
  using parent = TensorStorage<T, 3, Tensor3<T>>;

public:
  using value_type = typename parent::value_type;
  using proxy = Tensor3Proxy<T>;

public:
  Tensor3() : parent(){};

  Tensor3(UInt m, UInt n, UInt p, const T & def = T()) : parent(m, n, p, def) {}

  Tensor3(T * data, UInt m, UInt n, UInt p) : parent(data, m, n, p) {}

  Tensor3(const Tensor3 & src, bool deep_copy = true)
      : parent(src, deep_copy) {}

  Tensor3(const proxy & src) : parent(src) {}

public:
  /* ------------------------------------------------------------------------ */
  inline Tensor3 & operator=(const Tensor3 & src) {
    parent::operator=(src);
    return *this;
  }

  /* ---------------------------------------------------------------------- */
  inline T & operator()(UInt i, UInt j, UInt k) {
    AKANTU_DEBUG_ASSERT(
        (i < this->n[0]) && (j < this->n[1]) && (k < this->n[2]),
        "Access out of the tensor3! "
            << "You are trying to access the element "
            << "(" << i << ", " << j << ", " << k << ") in a tensor of size ("
            << this->n[0] << ", " << this->n[1] << ", " << this->n[2] << ")");
    return *(this->values + (k * this->n[0] + i) * this->n[1] + j);
  }

  inline const T & operator()(UInt i, UInt j, UInt k) const {
    AKANTU_DEBUG_ASSERT(
        (i < this->n[0]) && (j < this->n[1]) && (k < this->n[2]),
        "Access out of the tensor3! "
            << "You are trying to access the element "
            << "(" << i << ", " << j << ", " << k << ") in a tensor of size ("
            << this->n[0] << ", " << this->n[1] << ", " << this->n[2] << ")");
    return *(this->values + (k * this->n[0] + i) * this->n[1] + j);
  }

  inline MatrixProxy<T> operator()(UInt k) {
    AKANTU_DEBUG_ASSERT((k < this->n[2]),
                        "Access out of the tensor3! "
                            << "You are trying to access the slice " << k
                            << " in a tensor3 of size (" << this->n[0] << ", "
                            << this->n[1] << ", " << this->n[2] << ")");
    return MatrixProxy<T>(this->values + k * this->n[0] * this->n[1],
                          this->n[0], this->n[1]);
  }

  inline const MatrixProxy<T> operator()(UInt k) const {
    AKANTU_DEBUG_ASSERT((k < this->n[2]),
                        "Access out of the tensor3! "
                            << "You are trying to access the slice " << k
                            << " in a tensor3 of size (" << this->n[0] << ", "
                            << this->n[1] << ", " << this->n[2] << ")");
    return MatrixProxy<T>(this->values + k * this->n[0] * this->n[1],
                          this->n[0], this->n[1]);
  }

  inline MatrixProxy<T> operator[](UInt k) {
    return MatrixProxy<T>(this->values + k * this->n[0] * this->n[1],
                          this->n[0], this->n[1]);
  }

  inline const MatrixProxy<T> operator[](UInt k) const {
    return MatrixProxy<T>(this->values + k * this->n[0] * this->n[1],
                          this->n[0], this->n[1]);
  }
};

/* -------------------------------------------------------------------------- */
// support operations for the creation of other vectors
/* -------------------------------------------------------------------------- */
template <typename T>
Vector<T> operator*(const T & scalar, const Vector<T> & a) {
  Vector<T> r(a);
  r *= scalar;
  return r;
}

template <typename T>
Vector<T> operator*(const Vector<T> & a, const T & scalar) {
  Vector<T> r(a);
  r *= scalar;
  return r;
}

template <typename T>
Vector<T> operator/(const Vector<T> & a, const T & scalar) {
  Vector<T> r(a);
  r /= scalar;
  return r;
}

template <typename T>
Vector<T> operator*(const Vector<T> & a, const Vector<T> & b) {
  Vector<T> r(a);
  r *= b;
  return r;
}

template <typename T>
Vector<T> operator+(const Vector<T> & a, const Vector<T> & b) {
  Vector<T> r(a);
  r += b;
  return r;
}

template <typename T>
Vector<T> operator-(const Vector<T> & a, const Vector<T> & b) {
  Vector<T> r(a);
  r -= b;
  return r;
}

template <typename T>
Vector<T> operator*(const Matrix<T> & A, const Vector<T> & b) {
  Vector<T> r(b.size());
  r.template mul<false>(A, b);
  return r;
}

/* -------------------------------------------------------------------------- */
template <typename T>
Matrix<T> operator*(const T & scalar, const Matrix<T> & a) {
  Matrix<T> r(a);
  r *= scalar;
  return r;
}

template <typename T>
Matrix<T> operator*(const Matrix<T> & a, const T & scalar) {
  Matrix<T> r(a);
  r *= scalar;
  return r;
}

template <typename T>
Matrix<T> operator/(const Matrix<T> & a, const T & scalar) {
  Matrix<T> r(a);
  r /= scalar;
  return r;
}

template <typename T>
Matrix<T> operator+(const Matrix<T> & a, const Matrix<T> & b) {
  Matrix<T> r(a);
  r += b;
  return r;
}

template <typename T>
Matrix<T> operator-(const Matrix<T> & a, const Matrix<T> & b) {
  Matrix<T> r(a);
  r -= b;
  return r;
}

} // namespace akantu

#include <iterator>

namespace std {
template <typename R>
struct iterator_traits<::akantu::types::details::vector_iterator<R>> {
protected:
  using iterator = ::akantu::types::details::vector_iterator<R>;

public:
  using iterator_category = typename iterator::iterator_category;
  using value_type = typename iterator::value_type;
  using difference_type = typename iterator::difference_type;
  using pointer = typename iterator::pointer;
  using reference = typename iterator::reference;
};

template <typename Mat>
struct iterator_traits<::akantu::types::details::column_iterator<Mat>> {
protected:
  using iterator = ::akantu::types::details::column_iterator<Mat>;

public:
  using iterator_category = typename iterator::iterator_category;
  using value_type = typename iterator::value_type;
  using difference_type = typename iterator::difference_type;
  using pointer = typename iterator::pointer;
  using reference = typename iterator::reference;
};

} // namespace std

#endif /* __AKANTU_AKA_TYPES_HH__ */
