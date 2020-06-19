/**
 * @file   test_types.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri May 15 2015
 * @date last modification: Wed Jun 14 2017
 *
 * @brief  Test the types declared in aka_types.hh
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_types.hh"

#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace akantu;

const Real tolerance = 1e-15;

std::string itoa(UInt a) {
  std::stringstream sstr;
  sstr << a;
  return sstr.str();
}

UInt testcounter = 0;

struct wrap_error : std::runtime_error {
  wrap_error(const std::string & msg) : std::runtime_error(msg) {}
};

struct size_error : std::runtime_error {
  size_error(const std::string & msg) : std::runtime_error(msg) {}
};

struct data_error : std::runtime_error {
  data_error(const std::string & msg, UInt i)
      : std::runtime_error(msg), index(i) {}
  UInt index;
};

template <class type>
void compare_storages_with_ref(const type & a, Real * ref, UInt size, UInt line,
                               const std::string & txt) {
  std::cout << std::setw(3) << (testcounter++) << ": " << std::setw(10) << txt
            << " - " << a << " - wrapped: " << std::boolalpha << a.isWrapped()
            << std::endl;

  if (a.size() != size)
    throw size_error("the size is not correct " + itoa(a.size()) +
                     " instead of " + itoa(size) +
                     " [Test at line: " + itoa(line) + "]");

  Real * a_ptr = a.storage();
  for (UInt i = 0; i < a.size(); ++i) {
    if (!((std::abs(a_ptr[i]) < tolerance && std::abs(ref[i]) < tolerance) ||
          std::abs((a_ptr[i] - ref[i]) / a_ptr[i]) < tolerance)) {
      std::stringstream txt;
      txt << " std::abs(" << a_ptr[i] << " - " << ref[i]
          << " [= " << std::abs(a_ptr[i] - ref[i]) << "] ) > " << tolerance;
      throw data_error("storage differs at index " + itoa(i) +
                           " [Test at line: " + itoa(line) + "]" + txt.str(),
                       i);
    }
  }

  if (a_ptr == ref && !a.isWrapped())
    throw wrap_error(
        "the storage should be wrapped but it is not [Test at line: " +
        itoa(line) + "]");
  if (a_ptr != ref && a.isWrapped())
    throw wrap_error(
        "the storage should not be wrapped but it is [Test at line: " +
        itoa(line) + "]");
}

#define COMPARE(a, aref, txt)                                                  \
  compare_storages_with_ref(a, aref, sizeof(aref) / sizeof(aref[0]), __LINE__, \
                            txt)
#define COMPARE_STORAGE(a, aref, txt)                                          \
  compare_storages_with_ref(a, aref.storage(), aref.size(), __LINE__, txt)

const UInt ref_size = 10;

// clang-format off
/* -------------------------------------------------------------------------- */
void test_constructor() {
  std::cout << "=== Test constructors ===" << std::endl;
  Real ref1[ref_size] = { 0. };
  Real ref2[ref_size] = { 1563.58, 1563.58, 1563.58, 1563.58, 1563.58, 1563.58, 1563.58, 1563.58, 1563.58, 1563.58 };
  Real ref3[ref_size] = { 23.1594, 79.6184, 77.9052, 47.9922, 12.8674, 37.1445, 64.8991, 80.3364, 98.4064, 73.7858 };

  std::cout << "--  Vectors: " << std::endl;
  Vector<Real> v0 =  { 23.1594, 79.6184, 77.9052, 47.9922, 12.8674, 37.1445, 64.8991, 80.3364, 98.4064, 73.7858 };
                           ;          COMPARE        (    v0,   ref3, "init_list" );
  Vector<Real> v1(ref_size);          COMPARE        (    v1,   ref1, "normal"    );
  Vector<Real> v2(ref_size, 1563.58); COMPARE        (    v2,   ref2, "defval"    );
  Vector<Real> v3(ref3, ref_size);    COMPARE        (    v3,   ref3, "wrapped"   );
  Vector<Real> v3dcw(v3);             COMPARE        ( v3dcw,   ref3, "wdeepcopy" );
  Vector<Real> v3scw(v3, false);      COMPARE        ( v3scw,   ref3, "wshallow"  );
  Vector<Real> v3dc(v3dcw);           COMPARE_STORAGE(  v3dc,  v3dcw, "deepcopy"  );
  Vector<Real> v3sc(v3dcw, false);    COMPARE_STORAGE(  v3sc,  v3dcw, "shallow"   );
  VectorProxy<Real> vp1(ref3, ref_size);
  Vector<Real> v4(vp1);               COMPARE        (    v4,   ref3, "proxyptr"  );
  VectorProxy<Real> vp2(v3dcw);
  Vector<Real> v5(vp2);               COMPARE_STORAGE(    v5,  v3dcw, "proxyvdc"  );
  VectorProxy<Real> vp3(v3scw);
  Vector<Real> v6(vp3);               COMPARE        (    v6,   ref3, "proxyvsc"  );

  /* ------------------------------------------------------------------------ */
  std::cout << "--  Matrices: " << std::endl;
  Matrix<Real> m0  = {{23.1594, 37.1445},
                      {79.6184, 64.8991},
                      {77.9052, 80.3364},
                      {47.9922, 98.4064},
                      {12.8674, 73.7858}};
                                      COMPARE        (    m0, ref3  , "init_list" );
  Matrix<Real> m1(5, 2);              COMPARE        (    m1, ref1  , "normal"    );
  Matrix<Real> m1t(2, 5);             COMPARE        (   m1t, ref1  , "tnormal"   );
  Matrix<Real> m2(5, 2, 1563.58);     COMPARE        (    m2, ref2  , "defval"    );
  Matrix<Real> m2t(2, 5, 1563.58);    COMPARE        (   m2t, ref2  , "tdefval"   );
  Matrix<Real> m3(ref3, 5, 2);        COMPARE        (    m3, ref3  , "wrapped"   );
  Matrix<Real> m3t(ref3, 2, 5);       COMPARE        (   m3t, ref3  , "twrapped"  );
  Matrix<Real> m3dcw(m3);             COMPARE        ( m3dcw, ref3  , "wdeepcopy" );
  Matrix<Real> m3scw(m3, false);      COMPARE        ( m3scw, ref3  , "wshallow"  );
  Matrix<Real> m3dc(m3dcw);           COMPARE_STORAGE(  m3dc, m3dcw , "deepcopy"  );
  Matrix<Real> m3sc(m3dcw, false);    COMPARE_STORAGE(  m3sc, m3dcw , "shallow"   );
  Matrix<Real> m3tdcw(m3t);           COMPARE        (m3tdcw, ref3  , "twdeepcopy");
  Matrix<Real> m3tscw(m3t, false);    COMPARE        (m3tscw, ref3  , "twshallow" );
  Matrix<Real> m3tdc(m3tdcw);         COMPARE_STORAGE( m3tdc, m3tdcw, "tdeepcopy" );
  Matrix<Real> m3tsc(m3tdcw, false);  COMPARE_STORAGE( m3tsc, m3tdcw, "tshallow"  );
  MatrixProxy<Real> mp1(ref3, 5, 2);
  Matrix<Real> m4(mp1);               COMPARE        (    m4,   ref3, "proxyptr"  );
  MatrixProxy<Real> mp2(m3dcw);
  Matrix<Real> m5(mp2);               COMPARE_STORAGE(    m5,  m3dcw, "proxyvdc"  );
  MatrixProxy<Real> mp3(m3scw);
  Matrix<Real> m6(mp3);               COMPARE        (    m6,   ref3, "proxyvsc"  );
  MatrixProxy<Real> mp1t(ref3, 2, 5);
  Matrix<Real> m4t(mp1t);             COMPARE        (   m4t,   ref3, "tproxyptr" );
  MatrixProxy<Real> mp2t(m3tdcw);
  Matrix<Real> m5t(mp2t);             COMPARE_STORAGE(   m5t, m3tdcw, "tproxyvdc" );
  MatrixProxy<Real> mp3t(m3tscw);
  Matrix<Real> m6t(mp3t);             COMPARE        (   m6t,   ref3, "tproxyvsc" );
}

/* -------------------------------------------------------------------------- */
void test_equal_and_accessors() {
  std::cout << "=== Test operator=() ===" << std::endl;
  Real ref[ref_size] = { 23.1594, 79.6184, 77.9052, 47.9922, 12.8674, 37.1445, 64.8991, 80.3364, 98.4064, 73.7858 };
  Real mod[ref_size] = { 98.7982, 72.1227, 19.7815, 57.6722, 47.1088, 14.9865, 13.3171, 62.7973, 33.9493, 98.3052 };

  std::cout << "--  Vectors: " << std::endl;
  Vector<Real> v (ref, ref_size);
  Vector<Real> vm(mod, ref_size);
  Vector<Real> vref1(v);
  Vector<Real> v1;
  v1 = vref1;                                          COMPARE_STORAGE(v1, vref1, "simple="  );
  for (UInt i = 0; i < ref_size; ++i) v1 (i) = mod[i]; COMPARE        (v1,   mod, "s_acces"   );
  COMPARE_STORAGE(vref1, v, "refcheck1");

  Vector<Real> v2 = vref1;                             COMPARE_STORAGE(v2, vref1, "construc=");
  for (UInt i = 0; i < ref_size; ++i) v2 (i) = mod[i]; COMPARE        (v2,   mod, "c_acces"   );
  COMPARE_STORAGE(vref1, v, "refcheck2");

  Vector<Real> vref2(vref1, false);
  Vector<Real> v1w;
  v1w = vref2;                                         COMPARE_STORAGE(v1w, vref1, "w_simple=" );
  for (UInt i = 0; i < ref_size; ++i) v1w(i) = mod[i]; COMPARE        (v1w,   mod, "ws_acces"  );
  try { COMPARE(vref2, ref, "refcheck3"); } catch(wrap_error &) {}

  Vector<Real> v2w = vref2;                            COMPARE_STORAGE(v2w, vref1, "w_constru=");
  for (UInt i = 0; i < ref_size; ++i) v2w(i) = mod[i]; COMPARE        (v2w,   mod, "wc_acces"  );
  try { COMPARE(vref2, ref, "refcheck4"); } catch(wrap_error &) {}

  VectorProxy<Real> vp1(vref1);
  Vector<Real> v3;
  v3 = vp1;                                             COMPARE_STORAGE(v3, vref1, "p_simple=" );
  for (UInt i = 0; i < ref_size; ++i) v3(i) = mod[i];   COMPARE        (v3,   mod, "ps_acces"  );
  COMPARE_STORAGE(vref1, v, "refcheck5");

  Vector<Real> v4 = vp1;                                COMPARE_STORAGE(v4, vref1, "p_constru=");
  for (UInt i = 0; i < ref_size; ++i) v4(i) = mod[i];
  try { COMPARE(v4,   mod, "pc_acces" ); } catch (wrap_error &) {}

  COMPARE(vref1, mod, "refcheck6");
  try { COMPARE(vref2, mod, "refcheck7"); } catch(wrap_error &) {}

  vref2 = v;

  VectorProxy<Real> vp2(vref2);
  Vector<Real> v3w;
  v3w = vp2;                                           COMPARE_STORAGE(v3w, vref1, "pw_simpl=");
  for (UInt i = 0; i < ref_size; ++i) v3w(i) = mod[i]; COMPARE        (v3w,   mod, "pws_acces");
  try { COMPARE(vref2, ref, "refcheck8"); } catch(wrap_error &) {}

  Vector<Real> v4w = vp2;           COMPARE_STORAGE( v4w,  vref1, "pw_constr=");
  for (UInt i = 0; i < ref_size; ++i) v4w(i) = mod[i];
  try { COMPARE(v4w, mod, "pwc_acces"); } catch (wrap_error &) {}
  COMPARE_STORAGE(v4w, vref2, "refcheck9");
  try { COMPARE(vref2, mod, "refcheck10"); } catch(wrap_error &) {}

  vref1 = v;

  Real store[ref_size] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  Vector<Real> vs(store, 10);
  VectorProxy<Real> vp3(vs);
  vp3 = vref1;
  try { COMPARE(vref1, store, "vp_equal_v"); } catch(wrap_error &) {}

  // Vector<Real> vref3(vm);
  // VectorProxy<Real> vp4 = vref3;
  // vp3 = vp4;
  // try { COMPARE(vs, mod, "vp_equal_vp"); } catch(wrap_error &) {}

  /* ------------------------------------------------------------------------ */
  std::cout << "--  Matrices: " << std::endl;

  Matrix<Real> m (ref, 5, 2);
  Matrix<Real> mt(ref, 2, 5);

  Matrix<Real> m1 (5, 2);
  Matrix<Real> m1t(2, 5);

  for (UInt i = 0; i < 5; ++i) {
    for (UInt j = 0; j < 2; ++j) {
      m1(i, j) = ref[i + j*5];
      m1t(j, i) = ref[j + i*2];
    }
  }
  COMPARE_STORAGE( m1, m, "access"  );
  COMPARE_STORAGE(m1t, m, "t_access");

  Matrix<Real> mm (mod, 5, 2);
  Matrix<Real> mmt(mod, 2, 5);

  Matrix<Real> m2(m);
  Matrix<Real> m3(m);
  for (UInt j = 0; j < 2; ++j) {
    Vector<Real> v = m2(j);
    for (UInt i = 0; i < 5; ++i)
      v(i) = mm(i, j);
  }
  COMPARE_STORAGE(m2, mm, "slicing");

  for (UInt j = 0; j < 2; ++j)
    m3(j) = mm(j);

  COMPARE_STORAGE(m3,  mm, "slic_slic");
  COMPARE(mm, mod, "refcheck");


  Real mod_1[ref_size] = { 98.7982, 72.1227, 197.815, 57.6722, 47.1088, 14.9865, 13.3171, 627.973, 33.9493, 98.3052 };

  Matrix<Real> m4 (mm);
  m4 (2,0) = 197.815;
  m4 (2,1) = 627.973;
  COMPARE(m4,  mod_1, "partial");

  Matrix<Real> m4t(mmt);
  m4t(0,1) = 197.815;
  m4t(1,3) = 627.973;
  COMPARE(m4t, mod_1, "t_partial");
}

/* -------------------------------------------------------------------------- */
void test_simple_operators() {
  std::cout << "=== Test simple operation ===" << std::endl;
  Real ref[ref_size] = { 23.1594, 79.6184, 77.9052, 47.9922, 12.8674, 37.1445, 64.8991, 80.3364, 98.4064, 73.7858 };
  Real mod[ref_size] = { 98.7982, 72.1227, 19.7815, 57.6722, 47.1088, 14.9865, 13.3171, 62.7973, 33.9493, 98.3052 };

  Real ref_div[ref_size] = { 1.163905920192984e+00, 4.001326766509196e+00,
			     3.915227661071464e+00, 2.411910744798472e+00,
			     6.466680068348578e-01, 1.866745401547894e+00,
			     3.261589104432606e+00, 4.037410795054780e+00,
			     4.945542265554328e+00, 3.708201829329581e+00 };
  Real ref_tim[ref_size] = { 4.608257412000000e+02, 1.584246923200000e+03,
			     1.550157669600000e+03, 9.549487955999999e+02,
			     2.560355252000000e+02, 7.391012610000000e+02,
			     1.291362291800000e+03, 1.598533687200000e+03,
			     1.958090547200000e+03, 1.468189848400000e+03 };
  Real ref_p_mod[ref_size] = { 1.219576000000000e+02, 1.517411000000000e+02,
			       9.768670000000000e+01, 1.056644000000000e+02,
			       5.997620000000001e+01, 5.213100000000000e+01,
			       7.821620000000000e+01, 1.431337000000000e+02,
			       1.323557000000000e+02, 1.720910000000000e+02 };
  Real ref_m_mod[ref_size] = { -7.563879999999999e+01,  7.495699999999999e+00,
 			        5.812369999999999e+01, -9.680000000000000e+00,
			       -3.424140000000000e+01,  2.215800000000000e+01,
 			        5.158200000000001e+01,  1.753910000000000e+01,
 			        6.445710000000000e+01, -2.451940000000000e+01 };
  std::cout << "--  Vectors: " << std::endl;
  Vector<Real> v (ref, ref_size);
  Vector<Real> vm(mod, ref_size);
  Vector<Real> vref(v);
  Vector<Real> vmod(vm);

  Vector<Real> v1 = vref / 19.898;       COMPARE(v1, ref_div,   "v / s"   );
  Vector<Real> v2 = vref * 19.898;       COMPARE(v2, ref_tim,   "v * s"   );
  Vector<Real> v3 = 19.898 * vref;       COMPARE(v3, ref_tim,   "s * v"   );
  Vector<Real> v4 = vref + vmod;         COMPARE(v4, ref_p_mod, "v1 + v2" );
  Vector<Real> v5 = vref - vmod;         COMPARE(v5, ref_m_mod, "v1 - v2" );
  Vector<Real> v6 = vref; v6 *= 19.898;  COMPARE(v6, ref_tim,   "v *= s"  );
  Vector<Real> v7 = vref; v7 /= 19.898;  COMPARE(v7, ref_div,   "v /= s"  );
  Vector<Real> v8 = vref; v8 += vmod;    COMPARE(v8, ref_p_mod, "v1 += v2");
  Vector<Real> v9 = vref; v9 -= vmod;    COMPARE(v9, ref_m_mod, "v1 -= v2");

  std::cout << "--  Matrices: " << std::endl;
  Matrix<Real> m (ref, 5, 2);
  Matrix<Real> mm(mod, 5, 2);
  Matrix<Real> mref(m);
  Matrix<Real> mmod(mm);

  Matrix<Real> m1 = mref / 19.898;       COMPARE(m1, ref_div,   "m / s"   );
  Matrix<Real> m2 = mref * 19.898;       COMPARE(m2, ref_tim,   "m * s"   );
  Matrix<Real> m3 = 19.898 * mref;       COMPARE(m3, ref_tim,   "s * m"   );
  Matrix<Real> m4 = mref + mmod;         COMPARE(m4, ref_p_mod, "m1 + m2" );
  Matrix<Real> m5 = mref - mmod;         COMPARE(m5, ref_m_mod, "m1 - m2" );
  Matrix<Real> m6 = mref; m6 *= 19.898;  COMPARE(m6, ref_tim,   "m *= s"  );
  Matrix<Real> m7 = mref; m7 /= 19.898;  COMPARE(m7, ref_div,   "m /= s"  );
  Matrix<Real> m8 = mref; m8 += mmod;    COMPARE(m8, ref_p_mod, "m1 += m2");
  Matrix<Real> m9 = mref; m9 -= mmod;    COMPARE(m9, ref_m_mod, "m1 -= m2");
}
// clang-format on

/* -------------------------------------------------------------------------- */
int main() {
  test_constructor();
  test_equal_and_accessors();
  test_simple_operators();

  return 0;
}
