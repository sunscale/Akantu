/**
 * @file   resolution_augmented_lagrangian.cc
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date creation: Mon Sep 15 2014
 * @date last modification: Wed Sep 17 2014
 *
 * @brief  contact resolution classes
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

//#include "resolution_augmented_lagrangian.hh"

#include "resolution_augmented_lagrangian.hh"

#define COUT(name) std::cout << std::string(#name) << ": " << name << std::endl;

__BEGIN_AKANTU__

template <int Dim>
void ContactResolution<Dim, _static, _augmented_lagrangian>::initialize() {

  // give priority to command line arguments instead of those in a file
  if (contact_argparser.has("aka_penalty"))
    options_[Epsilon] = contact_argparser["aka_penalty"];
  else
    flags_[Automatic_penalty_parameter] = true;
  if (contact_argparser.has("aka_alpha"))
    options_[Alpha] = contact_argparser["aka_alpha"];
  if (contact_argparser.has("aka_utol"))
    options_[Multiplier_tol] = contact_argparser["aka_utol"];
  if (contact_argparser.has("aka_ntol"))
    options_[Newton_tol] = contact_argparser["aka_ntol"];
  if (contact_argparser.has("aka_usteps"))
    options_[Multiplier_max_steps] = contact_argparser["aka_usteps"];
  if (contact_argparser.has("aka_nsteps"))
    options_[Newton_max_steps] = contact_argparser["aka_nsteps"];
  if (contact_argparser.has("aka_verbose"))
    flags_[Verbose] = true;
}

template <int Dim>
ContactResolution<Dim, _static, _augmented_lagrangian>::ContactResolution(
    model_type &m)
    : Parsable(_st_contact), model_(m),
      multiplier_dumper_(m.getMesh().getNbNodes(), 3),
      pressure_dumper_(m.getMesh().getNbNodes(), 3) {
  // register dumpers
  m.getMesh().addDumpFieldExternal("multipliers", multiplier_dumper_);
  m.getMesh().addDumpFieldExternal("pressure", pressure_dumper_);

  // register parameters from the file
  registerParam("penalty", options_[Epsilon], _pat_parsable,
                "Penalty parameter for Augmented-Lagrangian formulation");
  registerParam("alpha", options_[Alpha], 1., _pat_parsable,
                "Multiplier for values of the penalty parameter");
  registerParam("utol", options_[Multiplier_tol], 1e-4, _pat_parsable,
                "Tolerance used for multipliers in the Uzawa method");
  registerParam("ntol", options_[Newton_tol], 1e-4, _pat_parsable,
                "Tolerance used in the Newton-Raphson inner convergence loop");
  registerParam("usteps", options_[Multiplier_max_steps], 100., _pat_parsable,
                "Maximum number of steps allowed in the Uzawa loop");
  registerParam("nsteps", options_[Newton_max_steps], 100., _pat_parsable,
                "Maximum number of steps allowed in the Newton-Raphson loop");

  // register parameters from the command line
  contact_argparser.addArgument(
      "--aka_penalty", "Penalty parameter for Augmented-Lagrangian formulation",
      1, cppargparse::_float);
  contact_argparser.addArgument(
      "--aka_alpha", "Multiplier for values of the penalty parameter", 1,
      cppargparse::_float);
  contact_argparser.addArgument(
      "--aka_utol", "Tolerance used for multipliers in the Uzawa method", 1,
      cppargparse::_float);
  contact_argparser.addArgument(
      "--aka_ntol",
      "Tolerance used in the Newton-Raphson inner convergence loop", 1,
      cppargparse::_float);
  contact_argparser.addArgument(
      "--aka_usteps", "Maximum number of steps allowed in the Uzawa loop", 1,
      cppargparse::_float);
  contact_argparser.addArgument(
      "--aka_nsteps",
      "Maximum number of steps allowed in the Newton-Raphson loop", 1,
      cppargparse::_float);
  contact_argparser.addArgument("--aka_verbose", "Verbose output flag", 0,
                                cppargparse::_boolean);
}

template <int Dim>
void
ContactResolution<Dim, _static, _augmented_lagrangian>::solveContactStepImpl(
    SearchBase *sp, Int2Type<_generalized_newton> gn) {

  ContactResolution &cd = *this;

  model_.implicitPred();
  model_.updateResidual();

  AKANTU_DEBUG_ASSERT(model_.stiffness_matrix != NULL,
                      "You should first initialize the implicit solver and "
                      "assemble the stiffness matrix");

  //** comments that start like this will comment on the work that has to be
  //** done for the implementation of the frictional terms in the code

  UInt k = 0;
  UInt ntotal = 0;

  std::list<int> nt;

  // get global stiffness matrix and force vector
  SparseMatrix &K = model_.getStiffnessMatrix();
  Array<Real> &F = model_.getResidual();

  // get size of the whole system
  UInt original = model_.increment->getSize() * Dim;
  UInt size = original + sm_.size();

  //** the size variable at this point is computed by adding only the number
  //** of slave nodes because for each slave node there's a lagrangian
  //** multiplier for the normal contact assigned to it. In the case of
  //** frictional contact, this variable will have to account for 2 (3)
  //** multipliers for each in slave node in 2D (3D), accounting for the
  //** tangential components.

  contact_status_ = std::map<UInt, bool>();
  status_change_ = std::map<UInt, int>();

  size_t ccc = 0;
  for (auto g : gaps_) {
    if (std::abs(g.second) <= 1.e-10) {
      contact_status_[g.first] = true;
      ++ccc;
    }
  }

  Array<Real> solution(size);
  Array<Real> rhs(size);
  rhs.clear();

  // extend data structures to consider Lagrange multipliers
  K.resize(size);

  //** with the right value of the 'size' variable, all these data structures
  //** will be resized accordingly

  cout << std::boolalpha;

  if (cd[Verbose])
    cout << "- Start Generalized Newton:" << endl;

  UInt j = 0;

  bool converged = false;
  bool converged_multiplier = false;

  Real fnorm = 0;

  do {
    Real nerror = 0.;

    cd.niter_ = j;

    // assemble material matrix
    model_.assembleStiffnessMatrix();

    // copy residual to rhs
    std::copy_n(F.storage(), original, rhs.storage());

    // compute contribution to tangent matrix and residual
    Real dummy;
    computeTangentAndResidual(solution, rhs, sp, dummy, gn);

    //** The computeTangentAndResidual is where the major development for the
    //** frictional part will take place.
    //** All the terms coded in this function into the stiffness matrix and the
    //** force vector take into account only the normal contact component. The
    //** implementation of the frictional part will include terms for the
    //** tangential multipliers. The process I followed for the implementation
    //** was to code the stiffness matrix terms from the book by Laursen
    //** (Computational Contact and Impact Mechanics). I then took the terms
    //** involving the lagrangian multiplier part from the thesis by Grzegorz
    //** Pietrzak (Continuum mechanics modelling and augmented Lagrangian
    //** formulation of large deformation frictional contact problems) as
    //** these terms were missing in the Laursen book. Some little changes had
    //** to be made into these terms, as Laursen and Pietrzak use different
    //** conventions for the sign of the gap function. In Pietrzak's thesis,
    //** the residual terms are given following Eqn. 6.20 (page 146), and the
    //** stiffness terms in following equation 6.25 (page 149).
    //** I suggest to start by the implementation of the Uzawa method, as it is
    //** not required to code the tangential lagrangian multiplier terms of
    //** Pietrzak

    // solve
    model_.template solve<IntegrationScheme2ndOrder::_displacement_corrector>(
        solution, 1., true, true, rhs);

    // copy the solution of the system to increment for the primal variable only
    std::copy_n(solution.begin(), model_.increment->getSize() * Dim,
                model_.increment->storage());

    // copy the solution of the system to increment for the lagrange multiplier
    size_t m = 0;
    std::vector<Real> multiplier_check(multipliers_.size());
    vector_type v_old(multipliers_.size());
    vector_type v_new(v_old);
    for (auto &pair : multipliers_) {

      v_old[m] = pair.second;
      v_new[m] = v_old[m] + solution[original + m];

      multiplier_check[m] = pair.second;
      multipliers_[pair.first] += solution[original + m];
      multiplier_check[m] -= pair.second;
      ++m;
    }

    Real sum_multiplier = 0.;
    for (auto m : multiplier_check)
      sum_multiplier += m * m;

    //** this check basically computes the L_2 norm of the lagrange multiplier
    //** difference, so this test may change when implementing the frictional
    //** part

    Real mnorm = sqrt(sum_multiplier);
    Real abs_tol = cd[Multiplier_tol];

    if (j == 0)
      fnorm = mnorm;

    converged_multiplier =
        (mnorm <= abs_tol || mnorm <= abs_tol * abs_tol * fnorm);

    model_.implicitCorr();
    model_.updateResidual();

    converged =
        model_.template testConvergence<_scc_increment>(cd[Newton_tol], nerror);

    if (cd[Verbose]) {

      size_t w = 10;
      cout << std::setw(2) << j << ": Primal: " << std::setw(w)
           << std::setprecision(4) << std::right << nerror
           << " <= " << cd[Newton_tol] << " = " << std::left << std::setw(5)
           << (nerror < cd[Newton_tol]) << " \tDual: " << std::setw(w)
           << std::setprecision(4) << std::right << mnorm
           << " <= " << cd[Multiplier_tol] << " = "
           << (mnorm <= cd[Multiplier_tol]) << ", " << std::setw(w) << mnorm
           << " <= " << abs_tol * abs_tol * fnorm << " = "
           << (mnorm <= abs_tol * abs_tol * fnorm) << endl;
    }

    ++j;
    AKANTU_DEBUG_INFO("[" << _scc_increment << "] Convergence iteration "
                          << std::setw(std::log10(cd[Newton_max_steps])) << j
                          << ": error " << nerror << (converged ? " < " : " > ")
                          << cd[Newton_tol] << std::endl);
  } while (!(converged && converged_multiplier) && j < cd[Newton_max_steps]);

  if (j == cd[Newton_max_steps]) {
    cout << "*** ERROR *** Newton-Raphson loop did not converge within max "
            "number of iterations: " << cd[Newton_max_steps] << endl;
    exit(1);
  }

  nt.push_back(j);
  ntotal += j;

  AKANTU_DEBUG_INFO("[" << _scc_increment << "] Uzawa convergence iteration "
                        << std::setw(std::log10(cd[Newton_max_steps])) << k
                        << std::endl);

  cout << "Generalized Newton iterations: " << j << endl;

  // dump vtk files
  this->dump();
}

template <int Dim>
void
ContactResolution<Dim, _static, _augmented_lagrangian>::solveContactStepImpl(
    SearchBase *sp, Int2Type<_uzawa> uz) {

  ContactResolution &cd = *this;

  model_.implicitPred();
  model_.updateResidual();

  AKANTU_DEBUG_ASSERT(model_.stiffness_matrix != NULL,
                      "You should first initialize the implicit solver and "
                      "assemble the stiffness matrix");

  //** comments that start like this will comment on the work that has to be
  //** done for the implementation of the frictional terms in the code

  // implementation of the Uzawa method for solving contact
  bool uzawa_converged = false;
  static UInt step = 0;
  UInt k = 0;
  UInt ntotal = 0;

  std::list<int> nt;

  std::ofstream ofs;
  ofs.open("iterations.out", std::ofstream::out | std::ofstream::app);

  // initialize Lagrange multipliers
  // NOTE: It doesn't make any difference to start from the previous
  // converged solution of Lagrange multipliers
  real_map lambda_new;

  cout << std::boolalpha;

  if (cd[Verbose])
    cout << "- Start Uzawa:" << endl;

  do {
    Real uerror = 0.;

    bool converged = false;
    UInt j = 0;

    cd.uiter_ = k;

    do {
      Real nerror = 0.;

      cd.niter_ = j;

      // assemble material matrix
      model_.assembleStiffnessMatrix();

      // compute contribution to tangent matrix and residual
      uzawa_converged = computeTangentAndResidual(lambda_new, sp, uerror, uz);

      //** The computeTangentAndResidual is where the major development for the
      //** frictional part will take place.
      //** All the terms coded in this function into the stiffness matrix and
      //** the force vector take into account only the normal contact component.
      //** The implementation of the frictional part will include terms for the
      //** tangential multipliers. The process I followed for the implementation
      //** was to code the stiffness matrix terms from the book by Laursen
      //** (Computational Contact and Impact Mechanics).
      //** I suggest you start by the implementing the Uzawa method first,
      //** before jumping to the more involved implementation of the tangential
      //** lagrangian multiplier terms of Pietrzak

      // solve
      model_.template solve<IntegrationScheme2ndOrder::_displacement_corrector>(
          *model_.increment, 1., true, true);

      model_.implicitCorr();
      model_.updateResidual();

      converged = model_.template testConvergence<_scc_increment>(
          cd[Newton_tol], nerror);

      if (cd[Verbose])
        cout << "    Newton: " << j << ", " << nerror << " < " << cd[Newton_tol]
             << " = " << (nerror < cd[Newton_tol]) << endl;

      ++j;
      AKANTU_DEBUG_INFO("[" << _scc_increment << "] Convergence iteration "
                            << std::setw(std::log10(cd[Newton_max_steps])) << j
                            << ": error " << nerror
                            << (converged ? " < " : " > ") << cd[Newton_tol]
                            << std::endl);
    } while (!converged && j < cd[Newton_max_steps]);

    if (cd[Verbose])
      cout << "  Uzawa: " << k << ", " << uerror << " < " << cd[Multiplier_tol]
           << " = " << (uerror < cd[Multiplier_tol]) << endl;

    if (j == cd[Newton_max_steps]) {
      cout << "*** ERROR *** Newton-Raphson loop did not converge within max "
              "number of iterations: " << cd[Newton_max_steps] << endl;
      exit(1);
    }

    nt.push_back(j);
    ntotal += j;

    // increment uzawa loop counter
    ++k;

    AKANTU_DEBUG_INFO("[" << _scc_increment << "] Uzawa convergence iteration "
                          << std::setw(std::log10(cd[Newton_max_steps])) << k
                          << std::endl);

    // update lagrange multipliers
    cd.multipliers_ = lambda_new;
  } while (!uzawa_converged && k < cd[Multiplier_max_steps]);

  if (k == cd[Multiplier_max_steps]) {
    cout << "*** ERROR *** Uzawa loop did not converge within max number of "
            "iterations: " << cd[Multiplier_max_steps] << endl;
    exit(1);
  }

  cout << "Summary: Uzawa [" << k << "]: Newton [" << ntotal << "]:";
  for (int n : nt)
    cout << " " << n;
  cout << endl;

  ofs << std::setw(10) << ++step << std::setw(10) << k << std::setw(10)
      << ntotal << endl;
  ofs.close();

  this->dump();
}

template <int Dim>
void ContactResolution<Dim, _static, _augmented_lagrangian>::dump() {

  multiplier_dumper_.clear();
  pressure_dumper_.clear();

  for (auto v : multipliers_) {
    element_type &el = sm_[v.first];

    if (el == element_type())
      continue;

    auto n = el.normal();

    Real lambda = v.second;

    for (size_t i = 0; i < n.size(); ++i)
      multiplier_dumper_(v.first, i) = lambda * n[i];

    // dump pressures only if area is associated with node
    auto it = areas_.find(v.first);
    if (it != areas_.end())
      for (size_t i = 0; i < n.size(); ++i) {
        Real a = it->second;
        assert(a != 0.);
        pressure_dumper_(v.first, i) = lambda * n[i] / a;
      }
    else
      cout << "*** WARNING *** Zero area for slave node " << v.first << endl;
  }
  model_.dump();
}

template <int Dim>
void
ContactResolution<Dim, _static, _augmented_lagrangian>::getPenaltyValues() {
  cout << "*** INFO *** Obtaining penalty parameters automatically. ";

  const SparseMatrix &Kconst = model_.getStiffnessMatrix();

  Real ave = 0.;
  size_t k = 0;

  // loop over pairs
  for (auto it = sm_.begin(); it != sm_.end(); ++it) {
    auto slave = it->first;
    auto master = it->second;

    if (master != element_type()) {
      std::vector<UInt> conn(master.numNodes() + 1); // 1 slave (not hardcoded)

      conn[0] = slave;
      for (UInt i = 0; i < master.numNodes(); ++i)
        conn[1 + i] = master.node(i);

      // compute normal
      vector_type nu = master.normal();

      // carry out stiffness multiplication with the normal
      // the product Kij*nj would give the force for a unit displacement
      // (i.e., the stiffness needed to move the node by 1)
      matrix_type r(Kconst.getSize(), master.numNodes() + 1);

      // loop over stifness matrix dimension
      for (size_t i = 0; i < Kconst.getSize(); ++i)
        // loop over problem dimensions
        for (int j = 0; j < Dim; ++j)
          // loop over nodes considered
          for (size_t k = 0; k < master.numNodes() + 1; ++k)
            r(i, k) += Kconst(i, conn[k] + j) * nu(j);

      // get results (norm of each column in r)
      vector_type rsum(master.numNodes() + 1);

      for (size_t i = 0; i < rsum.size(); ++i)
        for (size_t j = 0; j < r.rows(); ++j)
          rsum(i) += r(j, i) * r(j, i);

      // get average value as the penalty parameter
      Real epsilon = 0.;
      for (size_t i = 0; i < rsum.size(); ++i)
        epsilon += sqrt(rsum(i));

      epsilon /= master.numNodes() + 1;
      penalty_[slave] = epsilon;

      ave += penalty_[slave];
      ++k;
    }
    // dummy master
    else {
      // carry out stiffness multiplication with the normal
      // the product Kij*nj would give the force for a unit displacement
      // (i.e., the stiffness needed to move the node by 1)
      vector_type r(Kconst.getSize());

      // loop over stifness matrix dimension
      for (size_t i = 0; i < Kconst.getSize(); ++i)
        // loop over problem dimensions
        for (int j = 0; j < Dim; ++j)
          // loop over nodes considered
          r(i) += Kconst(i, slave + j) * 1. / Dim;

      // get results (norm of each column in r)
      Real epsilon = 0;
      for (size_t i = 0; i < r.size(); ++i)
        epsilon += r(i) * r(i);

      epsilon = sqrt(epsilon);
      penalty_[slave] = epsilon;

      ave += penalty_[slave];
      ++k;
    }
  }
  cout << "Average value: " << (*this)[Alpha] * ave / k << endl;
}

template <int dim> struct TangentTraits;

template <> struct TangentTraits<2> {

  constexpr static UInt dim = 2;
  constexpr static ElementType master_type = _segment_2;
  constexpr static InterpolationType interpolation_type =
      _itp_lagrange_segment_2;

  typedef Point<dim> point_type;
  typedef array::Array<1, Real> vector_type;
  typedef array::Array<2, Real> matrix_type;
  typedef SolidMechanicsModel model_type;

  template <class element_type>
  static bool projects(const point_type &s, const element_type &master,
                       const Array<Real> &position) {

    return has_projection(s, point_type(&position(master.node(0))),
                          point_type(&position(master.node(1))));
  }

  template <class real_tuple, class element_type, class vector_type>
  static std::tuple<matrix_type, vector_type>
  computeTangentAndResidual(model_type &model, real_tuple t,
                            element_type &master, const vector_type &sh,
                            const matrix_type &dsh, const vector_type &N) {

    const Array<Real> &position = model.getCurrentPosition();

    Real gap = std::get<0>(t);
    Real s1 = std::get<1>(t);
    Real s2 = std::get<2>(t);

    // compute the point on the surface
    point_type a(&position(master.node(0)));
    point_type b(&position(master.node(1)));

    vector_type nu = master.normal();

    // compute vector T
    point_type tau = dsh(0, 0) * a + dsh(0, 1) * b;
    vector_type T(dim * (master.numNodes() + 1));

    for (UInt i = 0; i < dim; ++i) {
      T[i] = tau[i];
      for (UInt j = 0; j < master.numNodes(); ++j)
        T[(1 + j) * dim + i] = -tau[i] * sh[j];
    }

    // compute N1
    vector_type N1(dim * (master.numNodes() + 1));

    for (UInt i = 0; i < dim; ++i) {
      for (UInt j = 0; j < master.numNodes(); ++j)
        N1[(1u + j) * dim + i] = -nu[i] * dsh(0u, j);
    }

    // compute m11
    Real m11 = tau * tau;

    // compute D1
    vector_type D1 = T + gap * N1;
    D1 *= 1. / m11;

    // Note: N1bar = N1 - k11*D1, but since k11 = 0 for 2D, then
    // N1bar = N1
    vector_type &N1bar = N1;

    // stiffness matrix (only non-zero terms for 2D implementation)
    matrix_type kc = s1 * N * transpose(N);            // first term
    kc += (s2 * gap * m11) * N1bar * transpose(N1bar); // second term
    kc -= s2 * D1 * transpose(N1);                     // sixth term
    kc -= s2 * N1 * transpose(D1);                     // eight term

    // residual vector
    vector_type fc = s2 * N;

    assert(kc.rows() == fc.size());

    return std::make_tuple(kc, fc);
  }

  template <class element_type>
  static std::tuple<point_type, vector_type>
  compute_projection(const point_type &s, element_type &master) {

    Distance_minimizator<dim, master_type> dm(s, master.coordinates());
    vector_type xi(1, dm.master_coordinates()[0]);
    return std::make_tuple(dm.point(), xi);
  }
};

template <> struct TangentTraits<3> {

  constexpr static UInt dim = 3;
  constexpr static ElementType master_type = _triangle_3;
  constexpr static InterpolationType interpolation_type =
      _itp_lagrange_triangle_3;

  typedef Point<dim> point_type;
  typedef array::Array<1, Real> vector_type;
  typedef array::Array<2, Real> matrix_type;
  typedef SolidMechanicsModel model_type;

  template <class element_type>
  static bool projects(const point_type &s, const element_type &master,
                       const Array<Real> &position) {

    return point_has_projection_to_triangle(
        s, point_type(&position(master.node(0))),
        point_type(&position(master.node(1))),
        point_type(&position(master.node(2))));
  }

  template <class real_tuple, class element_type, class vector_type>
  static std::tuple<matrix_type, vector_type>
  computeTangentAndResidual(model_type &model, real_tuple t,
                            element_type &master, const vector_type &sh,
                            const matrix_type &dsh, const vector_type &N) {

    const Array<Real> &position = model.getCurrentPosition();

    Real gap = std::get<0>(t);
    Real s1 = std::get<1>(t);
    Real s2 = std::get<2>(t);
    Real s3 = std::get<3>(t);

    // compute the point on the surface
    point_type a(&position(master.node(0)));
    point_type b(&position(master.node(1)));
    point_type c(&position(master.node(2)));

    vector_type nu = master.normal();

    point_type tau1 = dsh(0, 0) * a + dsh(0, 1) * b + dsh(0, 2) * c;
    point_type tau2 = dsh(1, 0) * a + dsh(1, 1) * b + dsh(1, 2) * c;

    vector_type nucheck(3);

    Math::vectorProduct3(&tau1[0], &tau2[0], &nucheck[0]);
    Math::normalize3(&nucheck[0]);

    if ((nucheck - nu)().norm() > 1.0e-10) {
      cout << "*** ERROR *** Normal failed" << endl;
      cout << "nu1: " << nu << endl;
      cout << "nu2: " << nucheck << endl;
      exit(1);
    }

    // compute vectors T1, T2, N1, N2
    size_t vsize = dim * (master.numNodes() + 1);
    vector_type T1(vsize), T2(vsize), N1(vsize), N2(vsize);
    for (UInt i = 0; i < dim; ++i) {
      T1[i] = tau1[i];
      T2[i] = tau2[i];
      for (UInt j = 0; j < master.numNodes(); ++j) {
        T1[(1 + j) * dim + i] = -tau1[i] * sh[j];
        T2[(1 + j) * dim + i] = -tau2[i] * sh[j];
        N1[(1 + j) * dim + i] = -nu[i] * dsh(0u, j);
        N2[(1 + j) * dim + i] = -nu[i] * dsh(1u, j);
      }
    }

    // compute matrix A = m + k*g  (but kappa is zero for linear elements)
    Real A11 = tau1 * tau1;
    Real A12 = tau1 * tau2;
    Real A22 = tau2 * tau2;
    Real detA = A11 * A22 - A12 * A12;

    // compute vectors D1, D2
    vector_type D1 =
        (1 / detA) * (A22 * (T1 + gap * N1)() - A12 * (T2 + gap * N2)())();
    vector_type D2 =
        (1 / detA) * (A11 * (T2 + gap * N2)() - A12 * (T1 + gap * N1)())();

    // Note: N1bar = N1 - k12*D2, but since k12 = 0 for linear elements, then
    // N1bar = N1, N2bar = N2
    vector_type &N1bar = N1;
    vector_type &N2bar = N2;

    // stiffness matrix (only non-zero terms for 3D implementation with linear
    // elements)

    // get covariant terms (det(A) = det(inv(A))
    Real m11 = A22 / detA;
    Real m12 = -A12 / detA;
    Real m22 = A11 / detA;

    // 1st term:
    //   epsilon * Heaviside(lambda + epsilon gap) * N * N' = s1 * N * N'
    matrix_type kc = s1 * N * transpose(N);
    // 2nd term:
    //  t_N * gap * m_11 * N1_bar * N1_bar', where t_N = <lambda + epsilon*gap>
    kc += (s3 * m11) * N1bar * transpose(N1bar);
    // 3rd and 4th terms:
    //   t_N * gap * m_12 * (N1_bar * N2_bar' + N2_bar * N1_bar')
    matrix_type tmp = N1bar * transpose(N2bar);
    tmp += N2bar * transpose(N1bar);
    kc += (s3 * m12) * tmp;
    // 5th term:
    //   t_N * gap * m_22 * N2_bar * N2_bar'
    kc += (s3 * m22) * N2bar * transpose(N2bar);
    // 6th term:
    //   - t_N * D1 * N1'
    kc -= s2 * D1 * transpose(N1);
    // 7th term:
    //   - t_N * D2 * N2'
    kc -= s2 * D2 * transpose(N2);
    // 8th term:
    //   - t_N * N1 * D1'
    kc -= s2 * N1 * transpose(D1);
    // 9th term:
    //   - t_N * N2 * D2'
    kc -= s2 * N2 * transpose(D2);

    // residual vector
    vector_type fc = s2 * N;
    assert(kc.rows() == fc.size());

    return std::make_tuple(kc, fc);
  }

  //! Function template specialization for inversion of a \f$ 3 \times 3 \f$
  // matrix.
  template <class matrix_type>
  static std::pair<matrix_type, Real> invert(matrix_type &A) {
    // obtain determinant of the matrix
    Real det = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
               A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
               A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);

    // compute inverse
    matrix_type inv(3, 3, 1. / det);

    inv[0][0] *= A[1][1] * A[2][2] - A[1][2] * A[2][1];
    inv[0][1] *= A[0][2] * A[2][1] - A[0][1] * A[2][2];
    inv[0][2] *= -A[0][2] * A[1][1] + A[0][1] * A[1][2];
    inv[1][0] *= A[1][2] * A[2][0] - A[1][0] * A[2][2];
    inv[1][1] *= A[0][0] * A[2][2] - A[0][2] * A[2][0];
    inv[1][2] *= A[0][2] * A[1][0] - A[0][0] * A[1][2];
    inv[2][0] *= -A[1][1] * A[2][0] + A[1][0] * A[2][1];
    inv[2][1] *= A[0][1] * A[2][0] - A[0][0] * A[2][1];
    inv[2][2] *= -A[0][1] * A[1][0] + A[0][0] * A[1][1];

    return std::make_pair(inv, det);
  }

  template <class vector_type, class point_type>
  static vector_type invert_map(const point_type &s, const point_type &a,
                                const point_type &b, const point_type &c) {
    typedef array::Array<2, Real> matrix_type;

    // matrix for inverse
    matrix_type A = { { b[0] - a[0], c[0] - a[0], a[0] },
                      { b[1] - a[1], c[1] - a[1], a[1] },
                      { b[2] - a[2], c[2] - a[2], a[2] } };

    std::pair<matrix_type, Real> Ainv = invert(A);
    vector_type x = { s[0], s[1], s[2] };
    vector_type r1 = Ainv.first * x;

    return vector_type{ r1[0],
                        r1[1] }; // return only the first two components of r1
  }

  template <class element_type>
  static std::tuple<point_type, vector_type>
  compute_projection(const point_type &s, element_type &master) {

    auto coord = master.coordinates();

    // compute the point on the surface
    point_type a(coord[0]);
    point_type b(coord[1]);
    point_type c(coord[2]);

    point_type p = closest_point_to_triangle(s, a, b, c);
    vector_type xi = invert_map<vector_type, point_type>(p, a, b, c);

    //      Distance_minimizator<dim, TangentTraits<dim>::master_type> dm(
    //          s, master.coordinates());
    //      xi = vector_type(dim - 1);
    //      for (int i = 0; i < dim - 1; ++i)
    //        xi[i] = dm.master_coordinates()[i];
    //      point_type p = dm.point();
    return std::make_tuple(p, xi);
  }
};

template <int dim>
bool ContactResolution<dim, _static, _augmented_lagrangian>::
    computeTangentAndResidual(real_map &lambda_new, SearchBase *cp, Real &error,
                              Int2Type<_uzawa>) {

  const Array<Real> &position = model_.getCurrentPosition();
  const Real tol = (*this)[Multiplier_tol];

  // get global stiffness matrix and force vector
  SparseMatrix &K = model_.getStiffnessMatrix();
  Array<Real> &F = model_.getResidual();
  const Array<Int> &eqnum =
      model_.getDOFSynchronizer().getLocalDOFEquationNumbers();

  static bool auto_flag = true;
  if (auto_flag) {
    auto_flag = false;
    if (!(*this)[Automatic_penalty_parameter]) {
      Real epsilon = (*this)[Epsilon];
      for (auto it = sm_.begin(); it != sm_.end(); ++it)
        penalty_[it->first] = epsilon;
      cout << "*** INFO *** Uniform penalty parameter used for all slaves: "
           << epsilon << endl;
      ;
    }
    // else get penalty values automatically
    else
      getPenaltyValues();
  }

  Real lm_diff = 0;
  Real lm_max = 0;

  auto it = sm_.begin();
  while (it != sm_.end()) {
    auto slave = it->first;

    Real epsilon = (*this)[Alpha] * penalty_[slave];
    AKANTU_DEBUG_ASSERT(epsilon != 0, "Penalty value cannot be zero");

    // get slave point
    point_type s(&position(slave));

    auto master = it->second;
    bool no_master = master == element_type();

    // if node lies outside triangle
    if (no_master || !TangentTraits<dim>::projects(s, master, position)) {

      auto r = cp->search(&position(slave));

      // try to find a new master
      if (r != -1) {
        it->second = master =
            element_type(model_, TangentTraits<dim>::master_type, r);
      }
      // else remove master-slave pair from simulation
      else {
        master = element_type();

        gaps_.erase(slave);
        lambda_new.erase(slave);
        ++it;
        continue;
      }
    }

    assert(master.type == TangentTraits<dim>::master_type);

    Distance_minimizator<dim, TangentTraits<dim>::master_type> dm(
        s, master.coordinates());

    vector_type xi = vector_type(dim - 1);
    for (int i = 0; i < dim - 1; ++i)
      xi[i] = dm.master_coordinates()[i];

    point_type p = dm.point();

    // compute normal
    vector_type nu = master.normal();
    point_type nup(static_cast<const Real *>(nu.data()));

    // compute and save gap
    Real gap = -(nup * (s - p));
    gaps_[slave] = gap;

    Real lambda_hat = multipliers_[slave] + epsilon * gap;

    if (lambda_hat < 0) {
      // increase iterator
      ++it;
      // save value of lambda
      lambda_new[slave] = 0;
      continue;
    }

    Real s1 = epsilon * Heaviside(lambda_hat);
    Real s2 = Macauley(lambda_hat); // max(0,lambda_hat)
    Real s3 = s2 * gap;

    std::vector<UInt> conn(master.numNodes() + 1); // 1 slave (not hardcoded)
    conn[0] = slave;
    for (UInt i = 0; i < master.numNodes(); ++i)
      conn[1 + i] = master.node(i);

    // evaluate shape functions at slave master coordinate
    vector_type sh(master.numNodes());
    InterpolationElement<TangentTraits<dim>::interpolation_type>::computeShapes(
        xi, sh);

    // compute vector N
    vector_type N(dim * (master.numNodes() + 1));
    for (UInt i = 0; i < dim; ++i) {
      N[i] = nu[i];
      for (UInt j = 0; j < master.numNodes(); ++j)
        N[(1 + j) * dim + i] = -nu[i] * sh[j];
    }

    matrix_type dsh(dim - 1, master.numNodes());
    InterpolationElement<TangentTraits<dim>::interpolation_type>::computeDNDS(
        xi, dsh);

    // obtain contribution to stiffness matrix and force vector depending on
    // the dimension
    auto t = TangentTraits<dim>::computeTangentAndResidual(
        model_, std::make_tuple(gap, s1, s2, s3), master, sh, dsh, N);

    matrix_type &kc = std::get<0>(t);
    vector_type &fc = std::get<1>(t);

    // assemble local components into global matrix and vector
    std::vector<UInt> eq;
    for (UInt i = 0; i < conn.size(); ++i)
      for (UInt j = 0; j < dim; ++j)
        eq.push_back(eqnum(conn[i] * dim + j));

    for (UInt i = 0; i < kc.rows(); ++i) {
      F[eq[i]] += fc(i);
      for (UInt j = i; j < kc.columns(); ++j) {
        K.addToProfile(eq[i], eq[j]);
        K(eq[i], eq[j]) += kc(i, j);
      }
    }

    // update multiplier
    lambda_new[slave] = s2;

    Real lm_old = multipliers_[slave];
    lm_max += lm_old * lm_old;
    lm_old -= s2;
    lm_diff += lm_old * lm_old;

    // increase iterator
    ++it;
  }

  if (lm_max < tol) {
    error = sqrt(lm_diff);
    return sqrt(lm_diff) < tol;
  }

  error = sqrt(lm_diff / lm_max);
  return sqrt(lm_diff / lm_max) < tol;
}

template <typename T>
Point<2, T>
closest_point_to_triangle(const Point<2, T> &p, const Point<2, T> &a,
                          const Point<2, T> &b, const Point<2, T> &c) {
  return Point<2, T>();
}

template <int dim>
bool ContactResolution<dim, _static, _augmented_lagrangian>::
    computeTangentAndResidual(Array<Real> &solution, Array<Real> &F,
                              SearchBase *cp, Real &error,
                              Int2Type<_generalized_newton>) {

  const Array<Real> &position = model_.getCurrentPosition();

  // get global stiffness matrix and force vector
  SparseMatrix &K = model_.getStiffnessMatrix();
  const Array<Int> &eqnum =
      model_.getDOFSynchronizer().getLocalDOFEquationNumbers();

  static bool auto_flag = true;
  if (auto_flag) {
    auto_flag = false;
    if (!(*this)[Automatic_penalty_parameter]) {
      Real epsilon = (*this)[Epsilon];
      for (auto it = sm_.begin(); it != sm_.end(); ++it)
        penalty_[it->first] = epsilon;
      cout << "*** INFO *** Uniform penalty parameter used for all slaves: "
           << epsilon << endl;
      ;

    }
    // else get penalty values automatically
    else
      getPenaltyValues();
  }

  // size of original system
  UInt original = model_.increment->getSize() * dim;

  // multiplier count
  size_t kk = 0;

  auto it = sm_.begin();
  while (it != sm_.end()) {
    auto slave = it->first;

    Real epsilon = (*this)[Alpha] * penalty_[slave];

    if (status_change_[slave] != 0) {

      ;
      epsilon *= (status_change_[slave] + 1.);
    }

    AKANTU_DEBUG_ASSERT(epsilon != 0, "Penalty value cannot be zero");

    // get slave point
    point_type s(&position(slave));

    auto master = it->second;
    bool no_master = master == element_type();

    static std::map<UInt, bool> excluded;

    // if node lies outside triangle
    if (no_master || !TangentTraits<dim>::projects(s, master, position)) {

      auto r = cp->search(&position(slave));

      // try to find a new master
      if (r != -1) {
        it->second = master =
            element_type(model_, TangentTraits<dim>::master_type, r);
        no_master = false;
      }
      // else remove master-slave pair from simulation
      else {
        master = element_type();
        no_master = true;
        excluded[slave] = true;
      }
    }

    Real gap;
    vector_type xi;

    if (!no_master) {

      assert(master.type == TangentTraits<dim>::master_type);

      auto tuple = TangentTraits<dim>::compute_projection(s, master);
      point_type &p = std::get<0>(tuple);
      xi = std::get<1>(tuple);

      // compute normal
      vector_type nu = master.normal();
      point_type nup(static_cast<const Real *>(nu.data()));

      //      Real old_gap = gaps_[slave];

      // compute and save gap
      gap = -(nup * (s - p));
      gaps_[slave] = gap;

      // track status
      // if node in contact
      if (contact_status_[slave]) {
        if (gap < -1.e-10) {
          contact_status_[slave] = false;
          ++status_change_[slave];
          //          cout<<"["<<status_change_[slave]<<"] changing to
          // non-contact status for node "<<slave<<". Gap from "<<old_gap<<" to
          // "<<gap<<endl;
        }

      } else {
        if (gap >= -1.e-10) {
          contact_status_[slave] = true;
          ++status_change_[slave];
          //          cout<<"["<<status_change_[slave]<<"] changing to contact
          // status for node "<<slave<<". Gap from "<<old_gap<<" to
          // "<<gap<<endl;
        }
      }
    }

    Real lambda_hat = multipliers_[slave] + epsilon * gap;

    // no contact
    if (lambda_hat < 0 || excluded[slave]) {

      size_t ii = original + kk;

      // add contribution to stiffness matrix and residual vector
      F[ii] = multipliers_[slave] / epsilon;

      K.addToProfile(ii, ii);
      K(ii, ii) += -1 / epsilon;
    }
    // contact
    else {

      Real s1 = epsilon * Heaviside(lambda_hat);
      Real s2 = Macauley(lambda_hat); // max(0,lambda_hat)
      Real s3 = s2 * gap;

      std::vector<UInt> conn(master.numNodes() + 1); // 1 slave (not hardcoded)
      conn[0] = slave;
      for (UInt i = 0; i < master.numNodes(); ++i)
        conn[1 + i] = master.node(i);

      // evaluate shape functions at slave master coordinate
      vector_type nu = master.normal();
      vector_type sh(master.numNodes());
      InterpolationElement<
          TangentTraits<dim>::interpolation_type>::computeShapes(xi, sh);

      // compute vector N
      vector_type N(dim * (master.numNodes() + 1));
      for (UInt i = 0; i < dim; ++i) {
        N[i] = nu[i];
        for (UInt j = 0; j < master.numNodes(); ++j)
          N[(1 + j) * dim + i] = -nu[i] * sh[j];
      }

      matrix_type dsh(dim - 1, master.numNodes());
      InterpolationElement<TangentTraits<dim>::interpolation_type>::computeDNDS(
          xi, dsh);

      // obtain contribution to stiffness matrix and force vector depending on
      // the dimension
      auto t = TangentTraits<dim>::computeTangentAndResidual(
          model_, std::make_tuple(gap, s1, s2, s3), master, sh, dsh, N);

      matrix_type &kc = std::get<0>(t);
      vector_type &fc = std::get<1>(t);

      Array<bool> &boundary = model_.getBlockedDOFs();

      // assemble local components into global matrix and vector not taking into
      // account fixed dofs
      std::vector<UInt> eq(conn.size() * dim);
      std::vector<bool> fixed(conn.size() * dim, false);
      for (UInt i = 0; i < conn.size(); ++i)
        for (UInt j = 0; j < dim; ++j) {
          eq.at(i *dim + j) = eqnum(conn[i] * dim + j);
          fixed.at(i *dim + j) = boundary(conn[i], j);
        }

      for (UInt i = 0; i < kc.rows(); ++i) {
        // if dof is blocked, don't add terms
        if (fixed.at(i))
          continue;
        F[eq[i]] += fc(i);
        for (UInt j = i; j < kc.columns(); ++j) {
          K.addToProfile(eq[i], eq[j]);
          K(eq[i], eq[j]) += kc(i, j);
        }
      }

      // terms corresponding to lagrangian multiplier contribution
      size_t ii = original + kk;

      // assemble contribution to force vector
      F[ii] = -gap;

      // assemble contribution to stiffness matrix (only upper-triangular)
      for (UInt i = 0; i < N.size(); ++i) {
        K.addToProfile(eq[i], ii);
        K(eq[i], ii) -= N[i];
      }
    }

    // increment multiplier counter
    ++kk;

    // increase iterator
    ++it;
  }

  return true;
}

template <int Dim>
std::ostream &operator<<(
    std::ostream &os,
    const ContactResolution<Dim, _static, _augmented_lagrangian> &cr) {

  typedef typename ContactResolution<
      Dim, _static, _augmented_lagrangian>::element_type element_type;

  os << "Augmented-Lagrangian resolution type. Parameters:" << endl;
  if (cr[Automatic_penalty_parameter])
    cout << "\tpenalty = auto" << endl;
  else
    cout << "\tpenalty = " << cr[Epsilon] << endl;

  cout << "\talpha   = " << cr[Alpha] << endl;
  cout << "\tutol    = " << cr[Multiplier_tol] << endl;
  cout << "\tntol    = " << cr[Newton_tol] << endl;
  cout << "\tusteps  = " << cr[Multiplier_max_steps] << endl;
  cout << "\tnsteps  = " << cr[Newton_max_steps] << endl;
  cout << "\tverbose = " << cr[Verbose] << endl;

  cout << "\n    Slave nodes: ";
  for (auto it = cr.sm_.begin(); it != cr.sm_.end(); ++it)
    os << it->first << " ";
  os << endl;

  // loop over pairs
  cout << "\n    Slave master pairs" << endl;
  for (auto it = cr.sm_.begin(); it != cr.sm_.end(); ++it) {
    auto slave = it->first;
    auto master = it->second;
    os << "\tslave: " << slave << ", Master: ";
    if (master == element_type())
      os << "none" << endl;
    else
      os << master << endl;
  }
  return os;
}

template std::ostream &operator<<(
    std::ostream &,
    const ContactResolution<2, _static, _augmented_lagrangian> &);

template std::ostream &operator<<(
    std::ostream &,
    const ContactResolution<3, _static, _augmented_lagrangian> &);

template class ContactResolution<2, _static, _augmented_lagrangian>;
template class ContactResolution<3, _static, _augmented_lagrangian>;

__END_AKANTU__
