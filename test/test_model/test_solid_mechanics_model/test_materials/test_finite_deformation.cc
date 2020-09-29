/* -------------------------------------------------------------------------- */
#include <solid_mechanics_model.hh>
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <type_traits>
/* -------------------------------------------------------------------------- */

using namespace akantu;


TEST(TestFiniteDeformation, NotUnit) {
  getStaticParser().parse("material_finite_deformation.dat");

  const double pi = std::atan(1) * 4;
  constexpr int dim = 3;

  Mesh mesh(dim);
  mesh.read("1_tetrahedron.msh");

  SolidMechanicsModel model(mesh);
  model.initFull(_analysis_method = _static);

#if DEBUG_TEST
  model.addDumpField("displacement");
  model.addDumpField("internal_force");
  model.addDumpField("stress");
  model.addDumpField("strain");
  model.dump();
#endif

  Matrix<Real> alpha{{0.00, 0.02, 0.03, 0.04},
                     {0.00, 0.06, 0.07, 0.08},
                     {0.00, 0.10, 0.11, 0.12}};

  auto impose_disp = [&] {
    model.getDisplacement().zero();
    for (auto data :  zip(make_view(mesh.getNodes(), dim),
                          make_view(model.getDisplacement(), dim),
                          make_view(model.getBlockedDOFs(), dim))) {
      auto & pos = std::get<0>(data);
      auto & dis = std::get<1>(data);
      auto & blocked = std::get<2>(data);

      blocked.set(true);
      
      dis += Vector<Real>(alpha(0));
      for (auto p : arange(dim)) {
        dis += Vector<Real>(alpha(1+p)) * pos(p);
      }
    }
  };

  impose_disp();
  model.solveStep();
#if DEBUG_TEST
  model.dump();
#endif

  auto stesses0 = model.getMaterial(0).getStress();
  auto displacement0 = model.getDisplacement();
  auto internal_force0 = model.getInternalForce();

  auto theta = pi / 4;
  Matrix<Real> R{{1., 0., 0.},
                 {0., std::cos(theta), -std::sin(theta)},
                 {0., std::sin(theta), std::cos(theta)}};

  impose_disp();
  for (auto data : zip(make_view(mesh.getNodes(), dim),
                       make_view(model.getDisplacement(), dim))) {
    auto & X = std::get<0>(data);
    auto & u = std::get<1>(data);
    u = R * (X + u) - X;
  }

  model.solveStep();
#if DEBUG_TEST
  model.dump();
#endif

  for (auto data : zip(make_view(mesh.getNodes(), dim),
                       make_view(model.getDisplacement(), dim),
                       make_view(displacement0, dim),
                       make_view(model.getInternalForce(), dim),
                       make_view(internal_force0, dim))) {
    auto pos = std::get<0>(data);
    Vector<Real> refdis(dim, 0.);

    refdis += Vector<Real>(alpha(0));
    for (auto p : arange(dim)) {
      refdis += Vector<Real>(alpha(1+p)) * pos(p);
    }

    auto dis = std::get<1>(data);
    auto dis0 = std::get<2>(data);

    auto err = refdis.distance(dis0);
    EXPECT_NEAR(err, 0,  1e-14);

    auto err1 = dis.distance(R * (pos + dis0) - pos);
    EXPECT_NEAR(err1, 0,  1e-14);

    auto f = std::get<3>(data);
    auto f0 = std::get<4>(data);

    auto err3 = f.distance(R * f0);
    EXPECT_NEAR(err3, 0,  1e-5);

  }
}
