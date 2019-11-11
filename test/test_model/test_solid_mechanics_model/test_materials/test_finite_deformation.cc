/* -------------------------------------------------------------------------- */
#include <solid_mechanics_model.hh>
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <type_traits>
/* -------------------------------------------------------------------------- */

using namespace akantu;


TEST(TestFiniteDeformation, NotUnit) {
  getStaticParser().parse("material_finite_deformation.dat");

  constexpr double pi = std::atan(1) * 4;
  constexpr int dim = 3;

  Mesh mesh(dim);
  mesh.read("cube.msh");

  SolidMechanicsModel model(mesh);
  model.initFull(_analysis_method = _static);

#if DEBUG_TEST
  model.addDumpField("displacement");
  model.addDumpField("internal_force");
  model.addDumpField("stress");
  model.addDumpField("strain");
  model.dump();
#endif

  Matrix<Real> alpha{{0.01, 0.02, 0.03, 0.04},
                     {0.05, 0.06, 0.07, 0.08},
                     {0.09, 0.10, 0.11, 0.12}};

  auto impose_disp = [&] {
    model.getDisplacement().clear();

    for (auto data : filter(mesh.getElementGroup("box").getNodeGroup(),
                            zip(make_view(mesh.getNodes(), dim),
                                make_view(model.getDisplacement(), dim),
                                make_view(model.getBlockedDOFs(), dim)))) {
      auto & pos = std::get<0>(data);
      auto & dis = std::get<1>(data);
      auto & blocked = std::get<2>(data);

      blocked.set(true);
      
      dis += Vector<Real>(alpha(0));
      for (auto p : arange(dim)) {
        dis += Vector<Real>(alpha(p)) * pos(p);
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
  auto internal_force = model.getInternalForce();

  auto theta = pi / 4;
  Matrix<Real> R{{1., 0., 0.},
                 {0., std::cos(theta), -std::sin(theta)},
                 {0., std::sin(theta), std::cos(theta)}};

  impose_disp();
  for (auto data : filter(mesh.getElementGroup("box").getNodeGroup(),
                          zip(make_view(mesh.getNodes(), dim),
                              make_view(model.getDisplacement(), dim)))) {
    auto & pos = std::get<0>(data);
    auto & dis = std::get<1>(data);
    pos = R * pos;
    dis = R * dis;
  }

  model.solveStep();
#if DEBUG_TEST
  model.dump();
#endif

  for (auto data : zip(make_view(mesh.getNodes(), dim),
                       make_view(model.getDisplacement(), dim),
                       make_view(displacement0, dim))) {
    auto pos = std::get<0>(data);
    Vector<Real> refdis(dim, 0.);

    refdis += Vector<Real>(alpha(0));
    for (auto p : arange(dim)) {
      refdis += Vector<Real>(alpha(p)) * pos(p);
    }

    auto dis = std::get<1>(data);
    auto dis0 = std::get<2>(data);

    auto norm = dis0.norm();

    auto err = refdis.distance(dis0);
    EXPECT_NEAR(err/norm, 0,  1e-15);

    auto err1 = dis.distance(R * dis0);
    EXPECT_NEAR(err1/norm, 0,  1e-15);
  }
}
