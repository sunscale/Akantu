/* -------------------------------------------------------------------------- */
#include "aka_iterators.hh"
#include "test_cohesive_fixture.hh"
/* -------------------------------------------------------------------------- */

TYPED_TEST(TestSMMCFixture, ExtrinsicModeI) {
  getStaticParser().parse("material_0.dat");
  this->is_extrinsic = true;
  this->analysis_method = _explicit_lumped_mass;

  this->testModeI();

  this->checkInsertion();

  auto & mat_co = this->model->getMaterial("insertion");
  Real G_c = mat_co.get("G_c");

  this->checkDissipated(G_c);
}

TYPED_TEST(TestSMMCFixture, ExtrinsicModeII) {
  getStaticParser().parse("material_0.dat");
  this->is_extrinsic = true;
  this->analysis_method = _explicit_lumped_mass;

  this->testModeII();

  this->checkInsertion();

  auto & mat_co = this->model->getMaterial("insertion");
  Real G_c = mat_co.get("G_c");

  this->checkDissipated(G_c);
}

TYPED_TEST(TestSMMCFixture, IntrinsicModeI) {
  getStaticParser().parse("material_1.dat");
  this->is_extrinsic = false;
  this->analysis_method = _explicit_lumped_mass;

  this->testModeI();

  this->checkInsertion();

  auto & mat_co = this->model->getMaterial("insertion");
  Real G_c = mat_co.get("G_c");

  this->checkDissipated(G_c);
}

TYPED_TEST(TestSMMCFixture, IntrinsicModeII) {
  getStaticParser().parse("material_1.dat");
  this->is_extrinsic = false;
  this->analysis_method = _explicit_lumped_mass;

  this->testModeII();

  this->checkInsertion();

  auto & mat_co = this->model->getMaterial("insertion");
  Real G_c = mat_co.get("G_c");

  this->checkDissipated(G_c);
}
