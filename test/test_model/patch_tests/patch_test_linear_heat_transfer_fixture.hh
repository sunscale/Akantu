/* -------------------------------------------------------------------------- */
#include "heat_transfer_model.hh"
#include "patch_test_linear_fixture.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_PATCH_TEST_LINEAR_HEAT_TRANSFER_FIXTURE_HH__
#define __AKANTU_PATCH_TEST_LINEAR_HEAT_TRANSFER_FIXTURE_HH__

/* -------------------------------------------------------------------------- */
template <typename type>
class TestPatchTestHTMLinear
    : public TestPatchTestLinear<type, HeatTransferModel> {
  using parent = TestPatchTestLinear<type, HeatTransferModel>;

public:
  void applyBC() override {
    parent::applyBC();
    auto & temperature = this->model->getTemperature();
    this->applyBConDOFs(temperature);
  }

  void checkAll() {
    auto & temperature = this->model->getTemperature();
    Matrix<Real> C = this->model->get("conductivity");
    this->checkDOFs(temperature);
    this->checkGradient(this->model->getTemperatureGradient(this->type),
                        temperature);
    this->checkResults([&](const Matrix<Real> & grad_T) { return C * grad_T.transpose(); },
                       this->model->getKgradT(this->type), temperature);
  }
};

using types = gtest_list_t<TestElementTypes>;

TYPED_TEST_CASE(TestPatchTestHTMLinear, types);

#endif /* __AKANTU_PATCH_TEST_LINEAR_HEAT_TRANSFER_FIXTURE_HH__ */
