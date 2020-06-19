/**
 * @file   patch_test_linear_heat_transfer_fixture.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Jan 30 2018
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  HeatTransfer patch tests fixture
 *
 *
 * Copyright (©) 2016-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "heat_transfer_model.hh"
/* -------------------------------------------------------------------------- */
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

  void initModel(const AnalysisMethod & method,
                 const std::string & material_file) override {
    TestPatchTestLinear<type, HeatTransferModel>::initModel(method,
                                                            material_file);
    if (method != _static)
      this->model->setTimeStep(0.5 * this->model->getStableTimeStep());
  }

  void checkAll() {
    auto & temperature = this->model->getTemperature();
    Matrix<Real> C = this->model->get("conductivity");
    this->checkDOFs(temperature);
    this->checkGradient(this->model->getTemperatureGradient(this->type),
                        temperature);
    this->checkResults(
        [&](const Matrix<Real> & grad_T) { return C * grad_T.transpose(); },
        this->model->getKgradT(this->type), temperature);
  }
};

using htm_types = gtest_list_t<TestElementTypes>;

TYPED_TEST_SUITE(TestPatchTestHTMLinear, htm_types);

#endif /* __AKANTU_PATCH_TEST_LINEAR_HEAT_TRANSFER_FIXTURE_HH__ */
