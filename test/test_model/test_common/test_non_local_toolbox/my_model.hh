/**
 * @file   my_model.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Sep 11 2017
 * @date last modification: Sat Feb 03 2018
 *
 * @brief  A Documented file.
 *
 * @section LICENSE
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
#include "integrator_gauss.hh"
#include "model.hh"
#include "non_local_manager.hh"
#include "non_local_manager_callback.hh"
#include "non_local_neighborhood_base.hh"
#include "shape_lagrange.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

class MyModel : public Model, public NonLocalManagerCallback {
  using MyFEEngineType = FEEngineTemplate<IntegratorGauss, ShapeLagrange>;

public:
  MyModel(Mesh & mesh, UInt spatial_dimension)
      : Model(mesh, ModelType::_model, spatial_dimension),
        manager(*this, *this) {
    registerFEEngineObject<MyFEEngineType>("FEEngine", mesh, spatial_dimension);
    manager.registerNeighborhood("test_region", "test_region");

    getFEEngine().initShapeFunctions();
    manager.initialize();
  }

  void initModel() override {}

  MatrixType getMatrixType(const ID &) override { return _mt_not_defined; }
  std::tuple<ID, TimeStepSolverType>
  getDefaultSolverID(const AnalysisMethod & /*method*/) {
    return std::make_tuple("test", TimeStepSolverType::_static);
  }

  void assembleMatrix(const ID &) override {}
  void assembleLumpedMatrix(const ID &) override {}
  void assembleResidual() override {}

  void onNodesAdded(const Array<UInt> &, const NewNodesEvent &) override {}

  void onNodesRemoved(const Array<UInt> &, const Array<UInt> &,
                      const RemovedNodesEvent &) override {}
  void onElementsAdded(const Array<Element> &,
                       const NewElementsEvent &) override {}
  void onElementsRemoved(const Array<Element> &,
                         const ElementTypeMapArray<UInt> &,
                         const RemovedElementsEvent &) override {}
  void onElementsChanged(const Array<Element> &, const Array<Element> &,
                         const ElementTypeMapArray<UInt> &,
                         const ChangedElementsEvent &) override {}

  void insertIntegrationPointsInNeighborhoods(
      const GhostType & ghost_type) override {
    ElementTypeMapArray<Real> quadrature_points_coordinates(
        "quadrature_points_coordinates_tmp_nl", this->id, this->memory_id);
    quadrature_points_coordinates.initialize(this->getFEEngine(),
                                             _nb_component = spatial_dimension,
                                             _ghost_type = ghost_type);

    IntegrationPoint q;
    q.ghost_type = ghost_type;
    q.global_num = 0;

    auto & neighborhood = manager.getNeighborhood("test_region");

    for (auto & type : quadrature_points_coordinates.elementTypes(
             spatial_dimension, ghost_type)) {
      q.type = type;
      auto & quads = quadrature_points_coordinates(type, ghost_type);
      this->getFEEngine().computeIntegrationPointsCoordinates(quads, type,
                                                              ghost_type);
      auto quad_it = quads.begin(quads.getNbComponent());
      auto quad_end = quads.end(quads.getNbComponent());
      q.num_point = 0;
      for (; quad_it != quad_end; ++quad_it) {
        neighborhood.insertIntegrationPoint(q, *quad_it);
        ++q.num_point;
        ++q.global_num;
      }
    }
  }

  void computeNonLocalStresses(const GhostType &) override {}

  void updateLocalInternal(ElementTypeMapReal &, const GhostType &,
                           const ElementKind &) override {}

  void updateNonLocalInternal(ElementTypeMapReal &, const GhostType &,
                              const ElementKind &) override {}

  const auto & getNonLocalManager() const { return manager; }

private:
  NonLocalManager manager;
};
