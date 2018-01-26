/* -------------------------------------------------------------------------- */
#include "aka_iterators.hh"
#include "communicator.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
/* -------------------------------------------------------------------------- */

using namespace akantu;

class TestSynchronizerFixture : public ::testing::Test {
public:
  virtual void SetUp() {
    const UInt spatial_dimension = 3;

    mesh = std::make_unique<Mesh>(spatial_dimension);

    const auto & comm = Communicator::getStaticCommunicator();
    Int prank = comm.whoAmI();

    if (prank == 0) {
      this->mesh->read("cube.msh");
    }
  }

  virtual void TearDown() { this->mesh.reset(nullptr); }

  void initBarycenters(ElementTypeMapArray<Real> & barycenters, Mesh & mesh) {
    auto spatial_dimension = mesh.getSpatialDimension();
    barycenters.initialize(mesh, _spatial_dimension = _all_dimensions,
                           _nb_component = spatial_dimension,
                           _with_nb_element = true);

    for (auto && ghost_type : ghost_types) {
      for (const auto & type : mesh.elementTypes(_all_dimensions, ghost_type)) {
        for (auto && data : enumerate(
                 make_view(barycenters(type, ghost_type), spatial_dimension))) {
          Element element{type, UInt(std::get<0>(data)), ghost_type};
          mesh.getBarycenter(element, std::get<1>(data));
        }
      }
    }
  }

  void distribute() { this->mesh->distribute(); }

protected:
  std::unique_ptr<Mesh> mesh;
};
