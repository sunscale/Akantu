/* -------------------------------------------------------------------------- */
#include "test_cohesive_fixture.hh"
#include "aka_iterators.hh"
/* -------------------------------------------------------------------------- */

TYPED_TEST(TestSMMCFixture, ModeI) {
  auto max_steps = 100;
  auto disp_inc = 0.1;
  for(auto _ [[gnu::unused]] : arange(max_steps)) {
    this->model->applyBC(BC::Dirichlet::IncrementValue(disp_inc, _x), "loading");
    if(this->is_extrinsic)
      this->model->checkCohesiveStress();

    this->model->solveStep();
  }

  auto nb_cohesive_element = this->mesh->getNbElement(TestFixture::cohesive_type);
  auto & mesh_facets = this->mesh->getMeshFacets();
  auto facet_type = mesh_facets.getFacetType(this->cohesive_type);
  const auto & group = mesh_facets.getElementGroup("insertion").getElements(facet_type);

  std::cout << nb_cohesive_element << " " << group.size() << std::endl;

  Real sigma = this->model->getMaterial("insertion").get("sigma_c");
  Real edis = this->model->getEnergy("dissipated");
  EXPECT_NEAR(this->surface * sigma, edis, 1e-8);
}

TYPED_TEST(TestSMMCFixture, ModeII) {
  auto max_steps = 100;
  auto disp_inc = 0.1;
  for(auto _ [[gnu::unused]] : arange(max_steps)) {
    this->model->applyBC(BC::Dirichlet::IncrementValue(disp_inc, _y ), "loading");
    if(this->is_extrinsic)
      this->model->checkCohesiveStress();

    this->model->solveStep();
  }

  auto nb_cohesive_element = this->mesh->getNbElement(TestFixture::cohesive_type);
  auto & mesh_facets = this->mesh->getMeshFacets();
  auto facet_type = mesh_facets.getFacetType(this->cohesive_type);
  const auto & group = mesh_facets.getElementGroup("insertion").getElements(facet_type);

  std::cout << nb_cohesive_element << " " << group.size() << std::endl;

  Real sigma = this->model->getMaterial("insertion").get("sigma_c");
  Real edis = this->model->getEnergy("dissipated");
  EXPECT_NEAR(this->surface * sigma, edis, 1e-8);
}
