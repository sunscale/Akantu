/* -------------------------------------------------------------------------- */
#include "non_linear_solver.hh"
#include "phase_field_model.hh"
#include "solid_mechanics_model.hh"
#include "coupler_solid_phasefield.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;
const UInt spatial_dimension = 2;

/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[]){
  
  initialize("material_notch.dat", argc, argv);

  // create mesh
  Mesh mesh(spatial_dimension);
  mesh.read("square_notch.msh");

  CouplerSolidPhaseField coupler(mesh);
  auto & model = coupler.getSolidMechanicsModel();
  auto & phase = coupler.getPhaseFieldModel();

  model.initFull(_analysis_method = _static);
  auto && mat_selector = std::make_shared<MeshDataMaterialSelector<std::string>>(
      "physical_names", model);
  model.setMaterialSelector(mat_selector);

  auto && selector = std::make_shared<MeshDataPhaseFieldSelector<std::string>>(
       "physical_names", phase);
  phase.setPhaseFieldSelector(selector);

  phase.initFull(_analysis_method = _static);

  model.applyBC(BC::Dirichlet::FixedValue(0., _y), "bottom");
  model.applyBC(BC::Dirichlet::FixedValue(0., _x), "left"); 
  
  model.setBaseName("phase_notch");
  model.addDumpField("stress");
  model.addDumpField("grad_u");
  model.addDumpFieldVector("displacement");
  model.addDumpField("damage");
  model.dump();

  UInt nbSteps = 1500;
  Real increment = 1e-5;

  for (UInt s = 0; s < nbSteps; ++s) {

    if (s >= 500) {
      increment = 1.e-6;
    }

    model.applyBC(BC::Dirichlet::IncrementValue(increment, _y), "top");

    coupler.solve();
    model.dump();
  }


  finalize();
  return EXIT_SUCCESS;
}
