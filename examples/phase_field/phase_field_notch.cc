/* -------------------------------------------------------------------------- */
#include "non_linear_solver.hh"
#include "phase_field_model.hh"
#include "solid_mechanics_model.hh"
#include "coupler_solid_phasefield.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
#include <fstream>
#include <chrono>
/* -------------------------------------------------------------------------- */

using namespace akantu;
using clk = std::chrono::high_resolution_clock;
using second = std::chrono::duration<double>;
using millisecond = std::chrono::duration<double, std::milli>;

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
  
  auto start_time = clk::now();


  for (UInt s = 1; s < nbSteps; ++s) {

    if (s >= 500) {
      increment = 1.e-6;
    }
   
    if (s % 10 == 0 ) {
      constexpr char wheel[] = "/-\\|";
      auto elapsed = clk::now() - start_time;
      auto time_per_step = elapsed / s;
      std::cout << "\r[" << wheel[(s / 10) % 4] << "] " << std::setw(5) << s
                << "/" << nbSteps << " (" << std::setprecision(2)
                << std::fixed << std::setw(8)
                << millisecond(time_per_step).count()
                << "ms/step - elapsed: " << std::setw(8)
                << second(elapsed).count() << "s - ETA: " << std::setw(8)
                << second((nbSteps - s) * time_per_step).count() << "s)"
                << std::string(' ', 20) << std::flush;
    }
    model.applyBC(BC::Dirichlet::IncrementValue(increment, _y), "top");
      
    coupler.solve();

    if ( s % 100 == 0) { 
      model.dump();
    }
    
  }


  finalize();
  return EXIT_SUCCESS;
}
