/**
 * @file   test_solid_mechanics_model_potential_energy.cc
 *
 * @author Tobias Brink <tobias.brink@epfl.ch>
 *
 * @date creation: Tue Nov 14 2017
 * @date last modification: Tue Nov 14 2017
 *
 * @brief  test potential energy of the linear elasticity model
 *
 * @section description                                                     
 *                                                                          
 * This test uses a linear elastic material with density = 1, Young's
 * modulus = 1, and Poisson's ratio = 0 and applies a linear
 * displacement from 0 to ε in x direction. The resulting potential
 * energy should be 0.5*Y*ε² = ε²/2. We test 3 different strains.
 *
 * @section LICENSE
 *
 * Copyright (©)  2017 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include <vector>

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {
  initialize("test_solid_mechanics_model_linear_elastic_potential_energy_material.dat", argc, argv);     

  UInt spatial_dimension = 2;     

  Mesh mesh(spatial_dimension);
  mesh.read("patch_tests/data/_triangle_3.msh");

  SolidMechanicsModel model(mesh);

  /// model initialization
  model.initFull(_analysis_method = _static); // _static, _explicit_lumped, _implicit_lumped

  auto & lower  = mesh.getLowerBounds();
  auto & upper  = mesh.getUpperBounds();
  auto length = upper(_x) - lower(_x);

  //std::cout << model.getMaterial(0) << std::endl;

  const Array<Real> & pos = mesh.getNodes();
  Array<Real> & disp = model.getDisplacement();
  Array<bool> & boun = model.getBlockedDOFs();

  std::vector<double> strains = {0.0, 0.1, 0.2, 0.3};
  for(auto && eps : strains) {

    /// boundary conditions
    for(auto && pair : zip(make_view(pos, spatial_dimension),
                           make_view(disp, spatial_dimension),
                           make_view(boun, spatial_dimension))) {
      auto & posv = std::get<0>(pair);
      auto & dispv = std::get<1>(pair);
      auto & bounv = std::get<2>(pair);
      auto reduced_x = (posv(_x) - lower(_x)) / length;
      dispv(_x) = reduced_x * eps;
      bounv(_x) = true;
    }

    /// "solve" a step (solution is imposed)
    model.solveStep();

    /// compare energy to analytical solution
    const Real E_ref = 0.5 * eps*eps;
    auto E_pot = model.getEnergy("potential");

    if (std::abs(E_ref - E_pot) > 1e-8) {
      std::cout << "FAIL for strain " << eps
                << "      reference energy: " << E_ref
                << "      calculated energy: " << E_pot << std::endl;
      return 1; // failure
    }
  }

  finalize();

  return EXIT_SUCCESS;
}
