/**
 * @file   test_multi_material_elastic.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Mar 03 2017
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Test with 2 elastic materials
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

#include "non_linear_solver.hh"
#include <solid_mechanics_model.hh>

using namespace akantu;

int main(int argc, char * argv[]) {
  initialize("test_multi_material_elastic.dat", argc, argv);

  UInt spatial_dimension = 2;
  Mesh mesh(spatial_dimension);

  mesh.read("test_multi_material_elastic.msh");

  SolidMechanicsModel model(mesh);

  auto && mat_sel = std::make_shared<MeshDataMaterialSelector<std::string>>(
      "physical_names", model);
  model.setMaterialSelector(mat_sel);

  model.initFull(_analysis_method = _static);

  model.applyBC(BC::Dirichlet::FlagOnly(_y), "ground");
  model.applyBC(BC::Dirichlet::FlagOnly(_x), "corner");

  Vector<Real> trac(spatial_dimension, 0.);
  trac(_y) = 1.;
  model.applyBC(BC::Neumann::FromTraction(trac), "air");

  model.addDumpField("external_force");
  model.addDumpField("internal_force");
  model.addDumpField("blocked_dofs");
  model.addDumpField("displacement");
  model.addDumpField("stress");
  model.addDumpField("grad_u");

  // model.dump();
  auto & solver = model.getNonLinearSolver("static");
  solver.set("max_iterations", 1);
  solver.set("threshold", 1e-8);
  solver.set("convergence_type", SolveConvergenceCriteria::_residual);

  model.solveStep();
  // model.dump();

  std::map<std::string, Matrix<Real>> ref_strain;
  ref_strain["strong"] = Matrix<Real>(spatial_dimension, spatial_dimension, 0.);
  ref_strain["strong"](_y, _y) = .5;

  ref_strain["weak"] = Matrix<Real>(spatial_dimension, spatial_dimension, 0.);
  ref_strain["weak"](_y, _y) = 1;

  Matrix<Real> ref_stress(spatial_dimension, spatial_dimension, 0.);
  ref_stress(_y, _y) = 1.;

  std::vector<std::string> mats = {"strong", "weak"};

  typedef Array<Real>::const_matrix_iterator mat_it;
  auto check = [](mat_it it, mat_it end, const Matrix<Real> & ref) -> bool {
    for (; it != end; ++it) {
      Real dist = (*it - ref).norm<L_2>();

      // std::cout << *it << " " << dist << " " << (dist < 1e-10 ? "OK" : "Not
      // OK") << std::endl;

      if (dist > 1e-10)
        return false;
    }

    return true;
  };

  for (auto & type : mesh.elementTypes(spatial_dimension)) {
    for (auto mat_id : mats) {
      auto & stress = model.getMaterial(mat_id).getStress(type);
      auto & grad_u = model.getMaterial(mat_id).getGradU(type);

      auto sit = stress.begin(spatial_dimension, spatial_dimension);
      auto send = stress.end(spatial_dimension, spatial_dimension);

      auto git = grad_u.begin(spatial_dimension, spatial_dimension);
      auto gend = grad_u.end(spatial_dimension, spatial_dimension);

      if (!check(sit, send, ref_stress))
        AKANTU_ERROR("The stresses are not correct");
      if (!check(git, gend, ref_strain[mat_id]))
        AKANTU_ERROR("The grad_u are not correct");
    }
  }

  finalize();

  return 0;
}
