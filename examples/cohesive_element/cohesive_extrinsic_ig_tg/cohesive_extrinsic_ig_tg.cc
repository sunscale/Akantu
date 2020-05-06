/**
 * @file   cohesive_extrinsic_ig_tg.cc
 *
 * @author Seyedeh Mohadeseh Taheri Mousavi <mohadeseh.taherimousavi@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Mon Jan 18 2016
 *
 * @brief  Test for considering different cohesive properties for intergranular
 * (IG) and
 * transgranular (TG) fractures in extrinsic cohesive elements
 *
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;

/* -------------------------------------------------------------------------- */
class Velocity : public BC::Dirichlet::DirichletFunctor {
public:
  explicit Velocity(SolidMechanicsModel & model, Real vel, BC::Axis ax = _x)
      : DirichletFunctor(ax), model(model), vel(vel) {
    disp = vel * model.getTimeStep();
  }

public:
  inline void operator()(UInt node, Vector<bool> & /*flags*/,
                         Vector<Real> & disp,
                         const Vector<Real> & coord) const {
    Real sign = std::signbit(coord(axis)) ? -1. : 1.;
    disp(axis) += sign * this->disp;
    model.getVelocity()(node, axis) = sign * vel;
  }

private:
  SolidMechanicsModel & model;
  Real vel, disp;
};

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  initialize("material.dat", argc, argv);

  const UInt spatial_dimension = 2;
  const UInt max_steps = 1000;

  Mesh mesh(spatial_dimension);
  mesh.read("square.msh");

  SolidMechanicsModelCohesive model(mesh);
  MaterialCohesiveRules rules{{{"top", "bottom"}, "tg_cohesive"},
                              {{"top", "top"}, "ig_cohesive"},
                              {{"bottom", "bottom"}, "ig_cohesive"}};

  /// model initialization
  auto material_selector =
      std::make_shared<MaterialCohesiveRulesSelector>(model, rules);
  auto && current_selector = model.getMaterialSelector();
  material_selector->setFallback(current_selector);
  current_selector.setFallback(
      std::make_shared<MeshDataMaterialSelector<std::string>>("physical_names",
                                                              model));
  model.setMaterialSelector(material_selector);

  model.initFull(_analysis_method = _explicit_lumped_mass,
                 _is_extrinsic = true);

  Real time_step = model.getStableTimeStep() * 0.05;
  model.setTimeStep(time_step);
  std::cout << "Time step: " << time_step << std::endl;

  model.assembleMassLumped();

  auto & position = mesh.getNodes();
  auto & velocity = model.getVelocity();

  model.applyBC(BC::Dirichlet::FlagOnly(_y), "top");
  model.applyBC(BC::Dirichlet::FlagOnly(_y), "bottom");

  model.applyBC(BC::Dirichlet::FlagOnly(_x), "left");
  model.applyBC(BC::Dirichlet::FlagOnly(_x), "right");

  model.setBaseName("extrinsic");
  model.addDumpFieldVector("displacement");
  model.addDumpField("velocity");
  model.addDumpField("acceleration");
  model.addDumpField("internal_force");
  model.addDumpField("stress");
  model.addDumpField("grad_u");
  model.addDumpField("material_index");
  model.dump();

  /// initial conditions
  Real loading_rate = 0.1;
  // bar_height  = 2
  Real VI = loading_rate * 2 * 0.5;
  for (auto && data : zip(make_view(position, spatial_dimension),
                          make_view(velocity, spatial_dimension))) {
    std::get<1>(data) = loading_rate * std::get<0>(data);
  }

  model.dump();

  Velocity vely(model, VI, _y);
  Velocity velx(model, VI, _x);

  /// Main loop
  for (UInt s = 1; s <= max_steps; ++s) {

    model.applyBC(vely, "top");
    model.applyBC(vely, "bottom");

    model.applyBC(velx, "left");
    model.applyBC(velx, "right");

    model.checkCohesiveStress();

    model.solveStep();

    if (s % 10 == 0) {
      model.dump();
      std::cout << "passing step " << s << "/" << max_steps << std::endl;
    }
  }

  return 0;
}
