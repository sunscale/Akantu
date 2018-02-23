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
 * @section LICENSE
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

class MultiGrainMaterialSelector : public DefaultMaterialCohesiveSelector {
public:
  MultiGrainMaterialSelector(const SolidMechanicsModelCohesive & model,
                             const ID & transgranular_id,
                             const ID & intergranular_id)
      : DefaultMaterialCohesiveSelector(model),
        transgranular_id(transgranular_id), intergranular_id(intergranular_id),
        model(model), mesh(model.getMesh()), mesh_facets(model.getMeshFacets()),
        spatial_dimension(model.getSpatialDimension()), nb_IG(0), nb_TG(0) {}

  UInt operator()(const Element & element) {
    if (mesh_facets.getSpatialDimension(element.type) ==
        (spatial_dimension - 1)) {
      const std::vector<Element> & element_to_subelement =
          mesh_facets.getElementToSubelement(element.type, element.ghost_type)(
              element.element);

      const Element & el1 = element_to_subelement[0];
      const Element & el2 = element_to_subelement[1];

      UInt grain_id1 =
          mesh.getData<UInt>("tag_0", el1.type, el1.ghost_type)(el1.element);
      if (el2 != ElementNull) {
        UInt grain_id2 =
            mesh.getData<UInt>("tag_0", el2.type, el2.ghost_type)(el2.element);
        if (grain_id1 == grain_id2) {
          // transgranular = 0 indicator
          nb_TG++;
          return model.getMaterialIndex(transgranular_id);
        } else {
          // intergranular = 1 indicator
          nb_IG++;
          return model.getMaterialIndex(intergranular_id);
        }
      } else {
        // transgranular = 0 indicator
        nb_TG++;
        return model.getMaterialIndex(transgranular_id);
      }
    } else {
      return DefaultMaterialCohesiveSelector::operator()(element);
    }
  }

private:
  ID transgranular_id, intergranular_id;
  const SolidMechanicsModelCohesive & model;
  const Mesh & mesh;
  const Mesh & mesh_facets;
  UInt spatial_dimension;

  UInt nb_IG;
  UInt nb_TG;
};

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  initialize("material.dat", argc, argv);

  const UInt spatial_dimension = 2;
  const UInt max_steps = 1000;

  Mesh mesh(spatial_dimension);
  mesh.read("square.msh");

  SolidMechanicsModelCohesive model(mesh);

  /// model initialization
  auto material_selector = std::make_shared<MultiGrainMaterialSelector>(
      model, "tg_cohesive", "ig_cohesive");

  model.setMaterialSelector(material_selector);
  model.initFull(_analysis_method = _explicit_lumped_mass,
                 _is_extrinsic = true);

  Real time_step = model.getStableTimeStep() * 0.05;
  model.setTimeStep(time_step);
  std::cout << "Time step: " << time_step << std::endl;

  model.assembleMassLumped();

  Array<Real> & position = mesh.getNodes();
  Array<Real> & velocity = model.getVelocity();
  Array<bool> & boundary = model.getBlockedDOFs();
  Array<Real> & displacement = model.getDisplacement();

  UInt nb_nodes = mesh.getNbNodes();

  /// boundary conditions
  for (UInt n = 0; n < nb_nodes; ++n) {
    if (position(n, 1) > 0.99 || position(n, 1) < -0.99)
      boundary(n, 1) = true;

    if (position(n, 0) > 0.99 || position(n, 0) < -0.99)
      boundary(n, 0) = true;
  }

  model.setBaseName("extrinsic");
  model.addDumpFieldVector("displacement");
  model.addDumpField("velocity");
  model.addDumpField("acceleration");
  model.addDumpField("internal_force");
  model.addDumpField("stress");
  model.addDumpField("grad_u");
  model.dump();

  /// initial conditions
  Real loading_rate = 0.1;
  // bar_height  = 2
  Real VI = loading_rate * 2 * 0.5;
  for (UInt n = 0; n < nb_nodes; ++n) {
    velocity(n, 1) = loading_rate * position(n, 1);
    velocity(n, 0) = loading_rate * position(n, 0);
  }

  model.dump();

  Real dispy = 0;

  /// Main loop
  for (UInt s = 1; s <= max_steps; ++s) {
    dispy += VI * time_step;
    /// update displacement on extreme nodes
    for (UInt n = 0; n < mesh.getNbNodes(); ++n) {
      if (position(n, 1) > 0.99) {
        displacement(n, 1) = dispy;
        velocity(n, 1) = VI;
      }
      if (position(n, 1) < -0.99) {
        displacement(n, 1) = -dispy;
        velocity(n, 1) = -VI;
      }
      if (position(n, 0) > 0.99) {
        displacement(n, 0) = dispy;
        velocity(n, 0) = VI;
      }
      if (position(n, 0) < -0.99) {
        displacement(n, 0) = -dispy;
        velocity(n, 0) = -VI;
      }
    }

    model.checkCohesiveStress();

    model.solveStep();

    if (s % 10 == 0) {
      model.dump();
      std::cout << "passing step " << s << "/" << max_steps << std::endl;
    }
  }

  finalize();
  return EXIT_SUCCESS;
}
