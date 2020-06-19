/**
 * @file   test_dumper.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Tue Sep 02 2014
 * @date last modification: Mon Jan 22 2018
 *
 * @brief  test dumper
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "dumper_iohelper_paraview.hh"
#include "dumper_nodal_field.hh"
#include "dumper_text.hh"
#include "dumper_variable.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {
  initialize("input_file.dat", argc, argv);

  UInt spatial_dimension = 3;
  Mesh mesh(spatial_dimension);
  mesh.read("test_dumper.msh");

  SolidMechanicsModel model(mesh);
  auto && mat_selector =
      std::make_shared<MeshDataMaterialSelector<std::string>>("physical_names",
                                                              model);
  model.setMaterialSelector(mat_selector);

  model.initFull();
  model.assembleInternalForces();

  Real time_step = 0.1;

  const Array<Real> & coord = mesh.getNodes();
  Array<Real> & disp = model.getDisplacement();
  Array<bool> & bound = model.getBlockedDOFs();

  for (UInt n = 0; n < mesh.getNbNodes(); ++n) {
    Real dist = 0.;
    for (UInt d = 0; d < spatial_dimension; ++d) {
      dist += coord(n, d) * coord(n, d);
    }
    dist = sqrt(dist);

    for (UInt d = 0; d < spatial_dimension; ++d) {
      disp(n, d) = (d + 1) * dist;
      bound(n, d) = bool((n % 2) + d);
    }
  }

  // dump boundary bottom as reference
  model.setGroupDirectory("paraview", "Bottom");
  model.setGroupBaseName("paraview_bottom", "Bottom");
  model.addDumpGroupField("displacement", "Bottom");
  model.addDumpGroupField("blocked_dofs", "Bottom");

  UInt nbp = 3;
  DumperParaview prvdumper("paraview_bottom_parallel", "paraview", false);
  iohelper::Dumper & prvdpr = prvdumper.getDumper();
  for (UInt p = 0; p < nbp; ++p) {
    prvdpr.setParallelContext(p, nbp, 0);
    if (p != 0) {
      prvdumper.unRegisterField("connectivities");
      prvdumper.unRegisterField("element_type");
      prvdumper.unRegisterField("positions");
      prvdumper.unRegisterField("displacement");
    }
    prvdumper.registerFilteredMesh(
        mesh, mesh.getElementGroup("Bottom").getElements(),
        mesh.getElementGroup("Bottom").getNodeGroup().getNodes());
    prvdumper.registerField(
        "displacement",
        std::make_shared<dumpers::NodalField<Real, true>>(
            model.getDisplacement(), 0, 0,
            &(mesh.getElementGroup("Bottom").getNodeGroup().getNodes())));
    prvdumper.dump(0);
  }

  DumperText txtdumper("text_bottom", iohelper::_tdm_csv);
  txtdumper.setDirectory("paraview");
  txtdumper.setPrecision(8);
  txtdumper.setTimeStep(time_step);
  txtdumper.registerFilteredMesh(
      mesh, mesh.getElementGroup("Bottom").getElements(),
      mesh.getElementGroup("Bottom").getNodeGroup().getNodes());
  txtdumper.registerField(
      "displacement",
      std::make_shared<dumpers::NodalField<Real, true>>(
          model.getDisplacement(), 0, 0,
          &(mesh.getElementGroup("Bottom").getNodeGroup().getNodes())));
  txtdumper.registerField(
      "blocked_dofs",
      std::make_shared<dumpers::NodalField<bool, true>>(
          model.getBlockedDOFs(), 0, 0,
          &(mesh.getElementGroup("Bottom").getNodeGroup().getNodes())));

  Real pot_energy = 1.2345567891;
  Vector<Real> gforces(2, 1.);
  txtdumper.registerVariable(
      "potential_energy", std::make_shared<dumpers::Variable<Real>>(pot_energy));
  txtdumper.registerVariable(
      "global_forces",
      std::make_shared<dumpers::Variable<Vector<Real>>>(gforces));

  // dump a first time before the main loop
  model.dumpGroup("Bottom");
  txtdumper.dump();

  Real time = 0.;
  for (UInt i = 1; i < 5; ++i) {
    pot_energy += 2.;
    gforces(0) += 0.1;
    gforces(1) += 0.2;

    // pre -> cor

    // increment time after all steps of integration
    time += time_step;

    // dump after time increment
    if (i % 2 == 0) {
      txtdumper.dump(time, i);
      model.dumpGroup("Bottom");

      // parallel test
      for (UInt p = 0; p < nbp; ++p) {
        prvdpr.setParallelContext(p, nbp, 0);
        prvdumper.dump(i);
      }
    }
  }

  finalize();
  return EXIT_SUCCESS;
}
