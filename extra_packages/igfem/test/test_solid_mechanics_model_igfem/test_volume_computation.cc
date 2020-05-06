/**
 * @file   test_volume_computation.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Thu Nov 26 12:20:15 2015
 *
 * @brief  test the volume computation for the different sub-materials
 *
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
#include "material_igfem_saw_tooth_damage.hh"
#include "material_iterative_stiffness_reduction.hh"
#include "solid_mechanics_model_igfem.hh"
/* -------------------------------------------------------------------------- */
using namespace akantu;

class TestMaterialSelector : public MaterialSelector {
public:
  TestMaterialSelector(SolidMechanicsModelIGFEM & model)
      : MaterialSelector(), model(model),
        spatial_dimension(model.getSpatialDimension()) {}

  UInt operator()(const Element & element) {
    if (Mesh::getKind(element.type) == _ek_igfem)
      return 2;
    else {
      /// regular elements
      const Mesh & mesh = model.getMesh();
      Vector<Real> barycenter(this->spatial_dimension);
      mesh.getBarycenter(element, barycenter);
      /// check if element belongs to ASR gel
      if (model.isInside(barycenter))
        return 1;
    }
    return 0;
  }

protected:
  SolidMechanicsModelIGFEM & model;
  UInt spatial_dimension;
};

/* -------------------------------------------------------------------------- */
/* Main                                                                       */
/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {

  Math::setTolerance(1e-13);
  debug::setDebugLevel(dblWarning);

  initialize("material_stiffness_reduction.dat", argc, argv);

  const UInt spatial_dimension = 2;
  StaticCommunicator & comm =
      akantu::StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  /// read the mesh and partion it
  Mesh mesh(spatial_dimension);
  akantu::MeshPartition * partition = NULL;

  if (prank == 0) {

    mesh.read("test_damage_transfer.msh");

    /// partition the mesh
    partition = new MeshPartitionScotch(mesh, spatial_dimension);

    partition->partitionate(psize);
  }

  /// model creation
  SolidMechanicsModelIGFEM model(mesh);
  model.initParallel(partition);
  delete partition;

  Math::setTolerance(1.e-14);
  /// intialize the geometry and set the material selector
  std::list<SK::Sphere_3> inclusions_list;
  model.registerGeometryObject(inclusions_list, "inclusion");
  Real val = 1000000000;
  Real radius_squared = (val - 0.1) * (val - 0.1);
  Vector<Real> center(spatial_dimension);
  center(0) = 0;
  center(1) = val;
  SK::Sphere_3 sphere(SK::Point_3(center(0), center(1), 0.), radius_squared);
  inclusions_list.push_back(sphere);
  TestMaterialSelector * mat_selector = new TestMaterialSelector(model);
  model.setMaterialSelector(*mat_selector);

  /// initialization of the model
  model.initFull();

  /// boundary conditions
  mesh.computeBoundingBox();
  const Vector<Real> & lowerBounds = mesh.getLowerBounds();
  const Vector<Real> & upperBounds = mesh.getUpperBounds();
  Real bottom = lowerBounds(1);
  Real top = upperBounds(1);
  Real left = lowerBounds(0);
  Real eps = std::abs((top - bottom) * 1e-6);
  const Array<Real> & pos = mesh.getNodes();
  Array<bool> & boun = model.getBlockedDOFs();
  Array<Real> & disp = model.getDisplacement();
  for (UInt n = 0; n < mesh.getNbNodes(); ++n) {
    if (std::abs(pos(n, 1) - bottom) < eps) {
      boun(n, 1) = true;
      disp(n, 1) = 0.;
    }
    if (std::abs(pos(n, 1) - top) < eps) {
      boun(n, 1) = true;
      disp(n, 1) = 1.e-3;
    }
    if (std::abs(pos(n, 0) - left) < eps) {
      boun(n, 0) = true;
      disp(n, 0) = 0.;
    }
  }

  /// add fields that should be dumped
  model.setBaseName("regular");
  model.addDumpField("material_index");
  model.addDumpFieldVector("displacement");
  ;
  model.addDumpField("stress");
  model.addDumpField("blocked_dofs");
  model.addDumpField("residual");
  model.addDumpField("grad_u");
  model.addDumpField("damage");
  model.addDumpField("partitions");
  model.addDumpField("Sc");
  model.addDumpField("force");
  model.addDumpField("equivalent_stress");
  model.addDumpField("ultimate_strain");
  model.setBaseNameToDumper("igfem elements", "igfem elements");
  model.addDumpFieldToDumper("igfem elements", "material_index");
  model.addDumpFieldVectorToDumper("igfem elements", "displacement");
  ;
  model.addDumpFieldToDumper("igfem elements", "stress");
  model.addDumpFieldToDumper("igfem elements", "blocked_dofs");
  model.addDumpFieldToDumper("igfem elements", "residual");
  model.addDumpFieldToDumper("igfem elements", "grad_u");
  model.addDumpFieldToDumper("igfem elements", "damage");
  model.addDumpFieldToDumper("igfem elements", "partitions");
  model.addDumpFieldToDumper("igfem elements", "Sc");
  model.addDumpFieldToDumper("igfem elements", "force");
  model.addDumpFieldToDumper("igfem elements", "equivalent_stress");
  model.addDumpFieldToDumper("igfem elements", "ultimate_strain");

  model.dump();
  model.dump("igfem elements");

  Real new_radius = (val - 0.1);
  model.moveInterface(new_radius);

  model.update("inclusion");

  model.dump();
  model.dump("igfem elements");

  /// get a reference to the all the materials
  const Material & standard_material_damage = model.getMaterial(0);
  const Material & standard_material_elastic = model.getMaterial(1);
  const Material & igfem_material = model.getMaterial(2);

  const ElementType standard_type = _triangle_3;
  const ElementType igfem_type = _igfem_triangle_5;

  /// compute the volume on both sides of the interface
  /// regular elements
  const Array<UInt> & material_filter_0 =
      standard_material_damage.getElementFilter(standard_type);
  const Array<UInt> & material_filter_1 =
      standard_material_elastic.getElementFilter(standard_type);
  const Array<UInt> & material_filter_2 =
      igfem_material.getElementFilter(igfem_type);

  Array<Real> Volume_0(
      material_filter_0.getSize() *
          model.getFEEngine().getNbIntegrationPoints(standard_type),
      1, 1.);
  Real volume_material_damage = model.getFEEngine().integrate(
      Volume_0, standard_type, _not_ghost, material_filter_0);

  Array<Real> Volume_1(
      material_filter_1.getSize() *
          model.getFEEngine().getNbIntegrationPoints(standard_type),
      1, 1.);
  Real volume_material_elastic = model.getFEEngine().integrate(
      Volume_1, standard_type, _not_ghost, material_filter_1);

  /// igfem elements
  const Array<UInt> & sub_mat =
      igfem_material.getInternal<UInt>("sub_material")(igfem_type, _not_ghost);
  Array<Real> sub_mat_to_real(sub_mat.getSize(), 1, 1.);
  for (UInt i = 0; i < sub_mat.getSize(); ++i)
    sub_mat_to_real(i) = Real(sub_mat(i));

  Real volume_outside = model.getFEEngine("IGFEMFEEngine")
                            .integrate(sub_mat_to_real, igfem_type, _not_ghost,
                                       material_filter_2);
  Array<Real> IGFEMVolume(sub_mat.getSize(), 1, 1.);
  Real total_igfem_volume =
      model.getFEEngine("IGFEMFEEngine")
          .integrate(IGFEMVolume, igfem_type, _not_ghost, material_filter_2);
  Real volume_inside = total_igfem_volume - volume_outside;

  Math::setTolerance(1.e-8);
  if (!Math::are_float_equal(volume_material_damage, 0.5) ||
      !Math::are_float_equal(volume_material_elastic, 0.25) ||
      !Math::are_float_equal(volume_outside, 0.1) ||
      !Math::are_float_equal(volume_inside, (0.15))) {
    std::cout << "the test failed!!!" << std::endl;
    finalize();
    return EXIT_FAILURE;
  }

  finalize();

  return EXIT_SUCCESS;
}
