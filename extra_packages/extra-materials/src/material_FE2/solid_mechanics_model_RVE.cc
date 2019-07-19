/**
 * @file   solid_mechanics_model_RVE.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Wed Jan 13 15:32:35 2016
 *
 * @brief  Implementation of SolidMechanicsModelRVE
 *
 * @section LICENSE
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
#include "solid_mechanics_model_RVE.hh"
#include "element_group.hh"
#include "material_damage_iterative.hh"
#include "node_group.hh"
#include "non_linear_solver.hh"
#include "non_local_manager.hh"
#include "parser.hh"
#include "sparse_matrix.hh"
/* -------------------------------------------------------------------------- */
#include <string>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
SolidMechanicsModelRVE::SolidMechanicsModelRVE(Mesh & mesh,
                                               bool use_RVE_mat_selector,
                                               UInt nb_gel_pockets, UInt dim,
                                               const ID & id,
                                               const MemoryID & memory_id)
    : SolidMechanicsModel(mesh, dim, id, memory_id), volume(0.),
      use_RVE_mat_selector(use_RVE_mat_selector),
      nb_gel_pockets(nb_gel_pockets), nb_dumps(0) {
  AKANTU_DEBUG_IN();
  /// find the four corner nodes of the RVE
  findCornerNodes();

  /// remove the corner nodes from the surface node groups:
  /// This most be done because corner nodes a not periodic
  mesh.getElementGroup("top").removeNode(corner_nodes(2));
  mesh.getElementGroup("top").removeNode(corner_nodes(3));
  mesh.getElementGroup("left").removeNode(corner_nodes(3));
  mesh.getElementGroup("left").removeNode(corner_nodes(0));
  mesh.getElementGroup("bottom").removeNode(corner_nodes(1));
  mesh.getElementGroup("bottom").removeNode(corner_nodes(0));
  mesh.getElementGroup("right").removeNode(corner_nodes(2));
  mesh.getElementGroup("right").removeNode(corner_nodes(1));

  const auto & bottom = mesh.getElementGroup("bottom").getNodeGroup();
  bottom_nodes.insert(bottom.begin(), bottom.end());

  const auto & left = mesh.getElementGroup("left").getNodeGroup();
  left_nodes.insert(left.begin(), left.end());

  // /// enforce periodicity on the displacement fluctuations
  // auto surface_pair_1 = std::make_pair("top", "bottom");
  // auto surface_pair_2 = std::make_pair("right", "left");
  // SurfacePairList surface_pairs_list;
  // surface_pairs_list.push_back(surface_pair_1);
  // surface_pairs_list.push_back(surface_pair_2);
  // TODO: To Nicolas correct the PBCs
  // this->setPBC(surface_pairs_list);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SolidMechanicsModelRVE::~SolidMechanicsModelRVE() = default;

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelRVE::initFullImpl(const ModelOptions & options) {
  AKANTU_DEBUG_IN();

  auto options_cp(options);
  options_cp.analysis_method = AnalysisMethod::_static;

  SolidMechanicsModel::initFullImpl(options_cp);

  // this->initMaterials();
  auto & fem = this->getFEEngine("SolidMechanicsFEEngine");

  /// compute the volume of the RVE
  GhostType gt = _not_ghost;
  for (auto element_type :
       this->mesh.elementTypes(spatial_dimension, gt, _ek_not_defined)) {
    Array<Real> Volume(this->mesh.getNbElement(element_type) *
                           fem.getNbIntegrationPoints(element_type),
                       1, 1.);
    this->volume = fem.integrate(Volume, element_type);
  }

  std::cout << "The volume of the RVE is " << this->volume << std::endl;

  /// dumping
  std::stringstream base_name;
  base_name << this->id; // << this->memory_id - 1;
  this->setBaseName(base_name.str());
  this->addDumpFieldVector("displacement");
  this->addDumpField("stress");
  this->addDumpField("grad_u");
  this->addDumpField("eigen_grad_u");
  this->addDumpField("blocked_dofs");
  this->addDumpField("material_index");
  this->addDumpField("damage");
  this->addDumpField("Sc");
  this->addDumpField("external_force");
  this->addDumpField("equivalent_stress");
  this->addDumpField("internal_force");
  this->addDumpField("delta_T");

  this->dump();
  this->nb_dumps += 1;
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelRVE::applyBoundaryConditions(
    const Matrix<Real> & displacement_gradient) {
  AKANTU_DEBUG_IN();
  /// get the position of the nodes
  const Array<Real> & pos = mesh.getNodes();
  /// storage for the coordinates of a given node and the displacement that will
  /// be applied
  Vector<Real> x(spatial_dimension);
  Vector<Real> appl_disp(spatial_dimension);

  /// fix top right node
  UInt node = this->corner_nodes(2);
  x(0) = pos(node, 0);
  x(1) = pos(node, 1);
  appl_disp.mul<false>(displacement_gradient, x);
  (*this->blocked_dofs)(node, 0) = true;
  (*this->displacement)(node, 0) = appl_disp(0);
  (*this->blocked_dofs)(node, 1) = true;
  (*this->displacement)(node, 1) = appl_disp(1);
  // (*this->blocked_dofs)(node,0) = true; (*this->displacement)(node,0) = 0.;
  // (*this->blocked_dofs)(node,1) = true; (*this->displacement)(node,1) = 0.;

  /// apply Hx at all the other corner nodes; H: displ. gradient
  node = this->corner_nodes(0);
  x(0) = pos(node, 0);
  x(1) = pos(node, 1);
  appl_disp.mul<false>(displacement_gradient, x);
  (*this->blocked_dofs)(node, 0) = true;
  (*this->displacement)(node, 0) = appl_disp(0);
  (*this->blocked_dofs)(node, 1) = true;
  (*this->displacement)(node, 1) = appl_disp(1);

  node = this->corner_nodes(1);
  x(0) = pos(node, 0);
  x(1) = pos(node, 1);
  appl_disp.mul<false>(displacement_gradient, x);
  (*this->blocked_dofs)(node, 0) = true;
  (*this->displacement)(node, 0) = appl_disp(0);
  (*this->blocked_dofs)(node, 1) = true;
  (*this->displacement)(node, 1) = appl_disp(1);

  node = this->corner_nodes(3);
  x(0) = pos(node, 0);
  x(1) = pos(node, 1);
  appl_disp.mul<false>(displacement_gradient, x);
  (*this->blocked_dofs)(node, 0) = true;
  (*this->displacement)(node, 0) = appl_disp(0);
  (*this->blocked_dofs)(node, 1) = true;
  (*this->displacement)(node, 1) = appl_disp(1);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelRVE::applyHomogeneousTemperature(
    const Real & temperature) {

  for (UInt m = 0; m < this->getNbMaterials(); ++m) {
    Material & mat = this->getMaterial(m);

    const ElementTypeMapArray<UInt> & filter_map = mat.getElementFilter();

    // Loop over all element types
    for (auto && type : filter_map.elementTypes(spatial_dimension)) {
      const Array<UInt> & filter = filter_map(type);
      if (filter.size() == 0)
        continue;

      Array<Real> & delta_T = mat.getArray<Real>("delta_T", type);
      Array<Real>::scalar_iterator delta_T_it = delta_T.begin();
      Array<Real>::scalar_iterator delta_T_end = delta_T.end();

      for (; delta_T_it != delta_T_end; ++delta_T_it) {
        *delta_T_it = temperature;
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelRVE::findCornerNodes() {
  AKANTU_DEBUG_IN();

  // find corner nodes
  const auto & position = mesh.getNodes();
  const auto & lower_bounds = mesh.getLowerBounds();
  const auto & upper_bounds = mesh.getUpperBounds();

  AKANTU_DEBUG_ASSERT(spatial_dimension == 2, "This is 2D only!");
  corner_nodes.resize(4);
  corner_nodes.set(UInt(-1));

  for (auto && data : enumerate(make_view(position, spatial_dimension))) {
    auto node = std::get<0>(data);
    const auto & X = std::get<1>(data);

    auto distance = X.distance(lower_bounds);
    // node 1
    if (Math::are_float_equal(distance, 0)) {
      corner_nodes(0) = node;
    }
    // node 2
    else if (Math::are_float_equal(X(_x), upper_bounds(_x)) &&
             Math::are_float_equal(X(_y), lower_bounds(_y))) {
      corner_nodes(1) = node;
    }
    // node 3
    else if (Math::are_float_equal(X(_x), upper_bounds(_x)) &&
             Math::are_float_equal(X(_y), upper_bounds(_y))) {
      corner_nodes(2) = node;
    }
    // node 4
    else if (Math::are_float_equal(X(_x), lower_bounds(_x)) &&
             Math::are_float_equal(X(_y), upper_bounds(_y))) {
      corner_nodes(3) = node;
    }
  }

  for (UInt i = 0; i < corner_nodes.size(); ++i) {
    if (corner_nodes(i) == UInt(-1))
      AKANTU_ERROR("The corner node " << i + 1 << " wasn't found");
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelRVE::advanceASR(const Matrix<Real> & prestrain) {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(spatial_dimension == 2, "This is 2D only!");

  /// apply the new eigenstrain
  for (auto element_type :
       mesh.elementTypes(spatial_dimension, _not_ghost, _ek_not_defined)) {
    Array<Real> & prestrain_vect =
        const_cast<Array<Real> &>(this->getMaterial("gel").getInternal<Real>(
            "eigen_grad_u")(element_type));
    auto prestrain_it =
        prestrain_vect.begin(spatial_dimension, spatial_dimension);
    auto prestrain_end =
        prestrain_vect.end(spatial_dimension, spatial_dimension);

    for (; prestrain_it != prestrain_end; ++prestrain_it)
      (*prestrain_it) = prestrain;
  }

  /// advance the damage
  MaterialDamageIterative<2> & mat_paste =
      dynamic_cast<MaterialDamageIterative<2> &>(*this->materials[1]);
  MaterialDamageIterative<2> & mat_aggregate =
      dynamic_cast<MaterialDamageIterative<2> &>(*this->materials[0]);
  UInt nb_damaged_elements = 0;
  Real max_eq_stress_aggregate = 0;
  Real max_eq_stress_paste = 0;

  auto & solver = this->getNonLinearSolver();
  solver.set("max_iterations", 2);
  solver.set("threshold", 1e-6);
  solver.set("convergence_type", SolveConvergenceCriteria::_solution);

  do {
    this->solveStep();

    /// compute damage
    max_eq_stress_aggregate = mat_aggregate.getNormMaxEquivalentStress();
    max_eq_stress_paste = mat_paste.getNormMaxEquivalentStress();

    nb_damaged_elements = 0;
    if (max_eq_stress_aggregate > max_eq_stress_paste)
      nb_damaged_elements = mat_aggregate.updateDamage();
    else if (max_eq_stress_aggregate < max_eq_stress_paste)
      nb_damaged_elements = mat_paste.updateDamage();
    else
      nb_damaged_elements =
          (mat_paste.updateDamage() + mat_aggregate.updateDamage());

    std::cout << "the number of damaged elements is " << nb_damaged_elements
              << std::endl;
  } while (nb_damaged_elements);

  if (this->nb_dumps % 10 == 0) {
    this->dump();
  }
  this->nb_dumps += 1;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Real SolidMechanicsModelRVE::averageTensorField(UInt row_index, UInt col_index,
                                                const ID & field_type) {
  AKANTU_DEBUG_IN();
  auto & fem = this->getFEEngine("SolidMechanicsFEEngine");
  Real average = 0;

  GhostType gt = _not_ghost;
  for (auto element_type :
       mesh.elementTypes(spatial_dimension, gt, _ek_not_defined)) {
    if (field_type == "stress") {
      for (UInt m = 0; m < this->materials.size(); ++m) {
        const auto & stress_vec = this->materials[m]->getStress(element_type);
        const auto & elem_filter =
            this->materials[m]->getElementFilter(element_type);
        Array<Real> int_stress_vec(elem_filter.size(),
                                   spatial_dimension * spatial_dimension,
                                   "int_of_stress");

        fem.integrate(stress_vec, int_stress_vec,
                      spatial_dimension * spatial_dimension, element_type,
                      _not_ghost, elem_filter);

        for (UInt k = 0; k < elem_filter.size(); ++k)
          average += int_stress_vec(k, row_index * spatial_dimension +
                                           col_index); // 3 is the value
                                                       // for the yy (in
                                                       // 3D, the value is
                                                       // 4)
      }
    } else if (field_type == "strain") {
      for (UInt m = 0; m < this->materials.size(); ++m) {
        const auto & gradu_vec = this->materials[m]->getGradU(element_type);
        const auto & elem_filter =
            this->materials[m]->getElementFilter(element_type);
        Array<Real> int_gradu_vec(elem_filter.size(),
                                  spatial_dimension * spatial_dimension,
                                  "int_of_gradu");

        fem.integrate(gradu_vec, int_gradu_vec,
                      spatial_dimension * spatial_dimension, element_type,
                      _not_ghost, elem_filter);

        for (UInt k = 0; k < elem_filter.size(); ++k)
          /// averaging is done only for normal components, so stress and strain
          /// are equal
          average +=
              0.5 *
              (int_gradu_vec(k, row_index * spatial_dimension + col_index) +
               int_gradu_vec(k, col_index * spatial_dimension + row_index));
      }
    } else if (field_type == "eigen_grad_u") {
      for (UInt m = 0; m < this->materials.size(); ++m) {
        const auto & eigen_gradu_vec =
            this->materials[m]->getInternal<Real>("eigen_grad_u")(element_type);
        const auto & elem_filter =
            this->materials[m]->getElementFilter(element_type);
        Array<Real> int_eigen_gradu_vec(elem_filter.size(),
                                        spatial_dimension * spatial_dimension,
                                        "int_of_gradu");

        fem.integrate(eigen_gradu_vec, int_eigen_gradu_vec,
                      spatial_dimension * spatial_dimension, element_type,
                      _not_ghost, elem_filter);

        for (UInt k = 0; k < elem_filter.size(); ++k)
          /// averaging is done only for normal components, so stress and strain
          /// are equal
          average +=
              int_eigen_gradu_vec(k, row_index * spatial_dimension + col_index);
      }
    } else {
      AKANTU_ERROR("Averaging not implemented for this field!!!");
    }
  }

  return average / this->volume;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelRVE::homogenizeStiffness(Matrix<Real> & C_macro) {
  AKANTU_DEBUG_IN();
  const UInt dim = 2;
  AKANTU_DEBUG_ASSERT(this->spatial_dimension == dim,
                      "Is only implemented for 2D!!!");

  /// apply three independent loading states to determine C
  /// 1. eps_el = (1;0;0) 2. eps_el = (0,1,0) 3. eps_el = (0,0,0.5)

  /// clear the eigenstrain
  Matrix<Real> zero_eigengradu(dim, dim, 0.);
  GhostType gt = _not_ghost;
  for (auto element_type : mesh.elementTypes(dim, gt, _ek_not_defined)) {
    auto & prestrain_vect =
        const_cast<Array<Real> &>(this->getMaterial("gel").getInternal<Real>(
            "eigen_grad_u")(element_type));
    auto prestrain_it =
        prestrain_vect.begin(spatial_dimension, spatial_dimension);
    auto prestrain_end =
        prestrain_vect.end(spatial_dimension, spatial_dimension);

    for (; prestrain_it != prestrain_end; ++prestrain_it)
      (*prestrain_it) = zero_eigengradu;
  }

  /// storage for results of 3 different loading states
  UInt voigt_size = VoigtHelper<dim>::size;
  Matrix<Real> stresses(voigt_size, voigt_size, 0.);
  Matrix<Real> strains(voigt_size, voigt_size, 0.);
  Matrix<Real> H(dim, dim, 0.);

  /// save the damage state before filling up cracks
  // ElementTypeMapReal saved_damage("saved_damage");
  // saved_damage.initialize(getFEEngine(), _nb_component = 1, _default_value =
  // 0);
  // this->fillCracks(saved_damage);

  /// virtual test 1:
  H(0, 0) = 0.01;
  this->performVirtualTesting(H, stresses, strains, 0);

  /// virtual test 2:
  H.clear();
  H(1, 1) = 0.01;
  this->performVirtualTesting(H, stresses, strains, 1);

  /// virtual test 3:
  H.clear();
  H(0, 1) = 0.01;
  this->performVirtualTesting(H, stresses, strains, 2);

  /// drain cracks
  // this->drainCracks(saved_damage);
  /// compute effective stiffness
  Matrix<Real> eps_inverse(voigt_size, voigt_size);
  eps_inverse.inverse(strains);
  C_macro.mul<false, false>(stresses, eps_inverse);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelRVE::performVirtualTesting(const Matrix<Real> & H,
                                                   Matrix<Real> & eff_stresses,
                                                   Matrix<Real> & eff_strains,
                                                   const UInt test_no) {
  AKANTU_DEBUG_IN();
  this->applyBoundaryConditions(H);

  auto & solver = this->getNonLinearSolver();
  solver.set("max_iterations", 2);
  solver.set("threshold", 1e-6);
  solver.set("convergence_type", SolveConvergenceCriteria::_solution);
  this->solveStep();

  /// get average stress and strain
  eff_stresses(0, test_no) = this->averageTensorField(0, 0, "stress");
  eff_strains(0, test_no) = this->averageTensorField(0, 0, "strain");
  eff_stresses(1, test_no) = this->averageTensorField(1, 1, "stress");
  eff_strains(1, test_no) = this->averageTensorField(1, 1, "strain");
  eff_stresses(2, test_no) = this->averageTensorField(1, 0, "stress");
  eff_strains(2, test_no) = 2. * this->averageTensorField(1, 0, "strain");
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelRVE::homogenizeEigenGradU(
    Matrix<Real> & eigen_gradu_macro) {
  AKANTU_DEBUG_IN();
  eigen_gradu_macro(0, 0) = this->averageTensorField(0, 0, "eigen_grad_u");
  eigen_gradu_macro(1, 1) = this->averageTensorField(1, 1, "eigen_grad_u");
  eigen_gradu_macro(0, 1) = this->averageTensorField(0, 1, "eigen_grad_u");
  eigen_gradu_macro(1, 0) = this->averageTensorField(1, 0, "eigen_grad_u");
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelRVE::initMaterials() {
  AKANTU_DEBUG_IN();

  // make sure the material are instantiated
  if (!are_materials_instantiated)
    instantiateMaterials();

  if (use_RVE_mat_selector) {
    const Vector<Real> & lowerBounds = mesh.getLowerBounds();
    const Vector<Real> & upperBounds = mesh.getUpperBounds();
    Real bottom = lowerBounds(1);
    Real top = upperBounds(1);
    Real box_size = std::abs(top - bottom);
    Real eps = box_size * 1e-6;

    auto tmp = std::make_shared<GelMaterialSelector>(*this, box_size, "gel",
                                                     this->nb_gel_pockets, eps);
    tmp->setFallback(material_selector);
    material_selector = tmp;
  }

  this->assignMaterialToElements();
  // synchronize the element material arrays
  this->synchronize(SynchronizationTag::_material_id);

  for (auto & material : materials) {
    /// init internals properties
    const auto type = material->getID();
    if (type.find("material_FE2") != std::string::npos)
      continue;
    material->initMaterial();
  }

  this->synchronize(SynchronizationTag::_smm_init_mat);

  if (this->non_local_manager) {
    this->non_local_manager->initialize();
  }
  // SolidMechanicsModel::initMaterials();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelRVE::fillCracks(ElementTypeMapReal & saved_damage) {
  const auto & mat_gel = this->getMaterial("gel");
  Real E_gel = mat_gel.get("E");
  Real E_homogenized = 0.;

  for (auto && mat : materials) {
    if (mat->getName() == "gel" || mat->getName() == "FE2_mat")
      continue;

    Real E = mat->get("E");
    auto & damage = mat->getInternal<Real>("damage");

    for (auto && type : damage.elementTypes()) {
      const auto & elem_filter = mat->getElementFilter(type);
      auto nb_integration_point = getFEEngine().getNbIntegrationPoints(type);

      auto sav_dam_it =
          make_view(saved_damage(type), nb_integration_point).begin();
      for (auto && data :
           zip(elem_filter, make_view(damage(type), nb_integration_point))) {
        auto el = std::get<0>(data);
        auto & dam = std::get<1>(data);
        Vector<Real> sav_dam = sav_dam_it[el];

        sav_dam = dam;

        for (auto q : arange(dam.size())) {
          E_homogenized = (E_gel - E) * dam(q) + E;
          dam(q) = 1. - (E_homogenized / E);
        }
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModelRVE::drainCracks(
    const ElementTypeMapReal & saved_damage) {
  for (auto && mat : materials) {
    if (mat->getName() == "gel" || mat->getName() == "FE2_mat")
      continue;
    auto & damage = mat->getInternal<Real>("damage");

    for (auto && type : damage.elementTypes()) {
      const auto & elem_filter = mat->getElementFilter(type);
      auto nb_integration_point = getFEEngine().getNbIntegrationPoints(type);

      auto sav_dam_it =
          make_view(saved_damage(type), nb_integration_point).begin();
      for (auto && data :
           zip(elem_filter, make_view(damage(type), nb_integration_point))) {
        auto el = std::get<0>(data);
        auto & dam = std::get<1>(data);
        Vector<Real> sav_dam = sav_dam_it[el];

        dam = sav_dam;
      }
    }
  }
}

} // namespace akantu
