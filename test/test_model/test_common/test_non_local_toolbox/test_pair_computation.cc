/**
 * @file   test_pair_computation.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Nov 25 2015
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  test the weight computation with and without grid
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "dumper_paraview.hh"
#include "non_local_manager.hh"
#include "non_local_neighborhood.hh"
#include "solid_mechanics_model.hh"
#include "test_material_damage.hh"
/* -------------------------------------------------------------------------- */
using namespace akantu;
typedef std::vector<std::pair<IntegrationPoint, IntegrationPoint>> PairList;

/* -------------------------------------------------------------------------- */
void computePairs(SolidMechanicsModel & model, PairList * pair_list);
int main(int argc, char * argv[]) {
  akantu::initialize("material_remove_damage.dat", argc, argv);

  // some configuration variables
  const UInt spatial_dimension = 2;

  const auto & comm = Communicator::getStaticCommunicator();
  Int prank = comm.whoAmI();

  // mesh creation and read
  Mesh mesh(spatial_dimension);
  if (prank == 0) {
    mesh.read("pair_test.msh");
  }
  mesh.distribute();

  /// model creation
  SolidMechanicsModel model(mesh);

  /// creation of material selector
  auto && mat_selector =
      std::make_shared<MeshDataMaterialSelector<std::string>>("physical_names",
                                                              model);
  model.setMaterialSelector(mat_selector);

  /// model initialization changed to use our material
  model.initFull();

  /// dump material index in paraview
  model.addDumpField("material_index");
  model.dump();

  /// compute the pairs by looping over all the quadrature points
  std::array<PairList, 2> pair_list;
  computePairs(model, pair_list.data());

  const auto & pairs_mat_1_not_ghost =
      model.getNonLocalManager().getNeighborhood("mat_1").getPairLists(_not_ghost);
  const auto & pairs_mat_1_ghost =
      model.getNonLocalManager().getNeighborhood("mat_1").getPairLists(_ghost);
  const auto & pairs_mat_2_not_ghost =
      model.getNonLocalManager().getNeighborhood("mat_2").getPairLists(_not_ghost);
  const auto & pairs_mat_2_ghost =
      model.getNonLocalManager().getNeighborhood("mat_2").getPairLists(_ghost);

  /// compare the number of pairs
  UInt nb_not_ghost_pairs_grid = pairs_mat_1_not_ghost.size() + pairs_mat_2_not_ghost.size();
  UInt nb_ghost_pairs_grid = pairs_mat_1_ghost.size() + pairs_mat_2_ghost.size();
  UInt nb_not_ghost_pairs_no_grid = pair_list[0].size();
  UInt nb_ghost_pairs_no_grid = pair_list[1].size();

  if ((nb_not_ghost_pairs_grid != nb_not_ghost_pairs_no_grid) ||
      (nb_ghost_pairs_grid != nb_ghost_pairs_no_grid)) {
    std::cout << "The number of pairs is not correct: TEST FAILED!!!"
              << std::endl;
    finalize();
    return EXIT_FAILURE;
  }

  for (UInt i = 0; i < pairs_mat_1_not_ghost.size(); ++i) {
    PairList::const_iterator it = std::find(
        pair_list[0].begin(), pair_list[0].end(), (pairs_mat_1_not_ghost)[i]);
    if (it == pair_list[0].end()) {
      std::cout << "The pairs are not correct" << std::endl;
      finalize();
      return EXIT_FAILURE;
    }
  }

  for (UInt i = 0; i < pairs_mat_2_not_ghost.size(); ++i) {
    PairList::const_iterator it = std::find(
        pair_list[0].begin(), pair_list[0].end(), (pairs_mat_2_not_ghost)[i]);
    if (it == pair_list[0].end()) {
      std::cout << "The pairs are not correct" << std::endl;
      finalize();
      return EXIT_FAILURE;
    }
  }

  for (UInt i = 0; i < pairs_mat_1_ghost.size(); ++i) {
    PairList::const_iterator it = std::find(
        pair_list[1].begin(), pair_list[1].end(), (pairs_mat_1_ghost)[i]);
    if (it == pair_list[1].end()) {
      std::cout << "The pairs are not correct" << std::endl;
      finalize();
      return EXIT_FAILURE;
    }
  }

  for (UInt i = 0; i < pairs_mat_2_ghost.size(); ++i) {
    PairList::const_iterator it = std::find(
        pair_list[1].begin(), pair_list[1].end(), (pairs_mat_2_ghost)[i]);
    if (it == pair_list[1].end()) {
      std::cout << "The pairs are not correct" << std::endl;
      finalize();
      return EXIT_FAILURE;
    }
  }

  finalize();

  return 0;
}

/* -------------------------------------------------------------------------- */
void computePairs(SolidMechanicsModel & model, PairList * pair_list) {
  ElementKind kind = _ek_regular;
  Mesh & mesh = model.getMesh();
  UInt spatial_dimension = model.getSpatialDimension();
  /// compute the quadrature points
  ElementTypeMapReal quad_coords("quad_coords");
  quad_coords.initialize(mesh, _nb_component = spatial_dimension,
                         _spatial_dimension = spatial_dimension,
                         _with_nb_element = true);
  model.getFEEngine().computeIntegrationPointsCoordinates(quad_coords);

  /// loop in a n^2 way over all the quads to generate the pairs
  Real neighborhood_radius = 0.5;

  IntegrationPoint q1;
  IntegrationPoint q2;
  GhostType ghost_type_1 = _not_ghost;
  q1.ghost_type = ghost_type_1;
  Vector<Real> q1_coords(spatial_dimension);
  Vector<Real> q2_coords(spatial_dimension);

  for (auto type_1 : mesh.elementTypes(spatial_dimension, _not_ghost, kind)) {
    q1.type = type_1;
    UInt nb_elements_1 = mesh.getNbElement(type_1, ghost_type_1);
    UInt nb_quads_1 = model.getFEEngine().getNbIntegrationPoints(type_1);
    Array<Real> & quad_coords_1 = quad_coords(q1.type, q1.ghost_type);
    auto coord_it_1 = quad_coords_1.begin(spatial_dimension);
    for (UInt e_1 = 0; e_1 < nb_elements_1; ++e_1) {
      q1.element = e_1;
      UInt mat_index_1 = model.getMaterialByElement(q1.type, q1.ghost_type)
                             .begin()[q1.element];
      for (UInt q_1 = 0; q_1 < nb_quads_1; ++q_1) {
        q1.global_num = nb_quads_1 * e_1 + q_1;
        q1.num_point = q_1;
        q1_coords = coord_it_1[q1.global_num];
        /// loop over all other quads and create pairs for this given quad
        for (auto ghost_type_2 : ghost_types) {
          q2.ghost_type = ghost_type_2;

          for (auto type_2 :
               mesh.elementTypes(spatial_dimension, ghost_type_2, kind)) {

            q2.type = type_2;
            UInt nb_elements_2 = mesh.getNbElement(type_2, ghost_type_2);
            UInt nb_quads_2 =
                model.getFEEngine().getNbIntegrationPoints(type_2);
            Array<Real> & quad_coords_2 = quad_coords(q2.type, q2.ghost_type);
            auto coord_it_2 = quad_coords_2.begin(spatial_dimension);
            for (UInt e_2 = 0; e_2 < nb_elements_2; ++e_2) {
              q2.element = e_2;
              UInt mat_index_2 =
                  model.getMaterialByElement(q2.type, q2.ghost_type)
                      .begin()[q2.element];
              for (UInt q_2 = 0; q_2 < nb_quads_2; ++q_2) {
                q2.global_num = nb_quads_2 * e_2 + q_2;
                q2.num_point = q_2;
                q2_coords = coord_it_2[q2.global_num];
                Real distance = q1_coords.distance(q2_coords);
                if (mat_index_1 != mat_index_2)
                  continue;
                else if (distance <=
                             neighborhood_radius + Math::getTolerance() &&
                         (q2.ghost_type == _ghost ||
                          (q2.ghost_type == _not_ghost &&
                           q1.global_num <=
                               q2.global_num))) { // storing only half lists
                  pair_list[q2.ghost_type].push_back(std::make_pair(q1, q2));
                }
              }
            }
          }
        }
      }
    }
  }
}
