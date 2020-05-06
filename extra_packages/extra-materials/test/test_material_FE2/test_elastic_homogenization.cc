/**
 * @file   test_elastic_homogenization.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Mon Jan 25 18:32:09 2016
 *
 * @brief  Test elastic homogenization of stiffness tensor
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
#include "material_elastic_linear_anisotropic.hh"
#include "solid_mechanics_model_RVE.hh"

using namespace akantu;

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {

  akantu::initialize("material_orthotropic.dat", argc, argv);

  const UInt spatial_dimension = 2;
  const ElementType element_type = _triangle_3;
  const GhostType ghost_type = _not_ghost;
  Mesh mesh(spatial_dimension);
  mesh.read("homogenized_plate.msh");
  SolidMechanicsModelRVE model(mesh, false);

  /// model initialization
  model.initFull();

  /// apply eigenstrain
  Array<Real> & prestrain_vect =
      const_cast<Array<Real> &>(model.getMaterial(0).getInternal<Real>(
          "eigen_grad_u")(element_type, ghost_type));
  Array<Real>::iterator<Matrix<Real>> prestrain_it =
      prestrain_vect.begin(spatial_dimension, spatial_dimension);
  Array<Real>::iterator<Matrix<Real>> prestrain_end =
      prestrain_vect.end(spatial_dimension, spatial_dimension);

  //(*prestrain_it)(0,0) = 0.2;
  //(*prestrain_it)(1,1) = 0.2;

  for (; prestrain_it != prestrain_end; ++prestrain_it)
    (*prestrain_it) += 1.0;

  /// storage for results of 3 different loading states
  UInt voigt_size = VoigtHelper<spatial_dimension>::size;

  MaterialElasticLinearAnisotropic<spatial_dimension> & mat =
      dynamic_cast<MaterialElasticLinearAnisotropic<spatial_dimension> &>(
          model.getMaterial(0));
  Matrix<Real> voigt_stiffness = mat.getVoigtStiffness();

  /// homogenize
  Matrix<Real> C(voigt_size, voigt_size);
  model.homogenizeStiffness(C);
  for (UInt i = 0; i < voigt_size; ++i) {
    for (UInt j = 0; j < voigt_size; ++j) {
      std::cout << "exact: " << voigt_stiffness(i, j)
                << " approximated: " << C(i, j) << std::endl;
      if (std::abs(voigt_stiffness(i, j) - C(i, j)) > 1.e-10) {
        std::cout << "The material homogenization failed" << std::endl;
        finalize();
        return EXIT_FAILURE;
      }
    }
  }

  finalize();
  return EXIT_SUCCESS;
}
