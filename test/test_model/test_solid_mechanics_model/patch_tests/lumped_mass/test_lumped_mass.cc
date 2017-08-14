/**
 * @file   test_lumped_mass.cc
 *
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Mar 27 2014
 * @date last modification: Thu Oct 15 2015
 *
 * @brief  test the lumping of the mass matrix
 *
 * @section LICENSE
 *
 * Copyright  (©)  2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de Lausanne)
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
#include "solid_mechanics_model.hh"

using namespace akantu;
/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[])
{
  UInt spatial_dimension = ElementClass<TYPE>::getSpatialDimension();
  const ElementType type = TYPE;

  akantu::initialize("material.dat", argc, argv);

  Mesh mesh(spatial_dimension);

  std::stringstream filename; filename << TYPE << ".msh";
  mesh.read(filename.str());

  SolidMechanicsModel model(mesh);

  /// model initialization
  model.initFull();
  model.assembleMassLumped();

  Real rho = model.getMaterial(0).getRho();

  FEEngine & fem = model.getFEEngine();
  UInt nb_element = mesh.getNbElement(type);
  UInt nb_quadrature_points = fem.getNbIntegrationPoints(type) * nb_element;
  Array<Real> rho_on_quad(nb_quadrature_points, 1, rho, "rho_on_quad");
  Real mass = fem.integrate(rho_on_quad, type);

  Array<Real>::const_vector_iterator mass_it = model.getMass().begin(spatial_dimension);
  Array<Real>::const_vector_iterator mass_end = model.getMass().end(spatial_dimension);

  Vector<Real> sum(spatial_dimension, 0.);
  for(; mass_it != mass_end; ++mass_it) {
    sum += *mass_it;
  }

  std::cout << mass << std::endl << sum << std::endl;

  if(!(std::abs((mass - sum[0])/mass) < 2e-15)) {
    std::cerr << "total mass is not correct" <<  std::endl;
    return EXIT_FAILURE;
  }

  for(UInt i = 1; i < spatial_dimension; ++i)
    if(!(std::abs((sum[0] - sum[i])/mass) < 1e-15)) {
      std::cerr << "total mass is not correct" <<  std::endl;
      return EXIT_FAILURE;
    }



  finalize();

  return EXIT_SUCCESS;
}
