/**
 * @file   test_solid_mechanics_model_boundary_condition.cc
 *
 * @author Dana Christen <dana.christen@epfl.ch>
 *
 * @date creation: Fri May 03 2013
 * @date last modification: Fri Sep 19 2014
 *
 * @brief  Test of the boundary condition class
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#include <iostream>
#include <sstream>
#include "aka_common.hh"
#include "aka_error.hh"
#include "mesh.hh"
#include "solid_mechanics_model.hh"
#include "boundary_condition_functor.hh"

using namespace akantu;

/* -------------------------------------------------------------------------- */

int main(int argc, char* argv[])
{
    UInt spatial_dimension(3);

    akantu::initialize("material.dat", argc, argv);

    Mesh mesh(spatial_dimension, "mesh_names");

    std::cout << "Loading the mesh." << std::endl;

    mesh.read("./cube_physical_names.msh");

    mesh.createGroupsFromMeshData<std::string>("physical_names");
    std::stringstream sstr;

    SolidMechanicsModel model(mesh);
    model.initFull();
    std::cout << model.getMaterial(0) << std::endl;

    Vector<Real> surface_traction(3);
    surface_traction(0)=0.0;
    surface_traction(1)=0.0;
    surface_traction(2)=-1.0;

    Matrix<Real> surface_stress(3, 3, 0.0);
    surface_stress(0,0)=0.0;
    surface_stress(1,1)=0.0;
    surface_stress(2,2)=-1.0;

    model.applyBC(BC::Dirichlet::FixedValue(13.0, _x), "Bottom");
    model.applyBC(BC::Dirichlet::FixedValue(13.0, _y), "Bottom");
    model.applyBC(BC::Dirichlet::FixedValue(13.0, _z), "Bottom");

    //model.applyBC(BC::Neumann::FromSameDim(surface_traction), "Top");
    model.applyBC(BC::Neumann::FromHigherDim(surface_stress), "Top");

    debug::setDebugLevel(dblTest);
    std::cout << model.getDisplacement();
    std::cout << model.getForce();
    debug::setDebugLevel(dblInfo);

    akantu::finalize();

    return EXIT_SUCCESS;
}


