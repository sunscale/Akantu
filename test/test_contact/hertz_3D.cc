/**
 * @file   hertz_3D.cc
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date creation: Tue May 13 2014
 * @date last modification: Fri Sep 19 2014
 *
 * @brief  This file tests for the Hertz solution in 3D
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

#include "contact_impl.hh"
#include "dumpable_inline_impl.hh"

using namespace akantu;

using std::cout;
using std::endl;
using std::setw;
using std::setprecision;


int main(int argc, char *argv[]) {
	// set dimension
	static const UInt dim = 3;

	// type definitions
	typedef Point <dim> point_type;
	typedef BoundingBox <dim> bbox_type;
	typedef SolidMechanicsModel model_type;
  typedef Contact <dim, MasterAssignator, SelectResolution <_static, _augmented_lagrangian> >
  contact_type;

	initialize("steel.dat", argc, argv);

	// create mesh
	Mesh mesh(dim);

	// read mesh
	mesh.read("hertz_3D.msh");

	// create model
	model_type model(mesh);
	SolidMechanicsModelOptions opt(_static);

	// initialize material
	model.initFull(opt);

	// create data structure that holds contact data
	contact_type cd(argc, argv, model);

	// optimal value of penalty multiplier
	cd[Alpha] = 0.125;
	cd[Multiplier_tol] = 1.e-2;
	cd[Newton_tol] = 1.e-2;

  // set Paraview output resluts
	model.setBaseName("contact");
	model.addDumpFieldVector("displacement");

	// call update current position to be able to call later
	// the function to get current positions
	model.updateCurrentPosition();


	// get physical names from Gmsh file
	mesh.createGroupsFromMeshData <std::string>("physical_names");

  // set-up bounding box to include slave nodes that lie inside it
	Real l1 = 1.;
	Real l2 = 0.2;
	Real l3 = 1.;
	point_type c1(-l1 / 2, -l2 / 2, -l3 / 2);
	point_type c2(l1 / 2,  l2 / 2,  l3 / 2);
	bbox_type bb(c1, c2);
	// get areas for the nodes of the circle
	// this is done by applying a unit pressure to the contact surface elements
	model.applyBC(BC::Neumann::FromHigherDim(Matrix <Real>::eye(3, 1.)), "contact_surface");
	Array <Real>& areas = model.getForce();

	// loop over contact surface nodes and store node areas
	ElementGroup &eg = mesh.getElementGroup("contact_surface");
  Array <Real> &coords = mesh.getNodes();
  
	cout << "- Adding areas to slave nodes. " << endl;
	for (auto nit = eg.node_begin(); nit != eg.node_end(); ++nit) {
    point_type p(&coords(*nit));
    
    // ignore slave node if it doesn't lie within the bounding box
		if (!(bb & p))
      continue;
    
    cd.addSlave(*nit);
		// compute area contributing to the slave node
		Real a = 0.;
		for (UInt i = 0; i < dim; ++i)
			a += pow(areas(*nit, i), 2.);
		cd.addArea(*nit, sqrt(a));
	}

	// set force value to zero
	areas.clear();

  // output contact data info
  cout<<cd;

  
  // add master surface to find pairs
	cd.searchSurface("rigid_surface");
	// apply boundary conditions for the rigid plane
	model.applyBC(BC::Dirichlet::FixedValue(0., _x), "bottom_body");
	model.applyBC(BC::Dirichlet::FixedValue(0., _y), "bottom_body");
	model.applyBC(BC::Dirichlet::FixedValue(0., _z), "bottom_body");

	// block z-disp in extreme points of top surface
	model.getBlockedDOFs()(1, 2) = true;
	model.getBlockedDOFs()(2, 2) = true;

	// block x-disp in extreme points of top surface
	model.getBlockedDOFs()(3, 0) = true;
	model.getBlockedDOFs()(4, 0) = true;

	const size_t steps = 20;
	Real data[3][steps]; // store results for printing
	Real step = 0.001;  // top displacement increment
	size_t k = 0;

	for (Real delta = 0; delta <= step * steps; delta += step) {
		// apply displacement to the top surface of the half-sphere
		model.applyBC(BC::Dirichlet::FixedValue(-delta, _y), "top_surface");

		// solve contact step, this function also dumps Paraview files
    solveContactStep<_uzawa>(cd);

		data[0][k] = delta;
		data[1][k] = cd.getForce();
		data[2][k] = cd.getMaxPressure();
		++k;
	}

	// print results
	size_t w = 14;
	cout << setprecision(4);
	cout << endl << setw(w) << "Disp." << setw(w) << "Force" << setw(w) << "Max pressure" << endl;
	for (size_t i = 0; i < steps; ++i)
		cout << setw(w) << data[0][i] << setw(w) << data[1][i] << setw(w) << data[2][i] << endl;

	// finalize simulation
	finalize();
	return EXIT_SUCCESS;
}
