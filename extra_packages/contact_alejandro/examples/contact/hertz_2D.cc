/**
 * @file   hertz_2D.cc
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date creation: Sun Oct 19 2014
 *
 * @brief  File used to obtain 2D implicit contact results and compare with
 * Hertz theory
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

#include <chrono>

#include "contact_impl.hh"

using namespace akantu;

using std::cout;
using std::endl;
using std::setw;
using std::setprecision;

int main(int argc, char *argv[]) {
	// set dimension
	static const UInt dim = 2;

	typedef Point <dim> point_type;
	typedef BoundingBox <dim> bbox_type;
	typedef SolidMechanicsModel model_type;
  typedef Contact <dim, MasterAssignator, SelectResolution <_static, _augmented_lagrangian> >
  contact_type;

	typedef std::chrono::high_resolution_clock clock;
	typedef std::chrono::seconds seconds;

	clock::time_point t0 = clock::now();

	initialize("steel.dat", argc, argv);

	// create mesh
	Mesh mesh(dim);

	// read mesh
	mesh.read("hertz_2D.msh");

	// create model
	model_type model(mesh);
	SolidMechanicsModelOptions opt(_static);

	// initialize material
	model.initFull(opt);
  model.updateCurrentPosition();

	// create data structure that holds contact data
	contact_type cd(argc, argv, model);

	// optimal value of penalty multiplier
	cd[Alpha] = 1.e-6;
	cd[Verbose] = true;
  
  // set Paraview output resluts
	model.setBaseName("contact");
	model.addDumpFieldVector("displacement");

	// use bounding box to minimize slave-master pairs
	Real r0 = 0.5;
	Real r1 = 0.15;
	point_type c1(-r0 / 2, -r1 / 2);
	point_type c2(r0 / 2, r1 / 2);
	bbox_type bb(c1, c2);
  
	// get physical names from mesh
	Array <Real> &coords = mesh.getNodes();
	mesh.createGroupsFromMeshData <std::string>("physical_names");

	// compute areas for slave nodes that are used for the computation of contact pressures
	model.applyBC(BC::Neumann::FromHigherDim(Matrix <Real>::eye(2, 1.)), "contact_surface");
	// NOTE: the areas are computed by assigning a unit pressure to the contact surface,
	// then the magnitude of the resulting force vector at nodes gives its associated area
	Array <Real>& areas = model.getForce();

	// add slave-master pairs and store slave node areas
	ElementGroup &eg = mesh.getElementGroup("contact_surface");
	for (auto nit = eg.node_begin(); nit != eg.node_end(); ++nit) {
		// get point of slave node
		point_type n(&coords(*nit));

		// process only if within bounding box
		if (bb & n) {
			cd.addSlave(*nit);

			// compute and add area to slave node
			Real a = 0.;
			for (UInt i = 0; i < dim; ++i)
				a += pow(areas(*nit, i), 2.);
			cd.addArea(*nit, sqrt(a));
		}
	}

	// clear force vector before the start of the simulation
	areas.clear();

	// add master surface to find pairs
	cd.searchSurface("rigid");

	// output contact data info
	cout << cd;

	// apply boundary conditions
	model.applyBC(BC::Dirichlet::FixedValue(0., _x), "rigid");
	model.applyBC(BC::Dirichlet::FixedValue(0., _y), "rigid");
	model.getBlockedDOFs()(7, 0) = true;

	Real data[3][50]; // store results for printing
	Real step = 0.001; // top displacement increment
	Real Delta = 0.05; // maximum imposed displacement
	size_t k = 0;

	// loop over displacement increments
	for (Real delta = step; delta <= Delta + step; delta += step) {
		// apply displacement at the top
		model.applyBC(BC::Dirichlet::FixedValue(-delta, _y), "top");

		// solve contact step (no need to call solve on the model object)
    solveContactStep<_generalized_newton>(cd);

		data[0][k] = delta;
		data[1][k] = cd.getForce();
		data[2][k] = cd.getMaxPressure();
		++k;
	}

	// print results
	size_t w = 10;
	cout << setprecision(2);
	cout << "\n\n" << setw(w) << "Disp." << setw(w) << "Force" << setw(w) << "Max pressure" << endl;
	for (int i = 0; i < 50; ++i)
		cout << setw(w) << data[0][i] << setw(w) << data[1][i] << setw(w) << data[2][i] << endl;

	clock::time_point t1 = clock::now();
	seconds total_s = std::chrono::duration_cast <seconds>(t1 - t0);
	cout << "*** INFO *** Simulation took " << total_s.count() << " s" << endl;

	// finalize simulation
	finalize();
	return EXIT_SUCCESS;
}
