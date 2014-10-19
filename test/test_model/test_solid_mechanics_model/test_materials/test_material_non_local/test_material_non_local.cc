/**
 * @file   test_material_non_local.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Aug 31 11:09:48 2011
 *
 * @brief  test of the main part of the non local materials
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

#include <iostream>
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.hh"


static void paraviewInit(iohelper::Dumper & dumper, const SolidMechanicsModel & model);
//static void paraviewDump(iohelper::Dumper & dumper);
#endif

ElementTypeMapArray<Real> quadrature_points_volumes("quadrature_points_volumes", "test");
const ElementType TYPE = _triangle_6;

int main(int argc, char *argv[]) {
  akantu::initialize("material_non_local.dat", argc, argv);
  debug::setDebugLevel(akantu::dblWarning);

  const UInt spatial_dimension = 2;
  Mesh mesh(spatial_dimension);
  MeshIOMSH mesh_io;
  mesh_io.read("mesh.msh", mesh);

  SolidMechanicsModel model(mesh);
  model.initFull();

  //  model.getFEEngine().getMesh().initElementTypeMapArray(quadrature_points_volumes, 1, 0);
  const MaterialNonLocal<spatial_dimension, BaseWeightFunction> & mat =
    dynamic_cast<const MaterialNonLocal<spatial_dimension, BaseWeightFunction> &>(model.getMaterial(0));
  //  mat.computeQuadraturePointsNeighborhoudVolumes(quadrature_points_volumes);
  Real radius = mat.getRadius();

  UInt nb_element  = mesh.getNbElement(TYPE);
  UInt nb_tot_quad = model.getFEEngine().getNbQuadraturePoints(TYPE) * nb_element;

  std::cout << mat << std::endl;

  Array<Real> quads(0, spatial_dimension);
  quads.resize(nb_tot_quad);

  model.getFEEngine().interpolateOnQuadraturePoints(mesh.getNodes(),
					       quads, spatial_dimension,
					       TYPE);

  Array<Real>::iterator< Vector<Real> > first_quad_1 = quads.begin(spatial_dimension);
  Array<Real>::iterator< Vector<Real> > last_quad_1 = quads.end(spatial_dimension);

  std::ofstream pout;
  pout.open("bf_pairs");
  UInt q1 = 0;

  Real R = mat.getRadius();

  for(;first_quad_1 != last_quad_1; ++first_quad_1, ++q1) {
    Array<Real>::iterator< Vector<Real> > first_quad_2 = quads.begin(spatial_dimension);
    //Array<Real>::iterator< Vector<Real> > last_quad_2 = quads.end(spatial_dimension);
    UInt q2 = 0;
    for(;first_quad_2 != last_quad_1; ++first_quad_2, ++q2) {
      Real d = first_quad_2->distance(*first_quad_1);
      if(d <= radius) {
	Real alpha = (1 - d*d/(R*R));
	alpha = alpha*alpha;
	pout << q1 << " " << q2 << " " << alpha << std::endl;
      }
    }
  }
  pout.close();

  mat.savePairs("cl_pairs");

  ElementTypeMapArray<Real> constant("constant_value", "test");
  mesh.initElementTypeMapArray(constant, 1, 0);
  Mesh::type_iterator it = mesh.firstType(spatial_dimension);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension);
  for(; it != last_type; ++it) {
    UInt nb_quadrature_points = model.getFEEngine().getNbQuadraturePoints(*it);
    UInt _nb_element = mesh.getNbElement(*it);

    Array<Real> & constant_vect = constant(*it);
    constant_vect.resize(_nb_element * nb_quadrature_points);

    std::fill_n(constant_vect.storage(), nb_quadrature_points * _nb_element, 1.);
  }

  ElementTypeMapArray<Real> constant_avg("constant_value_avg", "test");
  mesh.initElementTypeMapArray(constant_avg, 1, 0);

  mat.weightedAvergageOnNeighbours(constant, constant_avg, 1);

  debug::setDebugLevel(akantu::dblTest);
  std::cout << constant(TYPE) << std::endl;
  std::cout << constant_avg(TYPE) << std::endl;
  debug::setDebugLevel(akantu::dblWarning);

#ifdef AKANTU_USE_IOHELPER
  iohelper::DumperParaview dumper;
  paraviewInit(dumper, model);
#endif

  akantu::finalize();

  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
/* iohelper::Dumper vars                                                      */
/* -------------------------------------------------------------------------- */

#ifdef AKANTU_USE_IOHELPER

/* -------------------------------------------------------------------------- */
template <ElementType type> static iohelper::ElemType paraviewType();
template <> iohelper::ElemType paraviewType<_segment_2>()      { return iohelper::LINE1;     }
template <> iohelper::ElemType paraviewType<_segment_3>()      { return iohelper::LINE2;     }
template <> iohelper::ElemType paraviewType<_triangle_3>()     { return iohelper::TRIANGLE1; }
template <> iohelper::ElemType paraviewType<_triangle_6>()     { return iohelper::TRIANGLE2; }
template <> iohelper::ElemType paraviewType<_quadrangle_4>()   { return iohelper::QUAD1;     }
template <> iohelper::ElemType paraviewType<_tetrahedron_4>()  { return iohelper::TETRA1;    }
template <> iohelper::ElemType paraviewType<_tetrahedron_10>() { return iohelper::TETRA2;    }
template <> iohelper::ElemType paraviewType<_hexahedron_8>()   { return iohelper::HEX1;      }

/* -------------------------------------------------------------------------- */
void paraviewInit(iohelper::Dumper & dumper, const SolidMechanicsModel & model) {
  UInt spatial_dimension = ElementClass<TYPE>::getSpatialDimension();
  UInt nb_nodes   = model.getFEEngine().getMesh().getNbNodes();
  UInt nb_element = model.getFEEngine().getMesh().getNbElement(TYPE);

  std::stringstream filename; filename << "material_non_local_" << TYPE;

  dumper.SetMode(iohelper::TEXT);
  dumper.SetParallelContext(StaticCommunicator::getStaticCommunicator()->whoAmI(),
			    StaticCommunicator::getStaticCommunicator()->getNbProc());
  dumper.SetPoints(model.getFEEngine().getMesh().getNodes().storage(),
		   spatial_dimension, nb_nodes, filename.str().c_str());
  dumper.SetConnectivity((int *)model.getFEEngine().getMesh().getConnectivity(TYPE).storage(),
			 paraviewType<TYPE>(), nb_element, iohelper::C_MODE);
  dumper.AddElemDataField(quadrature_points_volumes(TYPE).storage(),
   			  1, "volume");
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
}

/* -------------------------------------------------------------------------- */
// void paraviewDump(iohelper::Dumper & dumper) {
//   dumper.Dump();
// }

/* -------------------------------------------------------------------------- */
#endif
