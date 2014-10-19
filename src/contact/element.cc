/**
 * @file   element.cc
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date creation: Fri Jan 04 2013
 * @date last modification: Tue May 13 2014
 *
 * @brief  contact element classes
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

#include "element.hh"
#include "aka_math.hh"
#include "aka_geometry.hh"

#if defined(AKANTU_BOOST_CHRONO) && !defined(AKANTU_NDEBUG)
#  include <boost/chrono.hpp>
#endif

__BEGIN_AKANTU__


template <>
bool check_penetration<2>(UInt node, const Element* el, SolidMechanicsModel& model) {
  
  typedef Point<2> point_type;
  
  Mesh& mesh = model.getMesh();
  const Array<Real> &x = model.getCurrentPosition();
  const Array<UInt> &conn = mesh.getConnectivity(el->type);
  
  point_type r(&x(node));
  
  // NOTE: switch on a enumerated type is a sign of a bad
  // object-oriented design
  switch (el->type) {
    case _segment_2:
    {
      // get element extreme points
      point_type p(&x(conn(el->element,0)));
      point_type q(&x(conn(el->element,1)));
      
      return left_turn(p, q, r) > 0;
    }
      break;
      
    default:
      cout<<"*** ERROR *** Function signed measure in file "<<__FILE__<<", line "<<__LINE__;
      cout<<", is not implemented for element of type "<<el->type<<endl;
      cout<<"*** ABORTING *** "<<endl;
      exit(1);
      break;
  }
}




template <>
bool check_penetration<3>(UInt node, const Element* el, SolidMechanicsModel& model) {
  
  typedef Point<3> point_type;
  
  Mesh& mesh = model.getMesh();
  const Array<Real> &x = model.getCurrentPosition();
  const Array<UInt> &conn = mesh.getConnectivity(el->type);
  
  point_type r(&x(node));
  
  // NOTE: switch on a enumerated type is a sign of a bad
  // object-oriented design
  switch (el->type) {
      
    case _triangle_3:
    {
      // get element extreme points
      point_type o(&x(conn(el->element,0)));
      point_type p(&x(conn(el->element,1)));
      point_type q(&x(conn(el->element,2)));
      
      // get signed volume
      point_type po = o - p;
      point_type pq = q - p;
      point_type pr = r - p;
      
      // cross product
      point_type c = cross(pq, po);

      Real v = pr*c;
      return v < 0;
    }
      break;
      
    default:
      cout<<"*** ERROR *** Function signed measure in file "<<__FILE__<<", line "<<__LINE__;
      cout<<", is not implemented for element of type "<<el->type<<endl;
      cout<<"*** ABORTING *** "<<endl;
      exit(1);
      break;
  }
}




template <>
Point<2> minimize_distance<2>(UInt node, const Element* el, SolidMechanicsModel& model) {
  
  const UInt d = 2;
  typedef Point<d> point_type;
  
  const Array<Real> &x = model.getCurrentPosition();
  
  point_type r(&x(node));
  
  // NOTE: switch on a enumerated type is a sign of a bad
  // object-oriented design
  switch (el->type) {
    case _segment_2:
    {
#if AKANTU_OPTIMIZATION
      Distance_minimzer<_segment_2> data(&x(node), el, model);
      return data.point();
#else 
      AKANTU_DEBUG_ERROR("To use this function you should activate the optimization at compile time");
      return Point<2>();
#endif
    }
      break;
      
    default:
      cout<<"*** ERROR *** Function signed measure in file "<<__FILE__<<", line "<<__LINE__;
      cout<<", is not implemented for element of type "<<el->type<<endl;
      cout<<"*** ABORTING *** "<<endl;
      exit(1);
      break;
  }
  
  assert(false);
  return point_type(); // avoid compiler warning: control reaches end of non-void function
}



template <>
Point<3> minimize_distance<3>(UInt node, const Element* el, SolidMechanicsModel& model) {

  const UInt d = 3;
  typedef Point<d> point_type;
  
  const Array<Real> &x = model.getCurrentPosition();
  
  point_type r(&x(node));
  
  // NOTE: switch on a enumerated type is a sign of a bad
  // object-oriented design
  switch (el->type) {
      
    case _triangle_3:
    {
#if AKANTU_OPTIMIZATION
      Distance_minimzer<_triangle_3> data(&x(node), el, model);

#if defined(AKANTU_BOOST_CHRONO) && !defined(AKANTU_NDEBUG)     
      typedef boost::chrono::high_resolution_clock clock_type;
      typedef typename clock_type::time_point time_type;
      time_type start = clock_type::now();
#endif
      data.optimize();

#if defined(AKANTU_BOOST_CHRONO) && !defined(AKANTU_NDEBUG)
      boost::chrono::nanoseconds ns = clock_type::now() - start;
      cout<<data.iterations()<<"\t"<<ns.count()<<endl;
#endif
//      cout<<data.point()<<endl;
      return data.point();
#else
      AKANTU_DEBUG_ERROR("To use this function you should activate the optimization at compile time");
      return Point<3>();
#endif

      
      
//      const UInt nb_nodes = 3;
//      
//      // get triangle coordinates from element and point coordinates
//      Mesh& mesh = model.getMesh();
//      std::vector<point_type> pts(nb_nodes);
//      
//      const Array<UInt> &conn = mesh.getConnectivity(el->type);
//      for (UInt i=0; i<nb_nodes; ++i)
//        for (UInt j=0; j<d; ++j)
//          pts.at(i)[j] = x(conn(el->element,i),j);
//      
//      // get closest point
//      time_type start = clock_type::now();
//      
//      Point<3> q = naive_closest_point_to_triangle(r,pts[0],pts[1],pts[2]);
//      
//      boost::chrono::nanoseconds ns = clock_type::now() - start;
//      cout<<ns.count()<<endl;
//      
////      static unsigned int k = 0;
////      if (sqrt((q - data.point()).sq_norm()) > 1e-2) {
////        cout<<++k<<endl;
////        cout<<"difference"<<endl;
////        cout<<data.point()<<endl;
////        cout<<q<<endl;
////        cout<<" iter -> "<<data.iterations()<<endl;
////      }
//      
//      
//      return q;
      
            
//      const UInt nb_nodes = 3;
//
//      // get triangle coordinates from element and point coordinates
//      Mesh& mesh = model.getMesh();
//      std::vector<point_type> pts(nb_nodes);
//      
//      const Array<UInt> &conn = mesh.getConnectivity(el->type);
//      for (UInt i=0; i<nb_nodes; ++i)
//        for (UInt j=0; j<d; ++j)
//          pts.at(i)[j] = x(conn(el->element,i),j);
//      
//      // get closest point
//      time_type start2 = clock_type::now();
//
//      Point<3> q = closest_point_to_triangle(r,pts[0],pts[1],pts[2]);
//      
//      boost::chrono::nanoseconds ns2 = clock_type::now() - start2;
//      cout<<ns2.count()<<endl;
//      cout<<q<<endl;
//      
//      static unsigned int k = 0;
//      Real diff = sqrt((q - data.point()).sq_norm());
//      if (diff > 1e-8) {
//        cout<<"diff -> "<<diff<<endl;
//        cout<<++k<<endl;
//        cout<<data.point()<<endl;
//        cout<<q<<endl;
//        cout<<" iter -> "<<data.iterations()<<endl;
//      }
//
//      
//      return q;
    }
      break;
      
      
    case _triangle_6:
    {
#if AKANTU_OPTIMIZATION
      Distance_minimzer<_triangle_6> data(&x(node), el, model);
#if defined(AKANTU_BOOST_CHRONO) && !defined(AKANTU_NDEBUG)      
      time_type start = clock_type::now();
#endif
      data.optimize();
      
#if defined(AKANTU_BOOST_CHRONO) && !defined(AKANTU_NDEBUG)
      boost::chrono::nanoseconds ns = clock_type::now() - start;
      cout<<data.iterations()<<"\t"<<ns.count()<<endl;
#endif
      return data.point();
#else
      AKANTU_DEBUG_ERROR("To use this function you should activate the optimization at compile time");
      return Point<3>();
#endif
    }
      break;
      
    default:
      cout<<"*** ERROR *** Function signed measure in file "<<__FILE__<<", line "<<__LINE__;
      cout<<", is not implemented for element of type "<<el->type<<endl;
      cout<<"*** ABORTING *** "<<endl;
      exit(1);
      break;
  }
  
  assert(false);
  return point_type(); // avoid compiler warning: control reaches end of non-void function
}



__END_AKANTU__
