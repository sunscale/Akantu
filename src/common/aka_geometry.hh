/**
 * @file   aka_geometry.hh
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date creation: Fri Jan 04 2013
 * @date last modification: Tue Sep 02 2014
 *
 * @brief  geometric operations
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

/* -------------------------------------------------------------------------- */


#ifndef __AKANTU_GEOMETRY_HH__
#define __AKANTU_GEOMETRY_HH__

#include <iostream>
#include <tuple>
#include "aka_point.hh"
#include "aka_plane.hh"
#include "aka_math.hh"

__BEGIN_AKANTU__

// predicates


// combined tolerance test, from Christer Ericson
template <typename T>
typename std::enable_if<std::is_floating_point<T>::value, bool>::type
equal(T x, T y, T tol = 2*std::numeric_limits<T>::epsilon()) {
  
  T absTol = tol;
  T relTol = absTol;

  // here both tolerances are equal, but the code is written
  // like this so that tolerance values can be assigned indepently
  // in the future
  return std::abs(x - y) <= std::max(absTol, relTol * std::max(std::abs(x), std::abs(y)));
}



// combined tolerance test, from Christer Ericson
template <typename T>
typename std::enable_if<std::is_integral<T>::value, bool>::type
equal(T x, T y)
{ return x == y; }


Real left_turn(const Point<2>& p, const Point<2>& q, const Point<2>& r);



// closest point computations


//! Computes the closest point laying on a segment to a point
/*! Given segment \c ab and point \c c, computes closest point \c d on ab.
 * Also returns \c t for the position of the point: a + t*(b - a)
 */
template <int d, typename T>
Point<d,T> closest_point_to_segment(const Point<d,T>& c,
                                    const Point<d,T>& a,
                                    const Point<d,T>& b) {
  
  Point<d,T> ab = b - a;
  
  
  // project c onto ab, computing parameterized position d(t) = a + t*(b – a)
  
  T t = (c - a)*ab / sqrt(ab*ab);
  
  // if outside segment, clamp t (and therefore d) to the closest endpoint
  if (t < 0.)
    t = 0.;
  else if (t > 1.)
    t = 1.;
  
  // compute projected position from the clamped t
  return a + t * ab;
}

//! Predicate that checks if a point has a projection on a line segment
/*! Given segment \c ab and point \c c, checks if the point has a projection in the segment.
 */
template <int d, typename T>
bool has_projection(const Point<d,T>& c,
                    const Point<d,T>& a,
                    const Point<d,T>& b) {
  
  Point<d,T> ab = b - a;
  // project c onto ab, computing parameterized position d(t) = a + t*(b – a)
  
  T t = (c - a)*ab / (ab*ab);
  return t > 0. && t < 1.;
}

//! Tests if a point has a projection to a triangle
/*! This function uses the concept of Voronoi regions to determine
 * if a point has a projection within a triangle defined by points
 * \c a, \c b, and \c c.
 */
template <typename T>
bool point_has_projection_to_triangle(const Point<3,T>& p,
                                      const Point<3,T>& a,
                                      const Point<3,T>& b,
                                      const Point<3,T>& c) {
  
  typedef Point<3,T> point_type;
  
  // obtain plane of the triangle
  Plane pi(a,b,c);
  
  // get point in the plane closest to p
  point_type q = closest_point_to_plane(p,pi);
  
  // return if point is within the triangle
  if (is_point_in_triangle(q, a, b, c))
    return true;
  return false;
}

//! Tests if point P lies inside a triangle
/*! The triangle is defined by points \c a, \c b and \c c.
 */

template <typename T>
bool is_point_in_triangle(const Point<3,T>& p,
                          const Point<3,T>& a,
                          const Point<3,T>& b,
                          const Point<3,T>& c) {

  typedef Point<3,T> point_type;
  
  point_type v0 = b-a, v1 = c-a, v2 = p-a;
  
  Real d00 = v0*v0;
  Real d01 = v0*v1;
  Real d11 = v1*v1;
  Real d20 = v2*v0;
  Real d21 = v2*v1;
  Real denom = d00*d11 - d01*d01;
  
  // compute parametric coordinates
  Real v = (d11 * d20 - d01 * d21) / denom;
  Real w = (d00 * d21 - d01 * d20) / denom;  
  return v >= 0. && w >= 0. && v + w <= 1.;
}



//! Compute the closest point to a triangle
/*! This function uses the concept of Voronoi regions to determine
 * the closest point \c p to a triangle defined by points \c a, \c b
 * \c c.
 */
template <typename T>
Point<3,T> closest_point_to_triangle(const Point<3,T>& p,
                                     const Point<3,T>& a,
                                     const Point<3,T>& b,
                                     const Point<3,T>& c) {
  
  typedef Point<3,T> point_type;
  
  // check if P in vertex region outside A
  point_type ab = b - a;
  point_type ac = c - a;
  point_type ap = p - a;
  
  // compute scalar products
  T d1 = ab * ap;
  T d2 = ac * ap;
  
  if (d1 <= 0. && d2 <= 0.)
    return a; // barycentric coordinates (1,0,0)
  
  // check if P in vertex region outside B
  point_type bp = p - b;
  
  T d3 = ab * bp;
  T d4 = ac * bp;
  
  if (d3 >= 0.0f && d4 <= d3)
    return b; // barycentric coordinates (0,1,0)
  
  // check if P in edge region of AB, if so return projection of P onto AB
  T vc = d1*d4 - d3*d2;
  if (vc <= 0. && d1 >= 0. && d3 <= 0.) {
    T v = d1 / (d1 - d3);
    return a + v * ab; // barycentric coordinates (1-v,v,0)
  }
  
  // check if P in vertex region outside C
  point_type cp = p - c;
  T d5 = ab * cp;
  T d6 = ac * cp;
  if (d6 >= 0.0f && d5 <= d6)
    return c; // barycentric coordinates (0,0,1)
  
  // check if P in edge region of AC, if so return projection of P onto AC
  T vb = d5*d2 - d1*d6;
  if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
    T w = d2 / (d2 - d6);
    return a + w * ac; // barycentric coordinates (1-w,0,w)
  }
  
  // Check if P in edge region of BC, if so return projection of P onto BC
  T va = d3*d6 - d5*d4;
  if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) {
    T w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    return b + w * (c - b); // barycentric coordinates (0,1-w,w)
  }
  
  // P inside face region. Compute Q through its barycentric coordinates (u,v,w)
  T denom = 1.0f / (va + vb + vc);
  T v = vb * denom;
  T w = vc * denom;
  
  return a + ab*v + ac*w; // = u*a + v*b + w*c, u = va*denom = 1.0f - v - w
}


template <typename T>
Point<3,T> closest_point_to_plane(const Point<3,T>& q, const Plane& p) {
  
  typedef Point<3,T> point_type;
  
  const point_type& n = p.normal();
  
  T t = (n*q - p.distance()) / (n*n);
  return q - t * n;
}

//! Compute the closest point to a triangle
/*! Obtains the plane of the triangle and checks if the point lies inside the
 * triangle. If not, it computes the closest point to each of the triangle 
 * edges.
 */
template <typename T>
Point<3,T> naive_closest_point_to_triangle(const Point<3,T>& p,
                                           const Point<3,T>& a,
                                           const Point<3,T>& b,
                                           const Point<3,T>& c) {
  
  typedef Point<3,T> point_type;
  
  // obtain plane of the triangle
  Plane pi(a,b,c);
  
  // get point in the plane closest to p
  point_type q = closest_point_to_plane(p,pi);
  
  // return if point is within the triangle
  if (is_point_in_triangle(q, a, b, c))
    return q;
  
  // else get the closest point taking into account all edges
  
  // first edge
  q = closest_point_to_segment(p, a, b);
  T d = (q-p).sq_norm();
  
  // second edge
  point_type r = closest_point_to_segment(p, b, c);
  
  T d2 = (r-p).sq_norm();
  if (d2 < d) {
    q = r;
    d = d2;
  }
  
  // third edge
  r = closest_point_to_segment(p,c,a);
  d2 = (r-p).sq_norm();
  if (d2 < d)
    q = r;
  
  // return closest point
  return q;
}


// intersect point p with velocity v with plane
// the function returns collision time and point of contact
// this function does not consider acceleration
template <typename T>
std::tuple<Real, Point<3,T> >
moving_point_against_plane(const Point<3,T>& p, const Point<3,T>& v, Plane& pi) {
  
  typedef Point<3,T> point_type;
  
  // compute distance of point to plane
  Real dist = pi.normal()*p - pi.distance();
  
  // if point already in the plane
  if (std::abs(dist) <= 1e-10)
    return std::make_tuple(0., p);
  
  else {
    
    Real denom = pi.normal()*v;
    // no intersection as poin moving parallel to or away from plane
    if (denom * dist >= 0.)
      return std::make_tuple(inf, point_type());
    
    // point moving towards the plane
    else {
      // point is moving towards the plane
      Real t = -dist/denom;
      return std::make_tuple(t, p + t*v);
    }
  }
}



template <int dim, typename T>
std::tuple<Real, Point<dim,T> >
moving_point_against_point(const Point<dim,T>& s1, const Point<dim,T>& s2, /* point centers */
                           const Point<dim,T>& v1, const Point<dim,T>& v2) /* point velocities */ {
  
  typedef Point<dim,T> point_type;
  typedef typename Point<dim,T>::value_type value_type;
  
  // vector between points
  point_type s = s2 - s1;
  // relative motion of s1 with respect to stationary s0
  point_type v = v2 - v1;
  value_type c = s*s;
  
  // if points within tolerance
  if (equal(s.sq_norm(), value_type()))
    return std::make_tuple(value_type(), s1);
    
  value_type epsilon = 2*std::numeric_limits<T>::epsilon();;
  
  value_type a = v*v;
  // if points not moving relative to each other
  if (a < epsilon)
    return std::make_tuple(inf, point_type());
  
  value_type b = v*s;
  // if points not moving towards each other
  if (b >= 0.)
    return std::make_tuple(inf, point_type());
  
  value_type d = b*b - a*c;
  // if no real-valued root (d < 0), points do not intersect
  if (d >= 0.) {
    value_type ts = (-b - sqrt(d))/a;
    point_type q = s1+v1*ts;
    return std::make_tuple(ts, q);
  }
  return std::make_tuple(inf, point_type());
}


__END_AKANTU__


#endif /* __AKANTU_GEOMETRY_HH__ */
