/**
 * @file   aka_bounding_box.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date creation: Fri Jan 04 2013
 * @date last modification: Mon Jun 02 2014
 *
 * @brief  class for a bounding box
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

#ifndef __AKANTU_AKA_BOUNDING_BOX_HH__
#define __AKANTU_AKA_BOUNDING_BOX_HH__

#include <iostream>
#include <iomanip>

#include "aka_common.hh"
#include "aka_point.hh"

__BEGIN_AKANTU__

using std::cout;
using std::endl;

template <int> class BoundingBox;

/// considers bounding box with respect to a list of points and adaptes it
template <int d, class point_container>
BoundingBox<d> computeBoundingBox(const point_container& points) {
  
  typedef typename point_container::const_iterator point_iterator;
  point_iterator it = points.begin();
  assert(it != points.end());
  
  BoundingBox<d> bbox(*it);
  for (++it; it != points.end(); ++it)
    bbox += *it;
  return bbox;  
}


template <int d, class nodes_container>
BoundingBox<d> createPointList(const nodes_container& nodes, const Array<Real>& coord);


template <int d>
class BoundingBox : public Bounding_volume<d> {

public:
  
  typedef Bounding_volume<d> base_type;
  typedef typename base_type::point_type point_type;
  typedef typename point_type::value_type value_type;

  static int dim()
  { return d; }

private:
  
  /// minimum point
  point_type min_;
  
  /// maximum point
  point_type max_;

public:
    
  /// default constructor, creates an inconsistent bounding box
  BoundingBox() : base_type(), min_(inf), max_(-inf) {}
  
  /// point constructor, sets the bounding box to the point
  BoundingBox(const point_type& p1) : base_type(), min_(p1), max_(p1) {}
  
  /// two-point constructor, calculates minimum and maximum points
  BoundingBox(const point_type& p1, const point_type& p2, bool compute = true) 
  : base_type(), min_(p1), max_(p2)
  {
    if (compute)
      for (Int i=0; i<d; ++i) {
        min_[i] = std::min(p1[i], p2[i]);
        max_[i] = std::max(p1[i], p2[i]);
      }
  }

  /// multiple-point constructor, creates a bounding box encompass multiple points
  template <class iterator>
  BoundingBox(iterator first, iterator last) : base_type(), min_(*first), max_(*first) {
    ++first;
    for (; first != last; ++first)
      this->operator+=(*first);
  }
  
  /// returns the measure of the bounding box: the measure is the volume of the box
  value_type measure() const {
    value_type v = 1;
    for (int i=0; i<d; ++i)
      v *= (max_[i] - min_[i]);
    return v;
  }
  
  /// 
  bool operator<(const BoundingBox& bbox) const {
    return min_ < bbox.min_ || (!(bbox.min_ < min_) && max_ < bbox.max_);
  }
  
  /**
     returns whether the two bounding boxes have min points with the same coordinates 
     and max points with the same coordinates
  **/
  bool operator==(const BoundingBox& bbox) const
  { return min_ == bbox.min_ && max_ == bbox.max_; }
  
  /**
     returns whether the two bounding boxes have min points with different coordinates 
     or max points with different coordinates
  **/
  bool operator!=(const BoundingBox& bbox) const
  { return !(*this == bbox); }
  
  /// extend the bounding box, if necessary, in order to encompass this additional point
  BoundingBox& operator+=(const point_type& point) {
    for (Int i=0; i<d; ++i) {
      min_[i] = std::min(min_[i], point[i]);
      max_[i] = std::max(max_[i], point[i]);
    }
    return *this;
  }
  
  /// extend the bounding box, if necessary, in order to encompass the entire given boundary box 
  BoundingBox& operator+=(const BoundingBox& bbox) {
    this->operator+=(bbox.min_);
    this->operator+=(bbox.max_);
    return *this;
  }
  
  /// is the point p geometrically inside of the bounding box?
  bool operator&(const point_type& p) const {
    
    Real e = 2*std::numeric_limits<Real>::epsilon();
    
    for (Int i=0; i<d; ++i)
      if (max_[i] < p[i] - e || p[i] < min_[i] - e)
        return false;
    return true;
  }
  
  /// is the bounding box bb geometrically entirely inside of the bounding box?
  bool operator&(const BoundingBox& bb) const {
    
    Real e = 2*std::numeric_limits<Real>::epsilon();
    
    for (Int i=0; i<d; ++i) {
      if (max_[i] < bb.min_[i] - e || bb.max_[i] < min_[i] - e)
        return false;
    }
    return true;
  }
  
  /// create the intersection bounding box
  BoundingBox operator&&(const BoundingBox& bb) const {
    
    BoundingBox intersection;
    for (Int i=0; i<d; ++i) {
      intersection.min_[i] = std::max(min_[i], bb.min_[i]);
      intersection.max_[i] = std::min(max_[i], bb.max_[i]);
    }
    return intersection;
  }
  
  /// get the point of the minimum corner
  const point_type& min()
  { return min_; }
  
  /// get the point of the maximum corner
  const point_type& max()
  { return max_; }
  
  /// get the point of the minimum corner
  point_type min() const
  { return min_; }
  
  /// get the point of the maximum corner
  point_type max() const
  { return max_; }
  
  /// get minimum coordinate of the bounding box in direction i
  Real min(size_t i) const
  { return min_[i]; }
  
  /// get maximum coordinate of the bounding box in direction i
  Real max(size_t i) const
  { return max_[i]; }
  
  virtual std::ostream& print(std::ostream& os) const {
    
    os<<*this;
    return os;
  }

public:
  
  /// directional increase of bounding box (if needed)
  void expand(Real coord, UInt dir) {
    AKANTU_DEBUG_ASSERT(dir < d, "");
  }
  
};

/// cumulate two bounding boxes and create one encompassing both
template <int d>
BoundingBox<d> operator+(const BoundingBox<d>& b1, const BoundingBox<d>& b2) {
  BoundingBox<d> r(b1);
  return r += b2;
}

template <int d>
std::ostream& operator<<(std::ostream&, const BoundingBox<d>&);



__END_AKANTU__

#endif /* __AKANTU_AKA_BOUNDING_BOX_HH__ */
