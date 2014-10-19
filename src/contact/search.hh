/**
 * @file   search.hh
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date creation: Fri Jan 04 2013
 * @date last modification: Mon Sep 15 2014
 *
 * @brief  global contact search
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

#ifndef __AKANTU_SEARCH_HH__
#define __AKANTU_SEARCH_HH__

#include <iostream>

#include "aka_common.hh"
#include "solid_mechanics_model.hh"
#include "aka_point.hh"
#include "aka_geometry.hh"


__BEGIN_AKANTU__



struct SearchBase {
	template <class model_type>
	SearchBase(model_type & model) {
	}

	virtual int search(const Real *) {
		static bool msg = false;
		if (!msg) {
			std::cout << " - Warning: calling default base searcher, type any key to continue." << std::endl;
			std::cin.ignore();
			msg = true;
		}
		return -1;
	}

	virtual ~SearchBase()
	{
	}

	//! Provide standard output of contact object
	friend std::ostream& operator << (std::ostream & os, const SearchBase &cr) {
		os << "Search base class (no search performed)" << std::endl;
		return os;
		}
};


template <int dim, class model_type>
class SearchTraits;

template <class model_type>
struct SearchTraits <3, model_type> {
	typedef Point <3> point_type;
	typedef ModelElement <model_type> element_type;

	static bool check_projection(const point_type& p, model_type& model, UInt id) {
		element_type m(model, _triangle_3, id);

		return point_has_projection_to_triangle(p, m.template point <3>(0), m.template point <3>(1), m.template point <3>(2));
	}
};


template <class model_type>
struct SearchTraits <2, model_type> {
	typedef Point <2> point_type;
	typedef ModelElement <model_type> element_type;

	static bool check_projection(const point_type& p, model_type& model, UInt id) {
		element_type m(model, _segment_2, id);

		return has_projection(p, m.template point <2>(0), m.template point <2>(1));
	}
};


template <int dim>
struct MasterAssignator : public SearchBase {
	typedef Point <dim> point_type;
	typedef SolidMechanicsModel model_type;

	std::list <std::string> surfaces_;
	model_type &model_;

	MasterAssignator(model_type & m) : SearchBase(m), surfaces_(), model_(m)
	{
	}


	virtual int search(const Real *ptr) {
		point_type p(ptr);

		for (auto surface: surfaces_) {
			ElementGroup &rs = model_.getMesh().getElementGroup(surface);

			for (ElementGroup::type_iterator tit = rs.firstType(); tit != rs.lastType(); ++tit)
				for (ElementGroup::const_element_iterator it = rs.element_begin(*tit);
				     it != rs.element_end(*tit); ++it)
					if (SearchTraits <dim, model_type>::check_projection(p, model_, *it))
						return *it;
		}
		return -1;
	}

	void searchSurface(const std::string& s) {
		surfaces_.push_back(s);
	}

	//! Provide standard output of contact object
	friend std::ostream& operator << (std::ostream & os, const MasterAssignator &s) {
		os << "MasterAssignator search class. Search surfaces: ";

		for (auto surface : s.surfaces_)
			os << " " << surface;
		os << std::endl;
		return os;
		}
};



__END_AKANTU__

#endif /* __AKANTU_SEARCH_HH__ */
