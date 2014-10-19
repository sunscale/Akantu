/**
  * @file   aka_debug_tools.hh
  *
  * @author Nicolas Richart <nicolas.richart@epfl.ch>
  *
  * @date   Mon Apr  8 16:33:02 2013
  *
  * @brief  Different tool to help to debug (compiled only in AKANTU_DEBUG mode)
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

/* -------------------------------------------------------------------------- */
#if defined(AKANTU_CORE_CXX11)
#  include <functional>
#endif
#include <iostream>
#include <fstream>

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "static_communicator.hh"
/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_AKA_DEBUG_TOOLS_HH__
#define __AKANTU_AKA_DEBUG_TOOLS_HH__

#if defined(AKANTU_DEBUG_TOOLS)

#define AKANTU_DEBUG_DBT(info)  AKANTU_DEBUG_("DBT", dblSecondary, info)

namespace akantu {
  namespace debug {
    enum DebugModule {
      _dm_all,
      _dm_test,
      _dm_synch,
      _dm_material,
      _dm_material_non_local,
      _dm_material_damage,
      _dm_material_cohesive,
      _dm_model_cohesive,
      _dm_debug_tools,
      _dm_integrator,
      _dm_end
    };

    struct DebugModule_def {
      typedef DebugModule type;
      static const type _begin_ = _dm_all;
      static const type _end_   = _dm_end;
    };

    typedef safe_enum<DebugModule_def> debug_module_t;


    class DebugModulesHandler {
    public:
      void addModule(DebugModule mod) {
	if(mod == _dm_all)
	  for(debug_module_t::iterator m(debug_module_t::begin()); m != debug_module_t::end(); ++m)
	    modules.insert(*m);
	else
	  modules.insert(mod);
      }

      void removeModule(DebugModule mod) {
	if(mod == _dm_all)
	  modules.clear();
	else
	  modules.erase(mod);
      }

      bool isActive(DebugModule mod) {
	return modules.find(mod) != modules.end();
      }
    private:
      std::set<DebugModule> modules;
    };


    /* ---------------------------------------------------------------------- */
    /* Debug Element Manager                                                  */
    /* ---------------------------------------------------------------------- */
    class DebugElementManager : public DebugModulesHandler, public MeshEventHandler {
    public:
      DebugElementManager() : barycenters("debug_element_manager", "debug_element_barycenters", 0), out(&std::cerr) {
      }

      void setOutStream(std::ostream & o) {
	out = &o;
      }

      void setMesh(const Mesh & imesh) {
	this->mesh = &imesh;
	UInt spatial_dimension = mesh->getSpatialDimension();
	mesh->initElementTypeMapArray(barycenters,
				     spatial_dimension,
				     _all_dimensions,
				     false,
				     _ek_not_defined);

	for (ghost_type_t::iterator gt = ghost_type_t::begin();  gt != ghost_type_t::end(); ++gt) {
	  Mesh::type_iterator it = mesh->firstType(_all_dimensions,
						   *gt,
						   _ek_not_defined);

	  Mesh::type_iterator last_type = mesh->lastType(_all_dimensions,
							 *gt,
							 _ek_not_defined);
	  for(; it != last_type; ++it) {
	    UInt nb_element = mesh->getNbElement(*it, *gt);
	    Array<Real> & barycenter = barycenters(*it, *gt);
	    barycenter.resize(nb_element);

	    Array<Real>::vector_iterator bary_it = barycenter.begin(spatial_dimension);
	    for (UInt elem = 0; elem < nb_element; ++elem) {
	      mesh->getBarycenter(elem, *it, bary_it->storage(), *gt);
	      ++bary_it;
	    }
	  }
	}
	const_cast<Mesh *>(mesh)->registerEventHandler(*this);
      }

    public:
      void clear() { element_to_debug.clear(); }

      // Add element
      bool addElement(const Element & element) {
	element_to_debug.insert(element);
	if(isActive(_dm_debug_tools))
	  AKANTU_DEBUG_DBT(element << ": registered in debugger");
	return true;
      }

      bool addElement(const Vector<Real> & pos, const ElementType & type,
		      const GhostType & ghost_type, Real tolerance) {
	Array<Real> bary_arr = barycenters(type, ghost_type);
	Array<Real>::const_iterator<Vector<Real> > bary_begin = bary_arr.begin(bary_arr.getNbComponent());
	Array<Real>::const_iterator<Vector<Real> > bary_it = bary_begin;
	Array<Real>::const_iterator<Vector<Real> > bary_end = bary_arr.end(bary_arr.getNbComponent());
	bool found = false;
	for(; bary_it != bary_end && !found; ++bary_it) {
	  found = pos.equal(*bary_it, tolerance);
	}

	if(!found) return false;

	UInt el_id = bary_it - bary_begin - 1;
	Element el(type, el_id, ghost_type, Mesh::getKind(type));
	return addElement(el);
      }

      bool addElement(const Vector<Real> & pos, Real tolerance) {
	bool finished = false;
	for (ghost_type_t::iterator gt = ghost_type_t::begin();  gt != ghost_type_t::end() && ! finished; ++gt) {
	  Mesh::type_iterator first = mesh->firstType(_all_dimensions, *gt, _ek_not_defined);
	  Mesh::type_iterator last  = mesh->lastType(_all_dimensions, *gt, _ek_not_defined);
	  for (; first != last && !finished; ++first) {
	    ElementType type = *first;
	    finished = addElement(pos, type, *gt, tolerance);
	  }
	}

	barycenters_of_elements_to_debug[pos] = tolerance;
	return finished;
      }

      inline void onElementsAdded(const Array<Element> & element_list,
				  __attribute__((unused)) const NewElementsEvent & event) {
	Array<Element>::const_iterator<Element> it  = element_list.begin();
	Array<Element>::const_iterator<Element> end = element_list.end();

	for (UInt el = 0; it != end; ++it, ++el) {
	  const Element & elem = *it;
	  Vector<Real> bary(mesh->getSpatialDimension(elem.type));
	  mesh->getBarycenter(elem, bary);
	  barycenters(elem.type, elem.ghost_type).push_back(bary);
	  std::map<Vector<Real>, Real>::iterator bit = barycenters_of_elements_to_debug.begin();
	  std::map<Vector<Real>, Real>::iterator bend = barycenters_of_elements_to_debug.end();
	  for (;bit != bend; ++bit) {
	    if(bary.equal(bit->first, bit->second)) addElement(elem);
	  }
	}
      }

      inline void onElementsRemoved(const Array<Element> & element_list,
				    const ElementTypeMapArray<UInt> & new_numbering,
				    __attribute__((unused)) const RemovedElementsEvent & event) {
	barycenters.onElementsRemoved(new_numbering);
	Array<Element>::const_iterator<Element> it  = element_list.begin();
	Array<Element>::const_iterator<Element> end = element_list.end();

	for (; it != end; ++it) {
	  std::set<Element>::iterator el = element_to_debug.find(*it);
	  if(el != element_to_debug.end()) {
	    if(isActive(_dm_debug_tools))
	      AKANTU_DEBUG_DBT(*it << ": unregistered from debuger");
	    element_to_debug.erase(el);
	  }
	}
      }

    public:
      // Actions
      template<typename T>
      void printData(DebugModule mod,
		     const std::string & ctxt,
		     const ElementTypeMapArray<T> & data) {
	if(!isActive(mod)) return;
	begin(ctxt);
	for(std::set<Element>::const_iterator it(element_to_debug.begin());
	    it != element_to_debug.end(); ++it) {
	  std::stringstream sout;
	  try {
	    const Array<T> & array = data(it->type, it->ghost_type);
	    UInt nb_data_per_element
	      = array.getSize() / mesh->getNbElement(it->type, it->ghost_type);

	    sout << " -> " << *it << " " << data.getID() << "(" << nb_data_per_element << ")" << ":";
	    typename Array<T>::template const_iterator< Vector<T> > data_it = array.begin(array.getNbComponent());
	    data_it += it->element * nb_data_per_element;
	    for (UInt d = 0; d < nb_data_per_element; ++d, ++data_it) {
	      sout << std::setprecision(15) << " " << *data_it;
	    }
	  } catch (...) {
	    sout << " -> " << *it << " has no data defined";
	  }

	  AKANTU_DEBUG_DBT(sout.str());
	}
	end(ctxt);
      }

      /* ---------------------------------------------------------------------- */
      template<typename T>
      void printData(DebugModule mod,
		     const std::string & ctxt,
		     const ElementTypeMapArray<T> & data,
		     const ElementTypeMapArray<UInt> & filter) {
	if(!isActive(mod)) return;

	begin(ctxt);
	for(std::set<Element>::const_iterator it(element_to_debug.begin());
	    it != element_to_debug.end(); ++it) {
	  std::stringstream sout;
	  sout << " -> " << *it << " " << data.getID() << ":";
	  UInt pos = filter(it->type, it->ghost_type).find(it->element);
	  if(pos == UInt(-1)) (*out) << " not in filter!";
	  else {
	    UInt nb_element = filter(it->type, it->ghost_type).getSize();
	    const Array<T> & array = data(it->type, it->ghost_type);
	    UInt nb_data_per_element = array.getSize() / nb_element;
	    typename Array<T>::template const_iterator< Vector<T> > data_it = array.begin(array.getNbComponent());
	    data_it += pos * nb_data_per_element;
	    for (UInt d = 0; d < nb_data_per_element; ++d, ++data_it) {
	      sout << " " << *data_it;
	    }
	  }
	  AKANTU_DEBUG_DBT(sout.str());
	}

	end(ctxt);
      }

#if defined(AKANTU_CORE_CXX11)
      void print(DebugModule mod,
		 std::function< std::string(const Element &) > funct) {
	if(!isActive(mod)) return;

	for(std::set<Element>::const_iterator it(element_to_debug.begin());
	    it != element_to_debug.end(); ++it) {
	  std::string func_str = funct(*it);
	  if(func_str != "")
	    AKANTU_DEBUG_DBT(" -> " << *it << ": " << func_str);
	}
      }
#endif

      // Actions
      void printElements(DebugModule mod, const std::string & ctxt) {
	if(!isActive(mod)) return;

	begin(ctxt);
	for(std::set<Element>::const_iterator it(element_to_debug.begin());
	    it != element_to_debug.end(); ++it) {
	  AKANTU_DEBUG_DBT(" -> " << *it);
	}
	end(ctxt);
      }

      void begin(const std::string & ctxt) {
	dbg_out = &debugger.getOutStream();
	debugger.setOutStream(*out);
	if(! element_to_debug.empty())
	  AKANTU_DEBUG_DBT("Begin ("<< ctxt <<")");
      }

      void end(const std::string & ctxt) {
	if(! element_to_debug.empty())
	  AKANTU_DEBUG_DBT("End ("<< ctxt <<")");
	debugger.setOutStream(*dbg_out);
      }

    protected:
      const Mesh * mesh;
      ElementTypeMapArray<Real> barycenters;
      std::set<Element> element_to_debug;
      std::map<Vector<Real>, Real> barycenters_of_elements_to_debug;
      std::ostream * out;
      std::ostream * dbg_out;
    };

    extern DebugElementManager element_manager;
  }

  /* ------------------------------------------------------------------------ */
  /* CSV Writer                                                               */
  /* ------------------------------------------------------------------------ */
  class csv_ofstream {
  public:
    csv_ofstream(char sep = ',') : sep(sep) {
    }

    ~csv_ofstream() {
      if(file.is_open()) file.close();
    }

    void open(std::string filename,
	      std::ios_base::openmode mode = std::ios_base::out) {
      file.open(filename.c_str(), mode);
    }

    template<class T>
    csv_ofstream & operator<<(T t) {
      file << t << sep;
      return *this;
    }

    csv_ofstream & operator<<(std::ostream& ( *pf )(std::ostream&)) {
      long pos = file.tellp();
      file.seekp(pos-1);
      file << pf;
      return *this;
    }

  private:
    std::ofstream file;
    char sep;
  };

}


#endif
#endif /* __AKANTU_AKA_DEBUG_TOOLS_HH__ */
