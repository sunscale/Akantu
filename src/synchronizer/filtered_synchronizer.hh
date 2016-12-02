/**
 * @file   filtered_synchronizer.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Mathilde Radiguet <mathilde.radiguet@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Dec 08 2015
 *
 * @brief  synchronizer that filters elements from another synchronizer
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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

#ifndef __AKANTU_FILTERED_SYNCHRONIZER_HH__
#define __AKANTU_FILTERED_SYNCHRONIZER_HH__

/* -------------------------------------------------------------------------- */
#include "element_synchronizer.hh"

/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
class SynchElementFilter {
public:
  virtual bool operator()(const Element &) = 0;
};

/* -------------------------------------------------------------------------- */
class FilteredSynchronizer : public ElementSynchronizer {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  FilteredSynchronizer(Mesh & mesh,
		       SynchronizerID id = "filtered_synchronizer",
		       MemoryID memory_id = 0);
  virtual ~FilteredSynchronizer() {};
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// get another synchronizer and filter its elements using a functor
  static FilteredSynchronizer * 
  createFilteredSynchronizer(const ElementSynchronizer & d_synchronizer,
			     SynchElementFilter & filter);

protected:
  /// set up the synchronizer
  void setupSynchronizer(const ElementSynchronizer & d_synchronizer,
			 SynchElementFilter & filter);
  /// push source elements into destination elements through the filter
  void updateElementList(Array<Element> * source_elements,
			 Array<Element> * destination_elements,
			 SynchElementFilter & filter);


protected:
  /// Define the receive list tag
  enum CommTags {
    RECEIVE_LIST_TAG  = 0
  };

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  
};

__END_AKANTU__

#endif /* __AKANTU_FILTERED_SYNCHRONIZER_HH__ */
