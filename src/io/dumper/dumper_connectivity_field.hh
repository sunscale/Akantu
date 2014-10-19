/**
 * @file   dumper_connectivity_field.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Tue Sep 02 2014
 * @date last modification: Thu Sep 04 2014
 *
 * @brief  Connectivity field dumper
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

#if defined(AKANTU_COHESIVE_ELEMENT)
/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__
__BEGIN_AKANTU_DUMPER__
/* -------------------------------------------------------------------------- */
template <class types>
class cohesive_connectivity_field_iterator : 
  public element_iterator<types, cohesive_connectivity_field_iterator> {

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */  
  
public:
  typedef element_iterator<types, dumper::cohesive_connectivity_field_iterator> parent;
  typedef typename types::return_type return_type;
  typedef typename types::field_type field_type;
  typedef typename types::array_iterator array_iterator;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

public:
  cohesive_connectivity_field_iterator(const field_type & field,
				       const typename field_type::type_iterator & t_it,
				       const typename field_type::type_iterator & t_it_end,
				       const array_iterator & array_it,
				       const array_iterator & array_it_end,
				       const GhostType ghost_type = _not_ghost) :
    parent(field, t_it, t_it_end, array_it, array_it_end,ghost_type) {

    write_order[_cohesive_3d_12].push_back(0);
    write_order[_cohesive_3d_12].push_back(1);
    write_order[_cohesive_3d_12].push_back(2);
    write_order[_cohesive_3d_12].push_back(6);
    write_order[_cohesive_3d_12].push_back(7);
    write_order[_cohesive_3d_12].push_back(8);
    write_order[_cohesive_3d_12].push_back(3);
    write_order[_cohesive_3d_12].push_back(4);
    write_order[_cohesive_3d_12].push_back(5);
    write_order[_cohesive_3d_12].push_back(9);
    write_order[_cohesive_3d_12].push_back(10);
    write_order[_cohesive_3d_12].push_back(11);

    write_order[_cohesive_3d_6].push_back(0);
    write_order[_cohesive_3d_6].push_back(1);
    write_order[_cohesive_3d_6].push_back(2);
    write_order[_cohesive_3d_6].push_back(3);
    write_order[_cohesive_3d_6].push_back(4);
    write_order[_cohesive_3d_6].push_back(5);

    write_order[_cohesive_2d_6].push_back(0);
    write_order[_cohesive_2d_6].push_back(2);
    write_order[_cohesive_2d_6].push_back(1);
    write_order[_cohesive_2d_6].push_back(4);
    write_order[_cohesive_2d_6].push_back(5);
    write_order[_cohesive_2d_6].push_back(3);

    write_order[_cohesive_2d_4].push_back(0);
    write_order[_cohesive_2d_4].push_back(1);
    write_order[_cohesive_2d_4].push_back(3);
    write_order[_cohesive_2d_4].push_back(2);

  }

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
  

  return_type operator*() {
    const ElementType & type = *this->tit;
    const Vector<UInt> & conn = *this->array_it;
    Vector<UInt> new_conn(conn.size());

    for (UInt n = 0; n < conn.size(); ++n)
      new_conn(n) = conn(write_order[type][n]);

    return new_conn;
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
  

protected:

  std::map<ElementType, std::vector<UInt> > write_order;

};


/* -------------------------------------------------------------------------- */
class CohesiveConnectivityField : 
  public GenericElementalField<SingleType<UInt,Vector,false>,
			       cohesive_connectivity_field_iterator> {

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */  
  
public:

  typedef SingleType<UInt,Vector,false> types;
  typedef cohesive_connectivity_field_iterator<types> iterator;
  typedef GenericElementalField<types,cohesive_connectivity_field_iterator> parent;


  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

public:

  CohesiveConnectivityField(const field_type & field,
			    UInt spatial_dimension = _all_dimensions,
			    GhostType ghost_type = _not_ghost) :
    parent(field, spatial_dimension, ghost_type, _ek_cohesive) { }


};

/* -------------------------------------------------------------------------- */


__END_AKANTU_DUMPER__
__END_AKANTU__

#endif
