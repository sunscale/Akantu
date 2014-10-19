/**
 * @file   contact_impl.hh
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date creation: Mon Sep 15 2014
 * @date last modification: Mon Sep 15 2014
 *
 * @brief  Contact interface class
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

//
//  contact.hh
//  akantu
//
//  Created by Alejandro Marcos Aragón on 6/24/14.
//
//

#ifndef __AKANTU_CONTACT_HH__
#define __AKANTU_CONTACT_HH__


#include "contact_common.hh"
#include "parser.hh"
#include "search.hh"
#include "resolution.hh"

__BEGIN_AKANTU__

template <int Dim, template <int> class Search_policy, class RP>
class Contact : public Search_policy<Dim>, public ContactResolution <Dim, RP::analysis, RP::method> {
  
  typedef SolidMechanicsModel model_type;
  typedef Search_policy<Dim> search_type;
  typedef ContactResolution <Dim, RP::analysis, RP::method> resolution_type;

public:
  Contact(int argc, char *argv[], model_type & m) : search_type(m), resolution_type(m)
  {
    // read parameters from file
    std::pair <Parser::const_section_iterator, Parser::const_section_iterator>
      sub_sect = getStaticParser().getSubSections(_st_contact);

    if (sub_sect.first != sub_sect.second)
      this->parseSection(*sub_sect.first);

    // read parameters the command line
    contact_argparser.parse(argc, argv, cppargparse::_remove_parsed);
    
    // finish initialization of resolution class
    this->initialize();
  }


  //! Provide standard output of contact object
  friend std::ostream& operator << (std::ostream & os, const Contact &cd) {
    const resolution_type& r(cd);
    const search_type& s(cd);

    os << "\nContact object info:" << endl;
    os << "  Search type: " << s << endl;
    os << "  Resolution type: " << r << endl;
    return os;
  }
};


template <ContactImplementationMethod i, class contact_type>
void solveContactStep(contact_type& c)
{ c.solveContactStep<i>(&c); }



__END_AKANTU__


#endif /* __AKANTU_CONTACT_HH__ */
