/**
 * @file   aka_ci_string.hh
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date creation: Fri Jan 04 2013
 * @date last modification: Fri Jan 04 2013
 *
 * @brief  Case insensitive string
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

#ifndef __AKANTU_CI_STRING_HH__
#define __AKANTU_CI_STRING_HH__

__BEGIN_AKANTU__


/// Traits class for a case insensitive string class
struct ci_char_traits : public std::char_traits<char> {
    
    static bool eq( char c1, char c2 )
    { return toupper(c1) == toupper(c2); }
    
    static bool ne( char c1, char c2 )
    { return toupper(c1) != toupper(c2); }
    
    static bool lt( char c1, char c2 )
    { return toupper(c1) <  toupper(c2); }
    
    static int compare(const char* s1, const char* s2, size_t n )
    { return memicmp( s1, s2, n ); }
    
private:
    static int memicmp(const void *s1, const void *s2, size_t n) {
        
        if (n != 0) {
            const unsigned char *p1 = (const unsigned char *)s1, *p2 = (const unsigned char *)s2;
            do {
                if (toupper(*p1) != toupper(*p2))
                    return (*p1 - *p2);
                p1++;
                p2++;
            } while (--n != 0);
        }
        return 0;
    }
};

/// case insensitive string type definition
typedef std::basic_string<char, ci_char_traits> ci_string;


/// provide standard output for case insensitive string
template<typename char_type, typename traits_type, typename allocator_type>
inline std::basic_ostream<char_type, traits_type>&
operator<<(std::basic_ostream<char_type, traits_type>& os,
	       const std::basic_string<char_type, ci_char_traits, allocator_type>& str) {
    return std::__ostream_insert(os, str.data(), str.size());
}

__END_AKANTU__

#endif /* __AKANTU_CI_STRING_HH__ */
