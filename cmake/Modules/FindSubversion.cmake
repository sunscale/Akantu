#===============================================================================
# @file   FindSubversion.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Sun Oct 19 2014
#
# @brief  find_package for subversion
#
# @section LICENSE
#
# Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
# (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free  software: you can redistribute it and/or  modify it under the
# terms  of the  GNU Lesser  General Public  License as  published by  the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
# details.
#
# You should  have received  a copy  of the GNU  Lesser General  Public License
# along with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================

# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

find_program(Subversion_SVN_EXECUTABLE svn
  DOC "subversion command line client")
mark_as_advanced(Subversion_SVN_EXECUTABLE)

if(Subversion_SVN_EXECUTABLE)
  # the subversion commands should be executed with the C locale, otherwise
  # the message (which are parsed) may be translated, Alex
  set(_Subversion_SAVED_LC_ALL "$ENV{LC_ALL}")
  set(ENV{LC_ALL} C)

  execute_process(COMMAND ${Subversion_SVN_EXECUTABLE} --version
    OUTPUT_VARIABLE Subversion_VERSION_SVN
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  # restore the previous LC_ALL
  set(ENV{LC_ALL} ${_Subversion_SAVED_LC_ALL})

  string(REGEX REPLACE "^(.*\n)?svn, version ([.0-9]+).*"
    "\\2" Subversion_VERSION_SVN "${Subversion_VERSION_SVN}")

  macro(Subversion_WC_INFO dir prefix)
    # the subversion commands should be executed with the C locale, otherwise
    # the message (which are parsed) may be translated, Alex
    set(_Subversion_SAVED_LC_ALL "$ENV{LC_ALL}")
    set(ENV{LC_ALL} C)

    execute_process(COMMAND ${Subversion_SVN_EXECUTABLE} info ${dir}
      OUTPUT_VARIABLE ${prefix}_WC_INFO
      ERROR_VARIABLE Subversion_svn_info_error
      RESULT_VARIABLE Subversion_svn_info_result
      OUTPUT_STRIP_TRAILING_WHITESPACE)

    if(NOT ${Subversion_svn_info_result} EQUAL 0)
      set(${prefix}_WC_FOUND FALSE)
      if(NOT "${ARGV2}" STREQUAL "ERROR_QUIET")
        message(SEND_ERROR "Command \"${Subversion_SVN_EXECUTABLE} info ${dir}\" failed with output:\n${Subversion_svn_info_error}")
      endif(NOT "${ARGV2}" STREQUAL "ERROR_QUIET")
     else(NOT ${Subversion_svn_info_result} EQUAL 0)
      set(${prefix}_WC_FOUND TRUE)

      string(REGEX REPLACE "^(.*\n)?URL: ([^\n]+).*"
        "\\2" ${prefix}_WC_URL "${${prefix}_WC_INFO}")
      string(REGEX REPLACE "^(.*\n)?Repository Root: ([^\n]+).*"
        "\\2" ${prefix}_WC_ROOT "${${prefix}_WC_INFO}")
      string(REGEX REPLACE "^(.*\n)?Revision: ([^\n]+).*"
        "\\2" ${prefix}_WC_REVISION "${${prefix}_WC_INFO}")
      string(REGEX REPLACE "^(.*\n)?Last Changed Author: ([^\n]+).*"
        "\\2" ${prefix}_WC_LAST_CHANGED_AUTHOR "${${prefix}_WC_INFO}")
      string(REGEX REPLACE "^(.*\n)?Last Changed Rev: ([^\n]+).*"
        "\\2" ${prefix}_WC_LAST_CHANGED_REV "${${prefix}_WC_INFO}")
      string(REGEX REPLACE "^(.*\n)?Last Changed Date: ([^\n]+).*"
        "\\2" ${prefix}_WC_LAST_CHANGED_DATE "${${prefix}_WC_INFO}")

    endif()

    # restore the previous LC_ALL
    set(ENV{LC_ALL} ${_Subversion_SAVED_LC_ALL})

  endmacro()

  macro(Subversion_WC_LOG dir prefix)
    # This macro can block if the certificate is not signed:
    # svn ask you to accept the certificate and wait for your answer
    # This macro requires a svn server network access (Internet most of the time)
    # and can also be slow since it access the svn server
    execute_process(COMMAND
      ${Subversion_SVN_EXECUTABLE} --non-interactive log -r BASE ${dir}
      OUTPUT_VARIABLE ${prefix}_LAST_CHANGED_LOG
      ERROR_VARIABLE Subversion_svn_log_error
      RESULT_VARIABLE Subversion_svn_log_result
      OUTPUT_STRIP_TRAILING_WHITESPACE)

    if(NOT ${Subversion_svn_log_result} EQUAL 0)
      message(SEND_ERROR "Command \"${Subversion_SVN_EXECUTABLE} log -r BASE ${dir}\" failed with output:\n${Subversion_svn_log_error}")
    endif()
  endmacro()

endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Subversion REQUIRED_VARS Subversion_SVN_EXECUTABLE
                                             VERSION_VAR Subversion_VERSION_SVN )

# for compatibility
set(Subversion_FOUND ${SUBVERSION_FOUND})
set(Subversion_SVN_FOUND ${SUBVERSION_FOUND})
