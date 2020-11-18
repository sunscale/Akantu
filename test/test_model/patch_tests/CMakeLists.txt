#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author David Simon Kammer <david.kammer@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Oct 22 2010
# @date last modification: Thu Feb 08 2018
#
# @brief  configuration for patch tests
#
# @section LICENSE
#
# Copyright (©) 2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
# @section DESCRIPTION
#
#===============================================================================

add_subdirectory(data)

register_gtest_sources(SOURCES patch_test_linear_elastic_explicit.cc
  PACKAGE solid_mechanics
  FILES_TO_COPY data/material_check_stress_plane_strain.dat
                data/material_check_stress_plane_stress.dat
                data/material_check_stress_plane_strain_finite_deformation.dat
                data/material_check_stress_plane_stress_finite_deformation.dat)
register_gtest_sources(SOURCES patch_test_linear_elastic_implicit.cc
  PACKAGE solid_mechanics implicit
  FILES_TO_COPY data/material_check_stress_plane_strain.dat
                data/material_check_stress_plane_stress.dat
                data/material_check_stress_plane_strain_finite_deformation.dat
                data/material_check_stress_plane_stress_finite_deformation.dat)

register_gtest_sources(SOURCES patch_test_linear_anisotropic_explicit.cc
  PACKAGE solid_mechanics lapack
  FILES_TO_COPY
    data/material_anisotropic_1.dat
    data/material_anisotropic_2.dat
    data/material_anisotropic_3.dat
  )

register_gtest_sources(SOURCES patch_test_linear_anisotropic_implicit.cc
  PACKAGE solid_mechanics lapack implicit
  FILES_TO_COPY
    data/material_anisotropic_1.dat
    data/material_anisotropic_2.dat
    data/material_anisotropic_3.dat
  )

register_gtest_sources(SOURCES test_lumped_mass.cc
  PACKAGE solid_mechanics
  FILES_TO_COPY data/material_lumped.dat)

register_gtest_sources(
  SOURCES patch_test_linear_heat_transfer_explicit.cc
  FILES_TO_COPY data/heat_transfer_input.dat
  PACKAGE heat_transfer)

register_gtest_sources(
  SOURCES patch_test_linear_heat_transfer_static.cc
  FILES_TO_COPY data/heat_transfer_input.dat
  PACKAGE heat_transfer implicit)

register_gtest_sources(
  SOURCES patch_test_linear_heat_transfer_implicit.cc
  FILES_TO_COPY data/heat_transfer_input.dat
  PACKAGE heat_transfer implicit
  )

register_gtest_test(patch_test_linear
  FILES_TO_COPY ${PATCH_TEST_MESHES})

register_test(test_linear_elastic_explicit_python
  SCRIPT test_patch_linear_elastic_explicit.py
  PYTHON
  PACKAGE python_interface
  FILES_TO_COPY patch_test_linear_fixture.py
  FILES_TO_COPY patch_test_linear_solid_mechanics_fixture.py
  FILES_TO_COPY ${PATCH_TEST_MESHES})

register_test(test_linear_elastic_static_python
  SCRIPT test_patch_linear_elastic_static.py
  PYTHON
  PACKAGE python_interface implicit
  FILES_TO_COPY patch_test_linear_fixture.py
  FILES_TO_COPY patch_test_linear_solid_mechanics_fixture.py)

register_test(test_linear_anisotropic_explicit_python
  SCRIPT test_patch_linear_anisotropic_explicit.py
  PYTHON
  UNSTABLE
  PACKAGE python_interface lapack
  FILES_TO_COPY patch_test_linear_fixture.py
  FILES_TO_COPY patch_test_linear_solid_mechanics_fixture.py
  FILES_TO_COPY data/material_anisotropic_3.dat)


register_test(patch_test_linear_heat_transfer_explicit_python
  SCRIPT test_patch_linear_heat_transfer_explicit.py
  PYTHON
  PACKAGE python_interface heat_transfer
  FILES_TO_COPY patch_test_linear_fixture.py
  FILES_TO_COPY patch_test_linear_heat_transfer_fixture.py
  FILES_TO_COPY data/heat_transfer_input.dat)

register_test(patch_test_linear_heat_transfer_static_python
  SCRIPT test_patch_linear_heat_transfer_static.py
  PYTHON
  PACKAGE python_interface heat_transfer implicit
  FILES_TO_COPY patch_test_linear_fixture.py
  FILES_TO_COPY patch_test_linear_heat_transfer_fixture.py
  FILES_TO_COPY data/heat_transfer_input.dat)

register_test(patch_test_linear_heat_transfer_implicit_python
  SCRIPT test_patch_linear_heat_transfer_implicit.py
  PYTHON
  PACKAGE python_interface heat_transfer implicit
  FILES_TO_COPY patch_test_linear_fixture.py
  FILES_TO_COPY patch_test_linear_heat_transfer_fixture.py
  FILES_TO_COPY data/heat_transfer_input.dat)
