##
 # Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 #
 # Author: Markus Breit
 # Creation date: 2014-06-13
 #
 # This file is part of NeuroBox, which is based on UG4.
 #
 # NeuroBox and UG4 are free software: You can redistribute it and/or modify it
 # under the terms of the GNU Lesser General Public License version 3
 # (as published by the Free Software Foundation) with the following additional
 # attribution requirements (according to LGPL/GPL v3 §7):
 #
 # (1) The following notice must be displayed in the appropriate legal notices
 # of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 #
 # (2) The following notice must be displayed at a prominent place in the
 # terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 #
 # (3) The following bibliography is recommended for citation and must be
 # preserved in all covered files:
 # "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 #   parallel geometric multigrid solver on hierarchically distributed grids.
 #   Computing and visualization in science 16, 4 (2013), 151-164"
 # "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 #   flexible software system for simulating PDE based models on high performance
 #   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 # "Stepniewski, M., Breit, M., Hoffer, M. and Queisser, G.
 #   NeuroBox: computational mathematics in multiscale neuroscience.
 #   Computing and visualization in science (2019).
 # "Breit, M. et al. Anatomically detailed and large-scale simulations studying
 #   synapse loss and synchrony using NeuroBox. Front. Neuroanat. 10 (2016), 8"
 #
 # This program is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 # GNU Lesser General Public License for more details.
 ##

################################################################################
# NEURO collection                                                             #
#                                                                              #
# This plugin aims to collect and unify any neuro-related functionality that   #
# is developed for ug4 and that the author thinks of as usable for others.     #
#                                                                              #
################################################################################

# name of plugin, path to plugin (relative to ug's root path) and sources
set(pluginName neuro_collection)

cmake_minimum_required(VERSION 2.6)

project(UG_PLUGIN_${pluginName})

# include the definitions and dependencies for ug-plugins.
include(${UG_ROOT_CMAKE_PATH}/ug_plugin_includes.cmake)
message(STATUS "")
message(STATUS "Additional options for the neuro_collection plugin (-DNC for setting options):")

# tests
option(NCTestsuite "Build NC Testsuite" ${NCTestsuite})
message(STATUS "      Testsuite:   " ${NCTestsuite} " (options are: ON, OFF)")


set(SOURCES neuro_collection_plugin.cpp
            buffer_fv1.cpp
            user_flux_bnd_fv1.cpp
            membrane_transport_fv1.cpp
            membrane_transporters/membrane_transporter_interface.cpp
            membrane_transporters/hh.cpp
            membrane_transporters/hh_charges.cpp
            membrane_transporters/hh_species.cpp
            membrane_transporters/leakage_ohmic.cpp
            membrane_transporters/ip3r.cpp
            membrane_transporters/ryr.cpp
            membrane_transporters/ryr_instat.cpp
            membrane_transporters/ryr_discrete.cpp
            membrane_transporters/ryr_implicit.cpp
            membrane_transporters/ryr_linearized.cpp
            membrane_transporters/serca.cpp
            membrane_transporters/leak.cpp
            membrane_transporters/pmca.cpp
            membrane_transporters/ncx.cpp
            membrane_transporters/vdcc_bg/vdcc_bg.cpp
            membrane_transporters/vdcc_bg/vdcc_bg_userdata.cpp
            membrane_transporters/mcu.cpp
            membrane_transporters/mncx.cpp
            membrane_transporters/nmdar.cpp
            stimulation/action_potential_train.cpp
            grid_generation/bouton_generator.cpp
            grid_generation/dendrite_generator.cpp
            grid_generation/spine_generation.cpp
            util/axon_util.cpp
            util/hh_util.cpp
            util/misc_util.cpp
            util/neurite_axial_refinement_marker.cpp
   )
   
set(SOURCES_TEST unit_tests/tests.cpp)

## add experimental neurite projector impl (but not in VRL)   
if (NOT buildForVRL)
	set(SOURCES ${SOURCES}
				grid_generation/neurites_from_swc.cpp
				grid_generation/polygonal_mesh_from_txt.cpp
				test/test_neurite_proj.cpp
				test/neurite_refMarkAdjuster.cpp
				test/neurite_util_impl.cpp
				test/neurite_grid_generation_impl.cpp
				test/neurite_math_util_impl.cpp
				test/tetrahedralize_util.cpp
				test/consistency_util_impl.cpp
      )
endif (NOT buildForVRL)

if (MembranePotentialMapping)
	set(NC_WITH_MPM 1)
    set(SOURCES ${SOURCES} membrane_transporters/vdcc_bg/vdcc_bg_vm2ug.cpp)

	if (MPMNEURON)
		set(NC_WITH_NEURON 1)
	    set(SOURCES ${SOURCES} membrane_transporters/vdcc_bg/vdcc_bg_neuron.cpp)
	endif(MPMNEURON)
endif (MembranePotentialMapping)

if (cable_neuron)
	set(NC_WITH_CABLENEURON 1)
	set(SOURCES ${SOURCES} membrane_transporters/vdcc_bg/vdcc_bg_cableneuron.cpp
	                       hybrid_neuron_communicator.cpp
	                       membrane_transport_1d.cpp
	)
endif (cable_neuron)

if (Parmetis)
	set(NC_WITH_PARMETIS 1)
endif (Parmetis)


if(buildEmbeddedPlugins)
   set(NCTestsuite OFF)
endif(buildEmbeddedPlugins)

if(${NCTestsuite} STREQUAL "ON")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${UG_ROOT_PATH}/bin/)
    add_executable(NCTestsuite ${SOURCES_TEST})
endif(${NCTestsuite} STREQUAL "ON")


################################################################################
# The code below doesn't have to be changed (usually)
################################################################################
# This will copy the configurable file config_cmake.h to
# the final configuration header nc_config.h,
# setting the #defines as given by the current CMake variables.
configure_file(config_cmake.h ${CMAKE_BINARY_DIR}/nc_config.h)


if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(${CMAKE_CURRENT_SOURCE_DIR} ${SOURCES})
	
else(buildEmbeddedPlugins)
	# create a shared library from the sources and link it against ug4.
	add_library(${pluginName} SHARED ${SOURCES})
	
	if (MembranePotentialMapping)
		# make plugin link against MPM
		set(linkLibraries MembranePotentialMapping ${linkLibraries})
	endif (MembranePotentialMapping)
	
	if (cable_neuron)
		# make plugin link against cable_neuron
		set(linkLibraries cable_neuron ${linkLibraries})
	endif (cable_neuron)
	
	if (Parmetis)
		# make plugin link against Parmetis
		set(linkLibraries Parmetis ${linkLibraries})
	endif (Parmetis)
	
	target_link_libraries (${pluginName} ug4 ${linkLibraries})
	
	if(${NCTestsuite} STREQUAL "ON")
		target_link_libraries (NCTestsuite ${pluginName} ug4)
	endif(${NCTestsuite} STREQUAL "ON")
endif(buildEmbeddedPlugins)

