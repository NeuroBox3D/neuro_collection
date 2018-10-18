################################################################################
# NEURO collection
#
# This plugin aims to collect and unify any neuro-related functionality that is
# developed for ug4 and that the author thinks of as usable for others.
#
# Created on: 13.06.2014
#     Author: mbreit
#
################################################################################

# name of plugin, path to plugin (relative to ug's root path) and sources
set(pluginName	neuro_collection)
set(SOURCES neuro_collection_plugin.cpp
            buffer_fv1.cpp
            surface_marking.cpp
            membrane_transport_fv1.cpp
            membrane_transporters/membrane_transporter_interface.cpp
            membrane_transporters/hh.cpp
            membrane_transporters/hh_charges.cpp
            membrane_transporters/hh_species.cpp
            membrane_transporters/leakage_ohmic.cpp
            membrane_transporters/ip3r.cpp
            membrane_transporters/ryr.cpp
            membrane_transporters/ryr2.cpp
            membrane_transporters/ryr_implicit.cpp
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
            grid_generation/bouton_generator/bouton_generator.cpp
   )

## add experimental neurite projector impl (but not in VRL)   
#if (NOT buildForVRL)
#	set(SOURCES ${SOURCES}
#			   test/test_neurite_proj.cpp
#			   test/neurite_refMarkAdjuster.cpp
#      )
#endif (NOT buildForVRL)


##message(STATUS "")
##message(STATUS "Additional options for the neuro_collection plugin (-DNC for setting options):")

if (MembranePotentialMapping)
	set(NC_WITH_VM2UG 1)
    set(SOURCES ${SOURCES} membrane_transporters/vdcc_bg/vdcc_bg_vm2ug.cpp)

	if (MPMNEURON)
		set(NC_WITH_NEURON 1)
	    set(SOURCES ${SOURCES} membrane_transporters/vdcc_bg/vdcc_bg_neuron.cpp)
	endif(MPMNEURON)
endif (MembranePotentialMapping)

if (cable_neuron)
	set(NC_WITH_CABLENEURON 1)
	set(SOURCES ${SOURCES} membrane_transporters/vdcc_bg/vdcc_bg_cableneuron.cpp
	                       hybrid_neuron_communicator.cpp)
endif (cable_neuron)


################################################################################
# The code below doesn't have to be changed (usually)
################################################################################
cmake_minimum_required(VERSION 2.6)

project(UG_PLUGIN_${pluginName})

# include the definitions and dependencies for ug-plugins.
include(${UG_ROOT_CMAKE_PATH}/ug_plugin_includes.cmake)


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
	
	## temporarily commented out until split synapse handler replaces synapse handler
	if (cable_neuron)
		# make plugin link against cable_neuron
		set(linkLibraries cable_neuron ${linkLibraries})
	endif (cable_neuron)
	
	target_link_libraries (${pluginName} ug4 ${linkLibraries})
endif(buildEmbeddedPlugins)

