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
set(SOURCES    neuro_collection_plugin.cpp
			   buffer_fv1.cpp
			   membrane_transport_fv1.cpp
			   membrane_transporters/membrane_transporter_interface.cpp
			   membrane_transporters/hh.cpp
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
if (NOT buildForVRL)
	set(SOURCES ${SOURCES}
			   test/test_neurite_proj.cpp
			   test/neurite_refMarkAdjuster.cpp
      )
endif (NOT buildForVRL)

# use QHULL if available
option(QHULL "Provide installation path to qhull implementation." OFF)
find_library(QHULL_LIB_CPP NAMES qhullcpp PATHS ${QHULL}/build/)
find_library(QHULL_LIB_STATIC NAMES qhullstatic PATHS ${QHULL}/build/)

if(${QHULL_LIB_CPP} STREQUAL "QHULL_LIB_CPP-NOTFOUND" OR ${QHULL_LIB_STATIC} STREQUAL "QHULL_LIB_CPP-NOTFOUND")
	message(WARNING "WARN: Couldn't find Qhull.org in the specified path: " ${QHULL}/build/libqhull/)
else(${QHULL_LIB_CPP} STREQUAL "QHULL_LIB_CPP-NOTFOUND" OR ${QHULL_LIB_STATIC} STREQUAL "QHULL_LIB_CPP-NOTFOUND")
	message(STATUS "Qhull.org library API is on (by -QHULL=${QHULL}).")
    message(STATUS "")
    set(NC_WITH_QHULL 1)
endif(${QHULL_LIB_CPP} STREQUAL "QHULL_LIB_CPP-NOTFOUND" OR ${QHULL_LIB_STATIC} STREQUAL "QHULL_LIB_CPP-NOTFOUND")
message(STATUS "${QHULL_LIB_CPP}")
message(STATUS "${QHULL_LIB_STATIC}")

include_directories(${QHULL}/src/libqhullcpp/ ${QHULL}/src/libqhull/ ${QHULL}/src/libqhull_r/ ${QHULL}/src/)

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
# the final configuration header config.h,
# setting the #defines as given by the current CMake variables.
configure_file(config_cmake.h config.h)


if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(${CMAKE_CURRENT_SOURCE_DIR} ${SOURCES})

	# add include path for config.h
	ExportIncludes(${CMAKE_CURRENT_BINARY_DIR})
	
else(buildEmbeddedPlugins)
	# create a shared library from the sources and link it against ug4.
	add_library(${pluginName} SHARED ${SOURCES})
	
	# add include path for tetgen_config.h
	include_directories(${CMAKE_CURRENT_BINARY_DIR})
	
	if (MembranePotentialMapping)
		# make plugin link against MPM
		set(linkLibraries MembranePotentialMapping ${linkLibraries} ${QHULL_LIB_CPP} ${QHULL_LIB_STATIC}) # /Users/stephan/Code/git/qhull/build/libqhullstatic.a)
	endif (MembranePotentialMapping)
	
	## temporarily commented out until split synapse handler replaces synapse handler
	if (cable_neuron)
		# make plugin link against cable_neuron
		set(linkLibraries cable_neuron ${linkLibraries})
	endif (cable_neuron)
	
	target_link_libraries (${pluginName} ug4 ${linkLibraries} ${QHULL_LIB_CPP} ${QHULL_LIB_STATIC}) # /Users/stephan/Code/git/qhull/build/libqhullstatic.a)
endif(buildEmbeddedPlugins)

