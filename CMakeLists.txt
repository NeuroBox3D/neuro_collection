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
set(pluginPath	plugins/experimental/neuro_collection)
set(SOURCES    neuro_collection_plugin.cpp
			   buffer_fv1.cpp
			   membrane_transport_fv1.cpp
			   membrane_transporters/membrane_transporter_interface.cpp
			   membrane_transporters/ip3r.cpp
			   membrane_transporters/ryr.cpp
			   membrane_transporters/serca.cpp
			   membrane_transporters/leak.cpp
			   membrane_transporters/pmca.cpp
			   membrane_transporters/ncx.cpp
			   membrane_transporters/vdcc_bg.cpp
			   stimulation/action_potential_train.cpp
   )


message(STATUS "")
message(STATUS "Info: Additional options for the neuro_collection (-DNC for setting options) plugin:")

if (NOT MembranePotentialMapping)
    message(FATAL_ERROR "ERROR: neuro_collection needs MembranePotentialMapping plugin. Please activate by setting -DMembranePotentialMapping=ON.")
endif (NOT MembranePotentialMapping)	

## ANN (requirement for vm2ug) #################################################

# get the path to annlib and set annlib include path and linker options
option(ANN "Sets the path in which ann shall be searched." OFF)

if(ANN)
    # the user can suply a ANN for his own compiled ANN   
    set(INTERNAL_ANN_PATH "${ANN}")
    message(STATUS "Info: ANN is ON, using user provided path ${ANN}")
else(ANN)
    UGInstall("ANN" INTERNAL_ANN_PATH)   	 
    message(STATUS "Info: ANN is ON, using uginstall provided path ${INTERNAL_ANN_PATH}")    
endif(ANN)


# find includes
FIND_PATH(ANN_INCLUDE_DIR ANN/ANN.h PATHS ${INTERNAL_ANN_PATH}/include)
IF (ANN_INCLUDE_DIR STREQUAL "ANN_INCLUDE_DIR-NOTFOUND")
	MESSAGE(FATAL_ERROR "ERROR: Could not find ANN header file ANN.h in ${INTERNAL_ANN_PATH}/include.")
ELSE (ANN_INCLUDE_DIR STREQUAL "ANN_INCLUDE_DIR-NOTFOUND")
	MESSAGE(STATUS "Info: ANN include:   ${ANN_INCLUDE_DIR}/ANN/ANN.h")
ENDIF (ANN_INCLUDE_DIR STREQUAL "ANN_INCLUDE_DIR-NOTFOUND")

# find libraries
SET(ANN_NAMES ${ANN_NAMES} ANN)
FIND_LIBRARY(ANN_LIBRARY NAMES ${ANN_NAMES} PATHS ${INTERNAL_ANN_PATH}/lib)
IF (ANN_LIBRARY STREQUAL "ANN_LIBRARY-NOTFOUND")
	MESSAGE(FATAL_ERROR "ERROR: Could not find ANN libraries.")
ELSE (ANN_LIBRARY STREQUAL "ANN_LIBRARY-NOTFOUND")
	MESSAGE(STATUS "Info: ANN library:   ${ANN_LIBRARY}")
ENDIF (ANN_LIBRARY STREQUAL "ANN_LIBRARY-NOTFOUND")

IF (ANN_LIBRARY AND ANN_INCLUDE_DIR)
	SET(ANN_LIBRARIES ${ANN_LIBRARY})
	SET(ANN_FOUND "YES")
ELSE (ANN_LIBRARY AND ANN_INCLUDE_DIR)
	SET(ANN_FOUND "NO")
ENDIF (ANN_LIBRARY AND ANN_INCLUDE_DIR)

IF (ANN_FOUND)
	include_directories(${ANN_INCLUDE_DIR})
	set(linkLibraries ${linkLibraries} ${ANN_LIBRARY})
endif(ANN_FOUND)

################################################################################

## NEURON interpreter ##########################################################
option(NCNEURON "Activates the NEURON interpreter." OFF)
if (NOT NCNEURON)
    message(STATUS "NEURON interpreter is deactivated. One can activate it with -DNCNEURON=ON.")
else(NOT NCNEURON)
    message(STATUS "NEURON interpreter is activated.")
    add_definitions(-DMPMNEURON)
    ExportDefinitions(-DMPMNEURON)
endif(NOT NCNEURON)

message(STATUS "")
################################################################################

################################################################################
# The code below doesn't have to be changed (usually)
################################################################################
cmake_minimum_required(VERSION 2.6)

project(UG_PLUGIN_${pluginName})

# include the definitions and dependencies for ug-plugins.
include("../../../cmake/ug_plugin_includes.cmake")

if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(${pluginPath} ${SOURCES})
	exportDependencies(${ANN_LIBRARY})
	exportIncludes(${ANN_INCLUDE_DIR})
else(buildEmbeddedPlugins)
	# create a shared library from the sources and link it against ug4.
	add_library(${pluginName} SHARED ${SOURCES})
	target_link_libraries (${pluginName} MembranePotentialMapping ug4 ${linkLibraries})
endif(buildEmbeddedPlugins)

