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
			   membrane_transporter_interface.cpp
			   ip3r.cpp
			   ryr.cpp
			   serca.cpp
			   leak.cpp
			   pmca.cpp
			   ncx.cpp
			   vdcc_bg.cpp
               dependent_neumann_boundary_fv1.cpp
   )


message(STATUS "")
message(STATUS "Info: Additional options for the neuro_collection (-DNC for setting options) plugin:")
## ANN (requirement for vm2ug) ###############################################

# get the path to annlib and set annlib include path and linker options

##############################################################################

# get the path to annlib and set annlib include path
option(ANNLIB "Sets the path in which annlib shall be searched." OFF)

if(NOT ANNLIB)
	message(FATAL_ERROR "ERROR: libann is required. Please set the option -DANNLIB=path_to_libann.")
endif(NOT ANNLIB)

# add include directoris for ANNLIB
find_library(ANNLIB_LIBS NAMES ANN PATHS ${ANNLIB}/lib)
include_directories(${ANNLIB}/include)

if(${ANNLIB_LIBS} STREQUAL "ANNLIB_LIBS-NOTFOUND")
	message(FATAL_ERROR "ERROR: Couldn't find ANNLIB in the sporecified path: " ${ANNLIB}/lib)
else(${ANNLIB_LIBS} STREQUAL "ANNLIB_LIBS-NOTFOUND")
	include_directories(${ANNLIB}/include)		
	message(STATUS "ANN Library API is on (by -ANNLIB=${ANNLIB}).")
endif(${ANNLIB_LIBS} STREQUAL "ANNLIB_LIBS-NOTFOUND")

######################
# NEURON interpreter #
######################
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
# The code below doesn't have to be changed (usually)
################################################################################
cmake_minimum_required(VERSION 2.6)

project(UG_PLUGIN_${pluginName})

# include the definitions and dependencies for ug-plugins.
include("../../../cmake/ug_plugin_includes.cmake")

if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(${pluginPath} ${SOURCES})
	exportDependencies(${ANNLIB_LIBS})
	exportIncludes(${ANNLIB}/include)
else(buildEmbeddedPlugins)
	# create a shared library from the sources and link it against ug4.
	add_library(${pluginName} SHARED ${SOURCES})
	target_link_libraries (${pluginName} ug4 ${linkLibraries})
endif(buildEmbeddedPlugins)

