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
			   membrane_transporters/ip3r.cpp
			   membrane_transporters/ryr.cpp
			   membrane_transporters/ryr2.cpp
			   membrane_transporters/serca.cpp
			   membrane_transporters/leak.cpp
			   membrane_transporters/pmca.cpp
			   membrane_transporters/ncx.cpp
			   membrane_transporters/vdcc_bg/vdcc_bg.cpp
			   membrane_transporters/vdcc_bg/vdcc_bg_userdata.cpp
			   membrane_transporters/mcu.cpp
			   membrane_transporters/mncx.cpp
			   stimulation/action_potential_train.cpp
			   grid_generation/bouton_generator/bouton_generator.cpp
   )


##message(STATUS "")
##message(STATUS "Info: Additional options for the neuro_collection (-DNC for setting options) plugin:")

if (MembranePotentialMapping)
    set(SOURCES ${SOURCES} membrane_transporters/vdcc_bg/vdcc_bg_vm2ug.cpp)

	if (MPMNEURON)
	    set(SOURCES ${SOURCES} membrane_transporters/vdcc_bg/vdcc_bg_neuron.cpp)
	endif(MPMNEURON)
endif (MembranePotentialMapping)



################################################################################
# The code below doesn't have to be changed (usually)
################################################################################
cmake_minimum_required(VERSION 2.6)

project(UG_PLUGIN_${pluginName})

# include the definitions and dependencies for ug-plugins.
include(${UG_ROOT_CMAKE_PATH}/ug_plugin_includes.cmake)

if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(${CMAKE_CURRENT_SOURCE_DIR} ${SOURCES})
	
	if (MembranePotentialMapping)
		# add definition to plugin_main.cpp to include Vm2Ug-dependent channels
		ExportSingleFileDefinitions(${CMAKE_CURRENT_SOURCE_DIR}/plugin_main.cpp VM2UG_ENABLED)
	
		if (MPNEURON)
			# add definition to plugin_main.cpp to include Neuron-dependent channels
			ExportSingleFileDefinitions(${CMAKE_CURRENT_SOURCE_DIR}/plugin_main.cpp NEURON_ENABLED)
		endif (MPNEURON)
		
	endif (MembranePotentialMapping)	
else(buildEmbeddedPlugins)
	# create a shared library from the sources and link it against ug4.
	add_library(${pluginName} SHARED ${SOURCES})
	
	if (MembranePotentialMapping)
		# add definition to plugin_main.cpp to include Vm2Ug-dependent channels
		set_source_files_properties(plugin_main.cpp PROPERTIES COMPILE_DEFINITIONS VM2UG_ENABLED)
		
		if (MPNEURON)
			# add definition to plugin_main.cpp to include Neuron-dependent channels
			set_source_files_properties(plugin_main.cpp PROPERTIES COMPILE_DEFINITIONS NEURON_ENABLED)
		endif (MPNEURON)
		
		# make plugin link against MPM
		set(linkLibraries MembranePotentialMapping ${linkLibraries})
	endif (MembranePotentialMapping)
	
	target_link_libraries (${pluginName} ug4 ${linkLibraries})
endif(buildEmbeddedPlugins)

