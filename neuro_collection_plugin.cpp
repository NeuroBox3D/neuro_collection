/*
 * NEURO collection
 *
 * This plugin aims to collect and unify any neuro-related functionality that is
 * developed for ug4 and that the author thinks of as usable for others.
 *
 *
 *  Created on: 13.06.2014
 *      Author: mbreit
 */

#include "bridge/util.h"

// replace this with util_domain_dependent.h or util_algebra_dependent.h
// to speed up compilation time
#include "bridge/util_domain_algebra_dependent.h"
#include "bridge/util_domain_dependent.h"

#include "lib_grid/global_attachments.h"  // for GlobalAttachments::declare_attachment

// configuration file for compile options
#include "nc_config.h"

#include "buffer_fv1.h"
#include "surface_marking.h"                                 // for SurfaceMarking
#include "membrane_transport_fv1.h"
#include "user_flux_bnd_fv1.h"
#include "membrane_transporters/membrane_transporter_interface.h"
#include "membrane_transporters/hh.h"
#include "membrane_transporters/hh_charges.h"
#include "membrane_transporters/hh_species.h"
#include "membrane_transporters/leakage_ohmic.h"
#include "membrane_transporters/ip3r.h"
#include "membrane_transporters/ryr.h"
#include "membrane_transporters/ryr2.h"
#include "membrane_transporters/ryr_discrete.h"
#include "membrane_transporters/ryr_implicit.h"
#include "membrane_transporters/serca.h"
#include "membrane_transporters/leak.h"
#include "membrane_transporters/pmca.h"
#include "membrane_transporters/ncx.h"
#include "membrane_transporters/vdcc_bg/vdcc_bg.h"
#include "membrane_transporters/vdcc_bg/vdcc_bg_userdata.h"

#ifdef NC_WITH_CABLENEURON
	#include "hybrid_neuron_communicator.h"
	#include "hybrid_synapse_current_assembler.h"
    #include "membrane_transporters/vdcc_bg/vdcc_bg_cableneuron.h"
#endif

#ifdef NC_WITH_VM2UG
	#include "membrane_transporters/vdcc_bg/vdcc_bg_vm2ug.h"
	#ifdef NC_WITH_NEURON
		#include "membrane_transporters/vdcc_bg/vdcc_bg_neuron.h"
	#endif
#endif

#include "membrane_transporters/mcu.h"
#include "membrane_transporters/mncx.h"
#include "membrane_transporters/nmdar.h"
#include "stimulation/action_potential_train.h"
#include "grid_generation/bouton_generator.h"
#include "grid_generation/spine_generation.h"

#include "lib_grid/refinement/projectors/neurite_projector.h"
#include "test/test_neurite_proj.h"

#include "util/measurement.h"
#include "util/ca_wave_util.h"
#include "util/axon_util.h"
#include "util/misc_util.h"
#include "util/neurite_axial_refinement_marker.h"
#include "lib_disc/function_spaces/grid_function.h"


using namespace std;
using namespace ug::bridge;

namespace ug{
namespace neuro_collection{

/** 
 *  \defgroup plugin_neuro_collection Plugin neuro_collection
 *  \ingroup plugins_experimental
 *  This plugin aims to collect and unify any neuro-related functionality that is
 *  developed for ug4 and that the author thinks of as usable for others.
 *  \{
 */

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

/**
 * Function called for the registration of Domain and Algebra dependent parts.
 * All Functions and Classes depending on both Domain and Algebra
 * are to be placed here when registering. The method is called for all
 * available Domain and Algebra types, based on the current build options.
 *
 * @param reg		registry
 * @param grp		group for sorting of functionality
 */
template <typename TDomain, typename TAlgebra>
static void DomainAlgebra(Registry& reg, string grp)
{
	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

	typedef GridFunction<TDomain, TAlgebra> TGridFunction;

	// extra commands for this plugin
	reg.add_function("take_measurement", static_cast<number (*)(SmartPtr<TGridFunction>, const number, const char*, const char*, const char*)>(&takeMeasurement<GridFunction<TDomain, TAlgebra> >), grp.c_str(),
					 "", "solution#time#subset names#function names#output file name",
					 "outputs average values of unknowns on subsets");

	// calcium wave examination functions
	reg.add_function("max_ryr_flux_density", &maxRyRFluxDensity<TGridFunction, RyRImplicit<TDomain> >, grp.c_str(),
						 "", "solution # function names for ca_cyt, ca_er, c1, c2 as c-style string #"
							 "RyR-carrying membrane subset names as c-style string # RyR channel",
						 "maximal flux density through RyR channel (mol/(m^2*s))");
	reg.add_function("wave_front_x", &waveFrontX<TGridFunction>, grp.c_str(),
					 "", "solution # function names for c1, c2 as c-style string #"
						 "RyR-carrying membrane subset names as c-style string # threshold open probability",
					 "rightmost vertex where threshold value is exceeded");

	// WaveProfileExporter
	{
		typedef WaveProfileExporter<TDomain, TAlgebra> T;
		string name = string("WaveProfileExporter").append(suffix);
		reg.add_class_<T>(name, grp)
			.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >,
				const char*, const char*, const std::string&)>
				("approximation space # function names (comma-separated c-string) # "
					"subset names (comma-separated c-string) # file base name")
			.add_method("exportWaveProfileX", &T::exportWaveProfileX, "", "", "")
			.set_construct_as_smart_pointer(true);

		reg.add_class_to_group(name, "WaveProfileExporter", tag);
	}

	// scaling of dimless solution vectors
	{
		reg.add_function("scale_dimless_vector", &scale_dimless_vector<TGridFunction>,
			grp.c_str(), "", "scaled output vector#dimless input vector#vector of scaling factors for each function",
			"Scales the dimensionless input vector using the given scaling factors for each function and writes "
			"the result to the output vector");
	}

	reg.add_function("take_measurement", static_cast<number (*)(SmartPtr<TGridFunction>, const number, const char*, const char*, const char*, const char*)>(&takeMeasurement<GridFunction<TDomain, TAlgebra> >), grp.c_str(),
					 "", "solution#time#subset names#function names#output file name#output file extension",
					 "outputs average values of unknowns on subsets");


	// export all template realizations of RyRImplicit::calculate_steady_state()
	{
		typedef RyRImplicit<TDomain> T;
		ClassGroupDesc* cgd = reg.get_class_group(std::string("RyRImplicit"));
		size_t numClasses = cgd->num_classes();
		size_t i = 0;
		for (; i < numClasses; ++i)
		{
			std::string classTag = cgd->get_class_tag(i);
			if (classTag == GetDomainTag<TDomain>())
			{
				ExportedClass<T>* expClass = dynamic_cast<ExportedClass<T>* >(cgd->get_class(i));
				UG_COND_THROW(!expClass, "Exported class can not be cast to the correct type.");

				expClass->add_method("calculate_steady_state",
					&T::template calculate_steady_state<typename TAlgebra::vector_type>, "", "", "");

				break;
			}
		}
		UG_COND_THROW(i == numClasses, "No class with domain tag '" << GetDomainTag<TDomain>()
			<< "' found in RyRImplicit class group to add algebra-dependent functionality to.");

		// same again for 1d special case
		typedef RyRImplicit_1drotsym<TDomain> T1;
		cgd = reg.get_class_group(std::string("RyRImplicit_1drotsym"));
		numClasses = cgd->num_classes();
		i = 0;
		for (; i < numClasses; ++i)
		{
			std::string classTag = cgd->get_class_tag(i);
			if (classTag == GetDomainTag<TDomain>())
			{
				ExportedClass<T1>* expClass = dynamic_cast<ExportedClass<T1>* >(cgd->get_class(i));
				UG_COND_THROW(!expClass, "Exported class can not be cast to the correct type.");

				expClass->add_method("calculate_steady_state",
					&T1::template calculate_steady_state<typename TAlgebra::vector_type>, "", "", "");

				break;
			}
		}
		UG_COND_THROW(i == numClasses, "No class with domain tag '" << GetDomainTag<TDomain>()
			<< "' found in RyRImplicit_1drotsym class group to add algebra-dependent functionality to.");
	}
}

/**
 * Function called for the registration of Domain dependent parts.
 * All Functions and Classes depending on the Domain
 * are to be placed here when registering. The method is called for all
 * available Domain types, based on the current build options.
 *
 * @param reg		registry
 * @param grp		group for sorting of functionality
 */
template <typename TDomain>
static void Domain(Registry& reg, string grp)
{
	const int dim = TDomain::dim;
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();


	// implementation of buffering reaction disc
	{
		typedef BufferFV1<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("BufferFV1").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(const char*)>("Subset(s)")
			.add_method("set_num_reactions", &T::set_num_reactions, "", "number of reactions | default | value=1",
				"set number of reactions about to be added")
			.add_method("add_reaction", static_cast<void (T::*)
				 (const char*, const char*, number, number, number)>(&T::add_reaction), "",
				  "buffering substance | selection | value=[\"clb\"] # "
				  "buffered substance | selection | value=[\"ca_cyt\"] # "
				  "total buffer | default |	value=160.0e-6D # "
				  "association rate | default | value=27.0e06D # "
				  "dissociation rate | default | value=19.0D",
				  "add a new reaction definition")
#ifdef UG_FOR_LUA
			.add_method("add_reaction", static_cast<void (T::*)
				 (const char*, const char*, number, const char*, const char*)>(&T::add_reaction), "",
				  "buffering substance | selection | value=[\"clb\"] # "
				  "buffered substance | selection | value=[\"ca_cyt\"] # "
				  "total buffer | default |	value=160.0e-6D # "
				  "association rate | default | value=27.0e06D # "
				  "dissociation rate | default | value=19.0D",
				  "add a new reaction definition")
#endif
			.add_method("add_reaction", static_cast<void (T::*)
				(const char*, const char*, SmartPtr<CplUserData<number, dim> >,
				 SmartPtr<CplUserData<number, dim> >, SmartPtr<CplUserData<number, dim> >)>(&T::add_reaction), "",
				 "buffering substance | selection | value=[\"clb\"] # "
				 "buffered substance | selection | value=[\"ca_cyt\"] # "
				 "total buffer | default |	value=160.0e-6D # "
				 "association rate | default | value=27.0e06D # "
				 "dissociation rate | default | value=19.0D",
				 "add a new reaction definition")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "BufferFV1", tag);
	}

	// SurfaceMarking
	{
		typedef SurfaceMarking<TDomain> T;
		typedef IElementMarkingStrategy<TDomain> TBase;
		string name = string("SurfaceMarking").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(ConstSmartPtr<TDomain>)>("domain")
			.template add_constructor<void (*)(number, size_t)>("tolerance#maximal level of refinement")
			.template add_constructor<void (*)(ConstSmartPtr<TDomain>, number, size_t)>("domain#tolerance#maximal level of refinement")
			.add_method("set_tolerance", &T::set_tolerance, "", "tolerance", "", "")
			.add_method("set_max_level", &T::set_max_level, "", "maximal refinement level", "", "")
			.add_method("add_surface", static_cast<void (T::*) (int, int)>(&T::add_surface),
				"", "surface subset index # adjacent element subset index to be refined", "", "")
			.add_method("remove_surface", static_cast<void (T::*) (int, int)>(&T::remove_surface),
				"", "surface subset index # adjacent element subset index to be refined", "", "")
			.add_method("add_surface", static_cast<void (T::*) (const std::string&, const std::string&)>(&T::add_surface),
				"", "surface subset name # adjacent element subset name to be refined", "", "")
			.add_method("remove_surface", static_cast<void (T::*) (const std::string&, const std::string&)>(&T::remove_surface),
				"", "surface subset name # adjacent element subset name to be refined", "", "")
			//.add_method("add_interface", &T::add_interface, "", "interfaces to take into account", "")
			.add_method("mark_without_error", &T::mark_without_error, "", "refiner#approximation space", "", "")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "SurfaceMarking", tag);
	}

	// implementation of two-sided membrane transport systems
	{
		typedef MembraneTransportFV1<TDomain> T;
		typedef FV1InnerBoundaryElemDisc<TDomain> TBase;
		string name = string("MembraneTransportFV1").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*, SmartPtr<IMembraneTransporter>)>("Subset(s) as comma-separated c-string#MembraneTransporter")
			.template add_constructor<void (*)(const std::vector<std::string>&, SmartPtr<IMembraneTransporter>)>("Subset(s) as vector#MembraneTransporter")
			.add_method("set_density_function", static_cast<void (T::*) (const number)> (&T::set_density_function),
						"", "", "add a constant density")
#ifdef UG_FOR_LUA
			.add_method("set_density_function", static_cast<void (T::*) (const char*)> (&T::set_density_function),
						"", "", "add a density function")
#endif
			.add_method("set_density_function", static_cast<void (T::*) (SmartPtr<CplUserData<number,dim> >)>
					(&T::set_density_function), "", "", "add a density function")
			.add_method("set_membrane_transporter", &T::set_membrane_transporter, "", "", "sets the membrane transport mechanism")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "MembraneTransportFV1", tag);
	}

	// implementation of two-sided membrane transport systems (1d "cable", fcts const in radius and angle)
	{
		typedef MembraneTransport1d<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("MembraneTransport1d").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*, SmartPtr<IMembraneTransporter>)>("Subset(s) as comma-separated c-string#MembraneTransporter")
			.template add_constructor<void (*)(const std::vector<std::string>&, SmartPtr<IMembraneTransporter>)>("Subset(s) as vector#MembraneTransporter")
			.add_method("set_density_function", static_cast<void (T::*) (const number)> (&T::set_density_function),
						"", "", "add a constant density")
#ifdef UG_FOR_LUA
			.add_method("set_density_function", static_cast<void (T::*) (const char*)> (&T::set_density_function),
						"", "", "add a density function")
#endif
			.add_method("set_density_function", static_cast<void (T::*) (SmartPtr<CplUserData<number,dim> >)>
					(&T::set_density_function), "", "", "add a density function")
			.add_method("set_radius", &T::set_radius, "", "", "sets the radius the membrane is located at")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "MembraneTransport1d", tag);
	}

	// user flux boundary
	{
		typedef UserFluxBoundaryFV1<TDomain> T;
		typedef FV1InnerBoundaryElemDisc<TDomain> TBase;
		string name = string("UserFluxBoundaryFV1").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(const char*, const char*)>("Function(s) as comma-separated c-string#Subset(s) as comma-separated c-string")
			.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s) as vector#Subset(s) as vector")
			.add_method("set_flux_function", static_cast<void (T::*) (SmartPtr<CplUserData<number, dim> >)> (&T::set_flux_function),
					"", "", "add a flux density function")
			.add_method("set_flux_function", static_cast<void (T::*) (number)> (&T::set_flux_function),
					"", "", "add a flux density function")
			.add_method("set_flux_function", static_cast<void (T::*) (const char*)> (&T::set_flux_function),
					"", "", "add a flux density function")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "UserFluxBoundaryFV1", tag);
	}

	// Hodgkin-Huxley channels
	{
		typedef HH<TDomain> T;
		typedef IMembraneTransporter TBase1;
		typedef IElemDisc<TDomain> TBase2;
		std::string name = std::string("HH").append(suffix);
		reg.add_class_<T, TBase1, TBase2>(name, grp)
			.template add_constructor<void (*)(const char*, const char*)>
				("Functions as comma-separated string with the following order: "
				 "\"inner potential\", \"outer potential\", \"gating param n\", \"gating param m\", \"gating param h\" # "
				 "subsets as comma-separated string")
			.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>
				("Function vector with the following order: "
				 "\"inner potential\", \"outer potential\", \"gating param n\", \"gating param m\", \"gating param h\" # "
				 "subsets vector")
			.add_method("set_conductances", &T::set_conductances, "", "g_K#g_Na", "")
			.add_method("set_reversal_potentials", &T::set_reversal_potentials, "", "E_K#E_Na", "")
			.add_method("set_reference_time", &T::set_reference_time, "", "reference time (in units of s)", "")
			.add_method("use_exact_gating_mode", &T::use_exact_gating_mode, "", "time step size", "")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "HH", tag);
	}

	// Hodgkin-Huxley channels with charge densities as source / target
	{
		typedef HHCharges<TDomain> T;
		typedef IMembraneTransporter TBase1;
		typedef IElemDisc<TDomain> TBase2;
		std::string name = std::string("HHCharges").append(suffix);
		reg.add_class_<T, TBase1, TBase2>(name, grp)
			.template add_constructor<void (*)(const char*, const char*)>
				("Functions as comma-separated string with the following order: "
				 "\"inner charge density\", \"outer charge density\", \"inner potential\", \"outer potential\","
				 "\"gating param n\", \"gating param m\", \"gating param h\" # "
				 "subsets as comma-separated string")
			.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>
				("Function vector with the following order: "
				 "\"inner charge density\", \"outer charge density\", \"inner potential\", \"outer potential\","
				 "\"gating param n\", \"gating param m\", \"gating param h\" # "
				 "subsets vector")
			.add_method("set_conductances", &T::set_conductances, "", "g_K#g_Na", "")
			.add_method("set_reversal_potentials", &T::set_reversal_potentials, "", "E_K#E_Na", "")
			.add_method("set_reference_time", &T::set_reference_time, "", "reference time (in units of s)", "")
			.add_method("use_exact_gating_mode", &T::use_exact_gating_mode, "", "time step size", "")
			.add_method("use_gating_explicit_current_mode", &T::use_gating_explicit_current_mode, "", "time step size", "")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "HHCharges", tag);
	}

	// Hodgkin-Huxley channels with Na and K unknowns
	{
		typedef HHSpecies<TDomain> T;
		typedef IMembraneTransporter TBase;
		std::string name = std::string("HHSpecies").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(const char*, const char*, ConstSmartPtr<ISubsetHandler>)>
				("Functions as comma-separated string with the following order: "
				 "\"inner [K+] \", \"outer [K+]\", \"inner [Na+] \", \"outer [Na+]\", "
				 "\"inner potential\", \"outer potential\" # "
				 "subsets as comma-separated string # SubsetHandler")
			.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&, ConstSmartPtr<ISubsetHandler>)>
				("Function vector with the following order: "
				 "\"inner charge density\", \"outer charge density\", \"inner potential\", \"outer potential\","
				 "\"gating param n\", \"gating param m\", \"gating param h\" # "
				 "subsets vector # SubsetHandler")
			.add_method("set_conductances", &T::set_conductances, "", "g_K#g_Na", "")
			.add_method("set_reversal_potentials", &T::set_reversal_potentials, "", "E_K#E_Na", "")
			.add_method("set_reference_time", &T::set_reference_time, "", "reference time (in units of s)", "")
			.add_method("set_temperature", &T::set_temperature, "", "temperature (in K)", "")
			.add_method("use_exact_gating_mode", &T::use_exact_gating_mode, "", "time step size", "")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "HHSpecies", tag);
	}

	// RyR2 (time-dep. RyR implementation)
	{
		typedef RyR2<TDomain> T;
		typedef IMembraneTransporter TBase;
		std::string name = std::string("RyR2").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(const char*, const char*, SmartPtr<ApproximationSpace<TDomain> >)>
				("Functions as comma-separated string with the order: "
				 "{\"cytosolic calcium\", \"endoplasmic calcium\"} # "
				 "subsets as comma-separated string # approximation space")
			.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&, SmartPtr<ApproximationSpace<TDomain> >)>
				("Function vector with the order: "
				 "{\"cytosolic calcium\", \"endoplasmic calcium\"} # "
				 "subsets vector, approximation space")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "RyR2", tag);
	}

	// fully implicit RyR
	{
		typedef RyRImplicit<TDomain> T;
		typedef IMembraneTransporter TBase1;
		typedef IElemDisc<TDomain> TBase2;
		std::string name = std::string("RyRImplicit").append(suffix);
		reg.add_class_<T, TBase1, TBase2>(name, grp)
			.template add_constructor<void (*)(const char*, const char*)>
				("Functions as comma-separated string with the order: "
				 "{\"cytosolic calcium\", \"endoplasmic calcium\", \"O2 channel state\", \"C1 channel state\", \"C2 channel state\"} # "
				 "subsets as comma-separated string")
			.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>
				("Function vector with the order: "
				 "{\"cytosolic calcium\", \"endoplasmic calcium\", \"O2 channel state\", \"C1 channel state\", \"C2 channel state\"} # "
				 "subsets vector,")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "RyRImplicit", tag);
	}

	// fully implicit RyR (special case 1d, rotationally symmetric "cable")
	{
		typedef RyRImplicit_1drotsym<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		std::string name = std::string("RyRImplicit_1drotsym").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(const char*, const char*)>
				("Functions as comma-separated string with the order: "
				 "{\"cytosolic calcium\", \"endoplasmic calcium\", \"O2 channel state\", \"C1 channel state\", \"C2 channel state\"} # "
				 "subsets as comma-separated string")
			.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>
				("Function vector with the order: "
				 "{\"cytosolic calcium\", \"endoplasmic calcium\", \"O2 channel state\", \"C1 channel state\", \"C2 channel state\"} # "
				 "subsets vector,")
			.add_method("set_calcium_scale", &T::set_calcium_scale, "", "cytosolic calcium scale", "")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "RyRImplicit_1drotsym", tag);
	}

	// VDCC base type
	{
		typedef VDCC_BG<TDomain> T;
		typedef IMembraneTransporter TBase;
		std::string name = std::string("VDCC_BG").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_method("set_channel_type_N", &T::template set_channel_type<T::BG_Ntype>,
						"", "", "set the channel type to N")
			.add_method("set_channel_type_L", &T::template set_channel_type<T::BG_Ltype>,
						"", "", "set the channel type to L")
			.add_method("set_channel_type_T", &T::template set_channel_type<T::BG_Ttype>,
						"", "", "set the channel type to T")
			.add_method("init", &T::init, "", "time", "initialize the Borg-Graham object");
		reg.add_class_to_group(name, "VDCC_BG", tag);
	}

	// VDCC with UserData
	{
		typedef VDCC_BG_UserData<TDomain> T;
		typedef VDCC_BG<TDomain> TBase;
		std::string name = std::string("VDCC_BG_UserData").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&, SmartPtr<ApproximationSpace<TDomain> >)>
				("function(s) as vector#subset(s) as vector#approxSpace")
			.template add_constructor<void (*)(const char*, const char*, SmartPtr<ApproximationSpace<TDomain> >)>
				("function(s) as comma-separated c-string#subset(s) as comma-separated c-string#approxSpace")
			.add_method("set_potential_function", static_cast<void (T::*) (const number)> (&T::set_potential_function),
						"", "", "add a potential function")
			.add_method("set_potential_function", static_cast<void (T::*) (const char*)> (&T::set_potential_function),
						"", "", "add a potential function")
			.add_method("set_potential_function", static_cast<void (T::*) (SmartPtr<CplUserData<number, dim> >)> (&T::set_potential_function),
						"", "", "add a potential function")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "VDCC_BG_UserData", tag);
	}

#ifdef NC_WITH_CABLENEURON
	// VDCC with cable_neuron
    {
        typedef VDCC_BG_CN<TDomain> T;
        typedef VDCC_BG<TDomain> TBase;
        std::string name = std::string("VDCC_BG_CN").append(suffix);
        reg.add_class_<T, TBase>(name, grp)
            .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&,
                SmartPtr<ApproximationSpace<TDomain> >, SmartPtr<ApproximationSpace<TDomain> >, const std::string&)>
                ("function(s) as vector#subset(s) as vector#approxSpace 1d#approxSpace 3d#potential function name")
            .add_method("set_domain_disc_1d", &T::set_domain_disc_1d, "", "domainDisc",
                "Set the 1d cable domain discretization.")
			.add_method("set_cable_disc", &T::set_cable_disc, "", "cableDisc",
				"Set the 1d cable element discretization.")
			.add_method("set_3d_neuron_ids", &T::set_3d_neuron_ids, "", "neuron ids as vector",
				"Set the 3d represented neuron IDs.")
            .add_method("set_initial_values", &T::set_initial_values, "", "initial value(s) as vector",
                "Set initial values for all unknowns in the 1d cable simulation.")
            .add_method("set_coordinate_scale_factor_3d_to_1d", &T::set_coordinate_scale_factor_3d_to_1d, "",
                "factor", "Set a factor for coordinate scaling from 3d to 1d representation.")
            .add_method("set_solver_output_verbose", &T::set_solver_output_verbose, "",
                "verbose", "Set whether the output of the 1d solver is to be verbose.")
            .add_method("set_vtk_output", &T::set_vtk_output, "",
                "file name#plot step", "Set a file name and a plotting interval for output to VTK file.")
            .add_method("set_time_steps_for_simulation_and_potential_update", &T::set_time_steps_for_simulation_and_potential_update, "",
                "simulation time step#potential update time step",
                "Set a time step size (maximum) for the 1d simulation and for the potential update.")
            // not necessary atm
            //.add_method("set_hybrid_neuron_communicator", &T::set_hybrid_neuron_communicator, "",
            //    "hybrid neuron communicator", "Set a hybrid neuron communicator.")
            .set_construct_as_smart_pointer(true);
        reg.add_class_to_group(name, "VDCC_BG_CN", tag);
    }
#endif

#ifdef NC_WITH_VM2UG
	// VDCC with Vm2UG
	{
		typedef VDCC_BG_VM2UG<TDomain> T;
		typedef VDCC_BG<TDomain> TBase;
		std::string name = std::string("VDCC_BG_VM2UG").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&,
				SmartPtr<ApproximationSpace<TDomain> >, const std::string, const char*, const std::string, const bool)>
				("function(s) as vector#subset(s) as vector#approxSpace#baseNameVmFile#timeFormat#extensionVmFile#fileInterval#fileOffset#vertexOrderOrPositionCanChange")
			.template add_constructor<void (*)(const char*, const char*,
				SmartPtr<ApproximationSpace<TDomain> >, const std::string, const char*, const std::string, const bool)>
				("function(s) as comma-separated c-string#subset(s) as comma-separated c-string#approxSpace#baseNameVmFile#timeFormat#extensionVmFile#fileInterval#fileOffset#vertexOrderOrPositionCanChange")
			.add_method("set_file_times", &T::set_file_times, "", "file interval#file offset (first file)", "set times for which files with potential values are available")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "VDCC_BG_VM2UG", tag);
	}

	// VDCC with Neuron
#ifdef NC_WITH_NEURON
	{
		typedef VDCC_BG_VM2UG_NEURON<TDomain> T;
		typedef VDCC_BG<TDomain> TBase;
		std::string name = std::string("VDCC_BG_VM2UG_NEURON").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&,
				SmartPtr<ApproximationSpace<TDomain> >, SmartPtr<membrane_potential_mapping::Transformator>, const std::string, const char*, const std::string, const bool)>
				("function(s) as vector#subset(s) as vector#approxSpace#baseNameVmFile#timeFormat#extensionVmFile#vertexOrderOrPositionCanChange")
			.template add_constructor<void (*)(const char*, const char*,
				SmartPtr<ApproximationSpace<TDomain> >, SmartPtr<membrane_potential_mapping::Transformator>, const std::string, const char*, const std::string, const bool)>
				("function(s) as comma-separated c-string#subset(s) as comma-separated c-string#approxSpace#baseNameVmFile#timeFormat#extensionVmFile#vertexOrderOrPositionCanChange")
			.add_method("set_transformator", static_cast<void (T::*) (SmartPtr<membrane_potential_mapping::Transformator>)> (&T::set_transformator), "", "", "")
			.add_method("set_provider", static_cast<void (T::*) (SmartPtr<membrane_potential_mapping::Mapper<TDomain::dim, number> >)> (&T::set_provider), "", "", "")
			.add_method("set_mapper", static_cast<void (T::*) (SmartPtr<membrane_potential_mapping::NeuronMPM>)> (&T::set_mapper), "", "", "")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "VDCC_BG_VM2UG_NEURON", tag);
	}
#endif // NC_WITH_NEURON
#endif // NC_WITH_VM2UG


	// mark for refinement functions
	{
		reg.add_function("mark_global", &mark_global<TDomain>, grp.c_str(), "", "refiner#domain", "");
		reg.add_function("mark_anisotropic", &mark_anisotropic<TDomain>, grp.c_str(), "", "refiner#domain#anisotropy threshold (<=1)", "");
		reg.add_function("mark_anisotropic_x", &mark_anisotropic_onlyX<TDomain>, grp.c_str(), "", "refiner#domain#anisotropy threshold (<=1)", "");
		reg.add_function("unmark_ranvier_areas", &unmark_ranvier_areas<TDomain>, grp.c_str(), "", "refiner#approx space#ranvier node subsets#unmark", "");
	}

	// extra commands for this plugin
	reg.add_function("compute_volume", static_cast<void (*) (ConstSmartPtr<ApproximationSpace<TDomain> >, const char*)>(&computeVolume<TDomain>), grp.c_str(),
					 "", "approxSpace#subsetNames", "outputs subset volumes");
	reg.add_function("compute_volume_of_subset", static_cast<number (*) (ConstSmartPtr<ApproximationSpace<TDomain> >, int)>(&computeVolume<TDomain>), grp.c_str(),
					 "volume of the subset", "approxSpace # subset index", "calculates subset volume");
}

/**
 * Function called for the registration of Dimension dependent parts.
 * All Functions and Classes depending on the Dimension
 * are to be placed here when registering. The method is called for all
 * available Dimension types, based on the current build options.
 *
 * @param reg		registry
 * @param grp		group for sorting of functionality
 */
template <int dim>
static void Dimension(Registry& reg, string grp)
{
	string suffix = GetDimensionSuffix<dim>();
	string tag = GetDimensionTag<dim>();

}

/**
 * Function called for the registration of Algebra dependent parts.
 * All Functions and Classes depending on Algebra
 * are to be placed here when registering. The method is called for all
 * available Algebra types, based on the current build options.
 *
 * @param reg		registry
 * @param grp		group for sorting of functionality
 */
template <typename TAlgebra>
static void Algebra(Registry& reg, string grp)
{
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();

#ifdef NC_WITH_CABLENEURON
	// HybridSynapseCurrentAssembler
	{
		typedef HybridSynapseCurrentAssembler<Domain3d, TAlgebra> T;
		typedef IDomainConstraint<Domain3d, TAlgebra> TBase;
		string name = string("HybridSynapseCurrentAssembler").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr<ApproximationSpace<Domain3d> >,
				SmartPtr<ApproximationSpace<Domain3d> >,
				SmartPtr<cable_neuron::synapse_handler::SynapseHandler<Domain3d> >,
				const std::vector<std::string>&, const std::string&)>
				("3d approximation space#1d approximation space#synapse handler#subset(s) vector#"
				 "calcium function name")
			.template add_constructor<void (*)(SmartPtr<ApproximationSpace<Domain3d> >,
				SmartPtr<ApproximationSpace<Domain3d> >,
				SmartPtr<cable_neuron::synapse_handler::SynapseHandler<Domain3d> >,
				const std::vector<std::string>&, const std::string&, const std::string&)>
				("3d approximation space#1d approximation space#synapse handler#subset(s) vector#"
				 "calcium function name#ip3 function name")
			.add_method("set_current_percentage", &T::set_current_percentage, "", "", "")
			.add_method("set_synaptic_radius", &T::set_synaptic_radius, "", "radius", "")
			.add_method("set_valency", &T::set_valency, "", "", "")
			.add_method("set_scaling_factors", &T::set_scaling_factors, "", "", "")
			.add_method("set_ip3_production_params", &T::set_ip3_production_params, "", "single synapse maximal production rate#decay rate", "")
			.add_method("set_3d_neuron_ids", &T::set_3d_neuron_ids, "", "", "")
			.set_construct_as_smart_pointer(true);

		reg.add_class_to_group(name, "HybridSynapseCurrentAssembler", tag);
	}
#endif

}

/**
 * Function called for the registration of Domain and Algebra independent parts.
 * All Functions and Classes not depending on Domain and Algebra
 * are to be placed here when registering.
 *
 * @param reg		registry
 * @param grp		group for sorting of functionality
 */
static void Common(Registry& reg, string grp)
{
	{
		typedef IMembraneTransporter T;
		std::string name = std::string("MembraneTransporter");
		reg.add_class_<T>(name, grp)
			.add_method("set_constant", &T::set_constant, "", "index i#value v",
						"sets a constant value v for the i-th passed unknown", "")
			.add_method("num_fluxes", &T::n_fluxes, "number of fluxes this transport mechanism calculates", "", "", "")
			.add_method("print_units", &T::print_units, "", "",
						"prints out the units used in the implementation of this membrane transport mechanism", "")
			.add_method("set_scale_inputs", &T::set_scale_inputs, "", "scaling factors (same number and order as for the constructor)",
						"Sets scaling factors for conversion of user's input variable units to the units of the implementation of "
						"this membrane transport mechanism.", "Default values: 1.0 (no scaling).")
			.add_method("set_scale_input", &T::set_scale_input, "", "index#scaling factor",
						"Sets a scaling factor for conversion of the user's input variable (specified by first parameter) units to the "
						"units of the implementation of this membrane transport mechanism.", "")
			.add_method("set_scale_fluxes", &T::set_scale_fluxes, "", "scaling factors",
						"Sets scaling factors for conversion of calculated fluxes to the units employed by the user.",
						"Default values: 1.0 (no scaling).")
			.add_method("set_scale_flux", &T::set_scale_flux, "", "index#scaling factor",
						"Sets a scaling factor for conversion of the calculated flux (specified by first parameter) to the unit employed "
						"by the user.", "");
			//.add_method("calc_flux", static_cast<number (T::*) (const std::vector<number>&, size_t) const>(&T::calc_flux), "", "input values#flux index#output flux",
			//		"calculates the specified flux through this mechanism", "");
			/* does not work, since vectors have to be const for exchange with lua
			.add_method("calc_flux", &T::calc_flux, "", "input values#output flux(es)",
						"calculates the flux(es) through this mechanism", "")
			.add_method("calc_flux_deriv", &T::calc_flux_deriv, "", "input values#output flux derivatives",
						"calculates the flux derivatives through this mechanism", "");
			*/
	}
	{
		typedef OhmicLeakage T;
		typedef IMembraneTransporter TBase;
		std::string name = std::string("OhmicLeakage");
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor<void (*)(const char*)>
				("Functions as comma-separated string with the following order: "
				 "\"inner charge density\", \"outer charge density\", \"inner potential\", \"outer potential\"")
			.add_constructor<void (*)(const std::vector<std::string>&)>
				("Function vector with the following order: "
				 "\"inner charge density\", \"outer charge density\", \"inner potential\", \"outer potential\"")
			.add_method("set_conductance", &T::set_conductance, "", "g_L", "")
			.add_method("set_reversal_potential", &T::set_reversal_potential, "", "E_L", "")
			.set_construct_as_smart_pointer(true);
	}
	{
		typedef OhmicLeakageCharges T;
		typedef IMembraneTransporter TBase;
		std::string name = std::string("OhmicLeakageCharges");
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor<void (*)(const char*)>
				("Functions as comma-separated string with the following order: "
				 "\"inner charge density\", \"outer charge density\", \"inner potential\", \"outer potential\"")
			.add_constructor<void (*)(const std::vector<std::string>&)>
				("Function vector with the following order: "
				 "\"inner charge density\", \"outer charge density\", \"inner potential\", \"outer potential\"")
			.add_method("set_conductance", &T::set_conductance, "", "g_L", "")
			.add_method("set_reversal_potential", &T::set_reversal_potential, "", "E_L", "")
			.set_construct_as_smart_pointer(true);
	}
	{
		typedef IP3R T;
		typedef IMembraneTransporter TBase;
		std::string name = std::string("IP3R");
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor<void (*)(const char*)>
				("Functions as comma-separated string with the following order: "
				 "\"cytosolic calcium, endoplasmic calcium, ip3\"")
			.add_constructor<void (*)(const std::vector<std::string>&)>
				("Function vector with the following order: "
				 "{\"cytosolic calcium\", \"endoplasmic calcium\", \"ip3\"}")
			.set_construct_as_smart_pointer(true);
	}
	{
		typedef RyR T;
		typedef IMembraneTransporter TBase;
		std::string name = std::string("RyR");
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor<void (*)(const char* )>
				("Functions as comma-separated string with the following order: "
				 "{\"cytosolic calcium\", \"endoplasmic calcium\"}")
			.add_constructor<void (*)(const std::vector<std::string>&)>
				("Function vector with the following order: "
				 "{\"cytosolic calcium\", \"endoplasmic calcium\"}")
			.set_construct_as_smart_pointer(true);
	}
	{
		typedef SERCA T;
		typedef IMembraneTransporter TBase;
		std::string name = std::string("SERCA");
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor<void (*)(const char*)>
				("Functions as comma-separated string with the following order: "
				 "{\"cytosolic calcium\", \"endoplasmic calcium\"}")
			.add_constructor<void (*)(const std::vector<std::string>&)>
				("Function vector with the following order: "
				 "{\"cytosolic calcium\", \"endoplasmic calcium\"}")
			.set_construct_as_smart_pointer(true);
	}
	{
		typedef Leak T;
		typedef IMembraneTransporter TBase;
		std::string name = std::string("Leak");
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor<void (*)(const char*)>
				("Functions as comma-separated string with the following order: "
				 "{\"source\", \"target\"}")
			.add_constructor<void (*)(const std::vector<std::string>&)>
				("Function vector with the following order: "
				 "{\"source\", \"target\"}")
			.set_construct_as_smart_pointer(true);
	}
	{
		typedef PMCA T;
		typedef IMembraneTransporter TBase;
		std::string name = std::string("PMCA");
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor<void (*)(const char*)>
				("Functions as comma-separated string with the following order: "
				 "{\"cytosolic calcium\", \"extracellular calcium\"}")
			.add_constructor<void (*)(const std::vector<std::string>&)>
				("Function vector with the following order: "
				 "{\"cytosolic calcium\", \"extracellular calcium\"}")
			.set_construct_as_smart_pointer(true);
	}
	{
		typedef NCX T;
		typedef IMembraneTransporter TBase;
		std::string name = std::string("NCX");
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor<void (*)(const char*)>
				("Functions as comma-separated string with the following order: "
				 "{\"cytosolic calcium\", \"extracellular calcium\"}")
			.add_constructor<void (*)(const std::vector<std::string>&)>
				("Function vector with the following order: "
				 "{\"cytosolic calcium\", \"extracellular calcium\"}")
			.set_construct_as_smart_pointer(true);
	}
	{
		typedef MCU T;
		typedef IMembraneTransporter TBase;
		std::string name = std::string("MCU");
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor<void (*)(const char*)>
				("Functions as comma-separated string with the following order: "
				 "{\"cytosolic calcium\", \"mitochondrial calcium\"}")
			.add_constructor<void (*)(const std::vector<std::string>&)>
				("Function vector with the following order: "
				 "{\"cytosolic calcium\", \"mitochondrial calcium\"}")
			.add_method("set_mit_volume", &T::set_mit_volume,
						"Sets mitochondrial volume.")
			.add_method("set_mit_surface", &T::set_mit_surface,
						"Sets mitochondrial surface.")
			.add_method("set_pi_cyt", &T::set_pi_cyt,
						"Sets cytosolic phosphate concentration.")
			.add_method("set_psi", &T::set_psi,
						"Sets mitochondrial membrane potential.")
			.add_method("set_mg_cyt", &T::set_mg_cyt,
						"Sets cytosolic Mg2+ concentration.")
			.add_method("set_mg_mit", &T::set_mg_mit,
						"Sets mitochondrial Mg2+ concentration.")
			.add_method("set_rate_constant", &T::set_rate_constant,
						"Sets rate constant k.")
			.add_method("get_flux", &T::get_flux, "Debug method.")
			.set_construct_as_smart_pointer(true);
	}
	{
		typedef MNCX T;
		typedef IMembraneTransporter TBase;
		std::string name = std::string("MNCX");
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor<void (*)(const char*)>
				("Functions as comma-separated string with the following order: "
				 "{\"cytosolic calcium\", \"mitochondrial calcium\", \"cytosolic sodium\", \"mitochondrial sodium\"}")
			.add_constructor<void (*)(const std::vector<std::string>&)>
				("Function vector with the following order: "
				 "{\"cytosolic calcium\", \"mitochondrial calcium\", \"cytosolic sodium\", \"mitochondrial sodium\"}")
			.add_method("set_mit_volume", &T::set_mit_volume,
						"Sets mitochondrial volume.")
			.add_method("set_mit_surface", &T::set_mit_surface,
						"Sets mitochondrial surface.")
			.add_method("set_psi", &T::set_psi,
						"Sets mitochondrial membrane potential.")
			.add_method("get_flux", &T::get_flux, "Debug method.")
			.set_construct_as_smart_pointer(true);
	}
	{
		typedef NMDAR T;
		typedef IMembraneTransporter TBase;
		std::string name = std::string("NMDAR");
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor<void (*)(const char*)>
				("Functions as comma-separated string with the following order: "
				"extracellular calcium, intracellular calcium")
			.add_constructor<void (*)(const std::vector<std::string>&)>
				("Function vector with the following order: {extracellular calcium, intracellular calcium}")
			.add_method("set_activation_time", &T::set_activation_time, "", "", "")
			.add_method("set_decay_time", &T::set_decay_time, "", "", "")
			.add_method("set_permeability", &T::set_permeability, "", "", "")
			.add_method("set_membrane_potential", &T::set_membrane_potential, "", "", "")
			.add_method("set_temperature", &T::set_temperature, "", "", "")
			.set_construct_as_smart_pointer(true);
	}
	{
		typedef ActionPotentialTrain T;
		std::string name = std::string("ActionPotentialTrain");
		reg.add_class_<T>(name, grp)
			.add_constructor()
			.add_constructor<void (*)(number, number, number, number)>
				("stimBegin#stimEnd#stimFreq#basicVoltage")
			.add_method("membrane_potential", &T::membrane_potential,
						"Returns membrane potential to given frequency stimulation interval.")
			.set_construct_as_smart_pointer(true);
	}

	// build bouton
	{
        reg.add_function("BuildBouton", &BuildBouton, grp,
                         "", "bExtSpace#radius#numRefinements#numReleaseSites#TbarHeight#TbarLegRadius#TbarTopRadius#TbarTopHeight#fileName",
                         "Generates a drosophila NMJ bouton volume grid.");
	}

	// build spine
	{
		// TODO: Rename "BuildSpine", remove ineffective parameters
        reg.add_function("BuildDendrite", &BuildSpine, grp,
                         "", "geometric param vector (cytosol radius, ER radius, dendrite length, spine position, "
                         "spine ER neck radius, spine ER neck length, spine ER head radius, spine ER head length, "
                         "spine neck radius, spine neck length, spine head radius, spine head length)"
                         "options vector (build a synapse? [ineffective], build ER?, build spine ER?, "
                         "synapse at different location? [ineffective], build spine ER head?)"
                         "#fileName",
                         "Generates a dendritic spine with a portion of the connected dendrite.");
	}

#ifndef UG_FOR_VRL
	// test neurite projector
	{
		reg.add_function("test_smoothing", &test_smoothing, "", "", "");
		reg.add_function("test_neurite_projector", &test_neurite_projector_with_four_section_tube, "", "", "");
		reg.add_function("test_neurite_projector_with_bp", &test_neurite_projector_with_four_section_tube_and_branch_point, "", "", "");
		reg.add_function("test_import_swc", &test_import_swc, "", "file name # anisotropy # refinements", "");
		reg.add_function("test_import_swc_with_er", &test_import_swc_with_er, "",
			"swc file name (input) # ugx file name (output) # ER scale factor # anisotropy # refinements", "");
		reg.add_function("test_import_swc_surf", &test_import_swc_surf, "", "file name#correction", "");
		reg.add_function("test_import_swc_1d", &test_import_swc_1d, "", "file name # anisotropy # refinements", "");
		//reg.add_function("apply_neurite_projector", &apply_neurite_projector, "", "multigrid, neurite projector", "");
	}
#endif

#ifdef UG_DIM_3
	{
		typedef NeuriteAxialRefinementMarker T;
		std::string name = std::string("NeuriteAxialRefinementMarker");
#ifdef NC_WITH_PARMETIS
		typedef parmetis::IUnificator<Volume> TBase;
		reg.add_class_<T, TBase>(name, grp)
#else
		reg.add_class_<T>(name, grp)
#endif
			.add_constructor<void (*)(SmartPtr<Domain3d>)>("domain")
			.add_method("mark", &T::mark, "", "refiner", "Marks neurites for axial refinement.")
			.set_construct_as_smart_pointer(true);
	}

	// mark for refinement functions
	{
		reg.add_function("MarkNeuriteForAxialRefinement", &MarkNeuriteForAxialRefinement, grp.c_str(), "", "refiner#domain", "");
	}
#endif
}

}; // end Functionality


struct NonBlockedFunctionality
{
	template <typename TDomain, typename TAlgebra>
	static void DomainAlgebra(Registry& reg, string grp)
	{
		string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
		string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

		typedef GridFunction<TDomain, TAlgebra> TGridFunction;

		// discrete RyR
		{
			typedef RyRDiscrete<TDomain, TAlgebra> T;
			typedef IDomainConstraint<TDomain, TAlgebra> TBase;
			string name = string("RyRDiscrete").append(suffix);
			reg.add_class_<T, TBase>(name, grp)
				.template add_constructor<void (*)(const char*, const char*)>
					("function names (comma-separated c-string) # subset names (comma-separated c-string)")
				.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>
					("function names (vector of strings) # subset names (vector of strings)")
				.add_method("calculate_steady_state", &T::calculate_steady_state, "", "solution", "")
				.add_method("set_cutoff_open_probability", &T::set_cutoff_open_probability, "", "cutoff probability",
					"set open probability below which the channel is supposed to be certainly closed")
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "RyRDiscrete", tag);
		}

		reg.add_function("max_ryr_flux_density", &maxRyRFluxDensity<TGridFunction, RyRDiscrete<TDomain, TAlgebra> >, grp.c_str(),
			"", "solution # function names for ca_cyt, ca_er, c1, c2 as c-style string #"
			"RyR-carrying membrane subset names as c-style string # RyR channel",
			"maximal flux density through RyR channel (mol/(m^2*s))");

	}
};


// end group plugin_neuro_collection
/// \}

} // end namespace neuro_collection


/**
 * This function is called when the plugin is loaded.
 */
extern "C" void
InitUGPlugin_neuro_collection(Registry* reg, string grp)
{
	grp.append("/neuro_collection");
	typedef neuro_collection::Functionality Functionality;

	GlobalAttachments::declare_attachment<ANumber>("diameter");
	typedef Attachment<NeuriteProjector::SurfaceParams> NPSurfParam;
	GlobalAttachments::declare_attachment<NPSurfParam>("npSurfParams", true);

	try
	{
		RegisterCommon<Functionality>(*reg,grp);
		//RegisterDimensionDependent<Functionality>(*reg,grp);
		RegisterDomainDependent<Functionality>(*reg,grp);
		RegisterAlgebraDependent<Functionality>(*reg,grp);
		RegisterDomainAlgebraDependent<Functionality>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);

	// some algebra-dependent code is only meant for non-blocked algebras
	typedef boost::mpl::list
	<
		#ifdef UG_CPU_1
		CPUAlgebra,
		#endif
		end_boost_list
	> CompileNonBlockedAlgebraList;

	typedef neuro_collection::NonBlockedFunctionality NonBlockedFunctionality;
	try
	{
		RegisterDomainAlgebraDependent<NonBlockedFunctionality, CompileDomainList, CompileNonBlockedAlgebraList>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

} // namespace ug
