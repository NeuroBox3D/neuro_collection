/**
 * NEURO collection
 *
 * This plugin aims to collect and unify any neuro-related functionality that is
 * developed for ug4 and that the author thinks of as usable for others.
 *
 *
 *  Created on: 13.06.2014
 *      Author: mbreit
**/

#include "bridge/util.h"

// replace this with util_domain_dependent.h or util_algebra_dependent.h
// to speed up compilation time
#include "bridge/util_domain_algebra_dependent.h"

#include "buffer_fv1.h"
#include "one_sided_borg_graham_fv1.h"
#include "dependent_neumann_boundary_fv1.h"
#include "one_sided_membrane_transport_fv1.h"
#include "two_sided_membrane_transport_fv1.h"


using namespace std;
using namespace ug::bridge;

namespace ug{
namespace neuro_collection{

/** 
 *  \defgroup plugin_neuro_collection Plugin neuro_collection
 *  \ingroup plugins_experimental
 *  This is a plugin for neuro-related functionality.
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
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain, typename TAlgebra>
static void DomainAlgebra(Registry& reg, string grp)
{
	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

}

/**
 * Function called for the registration of Domain dependent parts.
 * All Functions and Classes depending on the Domain
 * are to be placed here when registering. The method is called for all
 * available Domain types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain>
static void Domain(Registry& reg, string grp)
{
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();

	const int dim = TDomain::dim;

	// implementation of buffering reaction disc
	{
		typedef BufferFV1<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("BufferFV1").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(const char*)>("Subset(s)")
			.add_method("add_reaction", static_cast<void (T::*)
				 (const char*, const char*, number, number, number)>(&T::add_reaction), "",
				  "buffering substance | selection | value=[\"clb\"] # "
				  "buffered substance | selection | value=[\"ca_cyt\"] # "
				  "total buffer | default |	value=160.0e-6D # "
				  "association rate | default | value=27.0e06D # "
				  "dissociation rate | default | value=19.0D",
				  "add a new reaction definition")
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

	// implementation of a Neumann boundary depending on the unknowns
	{
		typedef DependentNeumannBoundaryFV1<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("DependentNeumannBoundaryFV1").append(suffix);
		reg.add_class_<T, TBase >(name, grp);
		reg.add_class_to_group(name, "DependentNeumannBoundaryFV1", tag);
	}


	// one-sided membrane transport systems
	{
		typedef OneSidedMembraneTransportFV1<TDomain> T0;
		typedef DependentNeumannBoundaryFV1<TDomain> TBase;
		string name = string("OneSidedMembraneTransportFV1").append(suffix);
		reg.add_class_<T0, TBase >(name, grp)
			.add_method("set_density_function", static_cast<void (T0::*) (const char*)> (&T0::set_density_function),
							"", "", "add a density function")
			.add_method("set_density_function", static_cast<void (T0::*) (SmartPtr<CplUserData<number,dim> >)>
					(&T0::set_density_function), "", "", "add a density function");
		reg.add_class_to_group(name, "OneSidedMembraneTransportFV1", tag);
	}

	// one-sided plasma membrane transport
	{
		typedef OneSidedMembraneTransportFV1<TDomain> T0;
		typedef OneSidedPMCAFV1<TDomain> T1;
		string name = string("OneSidedPMCAFV1").append(suffix);
		reg.add_class_<T1, T0>(name, grp)
			.template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "OneSidedPMCAFV1", tag);

		typedef OneSidedNCXFV1<TDomain> T2;
		name = string("OneSidedNCXFV1").append(suffix);
		reg.add_class_<T2, T0>(name, grp)
			.template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "OneSidedNCXFV1", tag);

		typedef OneSidedPMCalciumLeakFV1<TDomain> T3;
		name = string("OneSidedPMCalciumLeakFV1").append(suffix);
		reg.add_class_<T3, T0>(name, grp)
			.template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "OneSidedPMCalciumLeakFV1", tag);
	}

	// one-sided Borg-Graham implementation
	{
		typedef OneSidedBorgGrahamFV1<TDomain> TBG0;
		typedef OneSidedMembraneTransportFV1<TDomain> TBase;
		string name = string("OneSidedBorgGrahamFV1").append(suffix);
		reg.add_class_<TBG0, TBase >(name, grp)
			.add_method("set_channel_type_N", &TBG0::template set_channel_type<TBG0::BG_Ntype>,
						"", "", "set the channel type to N")
			.add_method("set_channel_type_L", &TBG0::template set_channel_type<TBG0::BG_Ltype>,
						"", "", "set the channel type to L")
			.add_method("set_channel_type_T", &TBG0::template set_channel_type<TBG0::BG_Ttype>,
						"", "", "set the channel type to T")
			.add_method("init", &TBG0::init, "", "time", "initialize the Borg-Graham object")
			.add_method("update_gating", &TBG0::update_gating, "", "time",
						"update gating \"particles\" to new time");
		reg.add_class_to_group(name, "OneSidedBorgGrahamFV1", tag);

		typedef OneSidedBorgGrahamFV1WithVM2UG<TDomain> TBG1;
		name = string("OneSidedBorgGrahamFV1WithVM2UG").append(suffix);
		reg.add_class_<TBG1, TBG0 >(name, grp)
			.template add_constructor<void (*)(const char*, const char*, ApproximationSpace<TDomain>&,
											   const std::string, const char*, const std::string, const bool)>
				("function(s)#subset(s)#approxSpace#baseNameVmFile#timeFormat#extensionVmFile#vertexOrderOrPositionCanChange")
			.add_method("update_potential", &TBG1::update_potential, "", "time", "update membrane potential to new time")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "OneSidedBorgGrahamFV1WithVM2UG", tag);
	}


	// implementation of two-sided membrane transport systems
	{
		typedef TwoSidedMembraneTransportFV1<TDomain> T;
		typedef FV1InnerBoundaryElemDisc<TDomain> TBase;
		string name = string("TwoSidedMembraneTransportFV1").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.add_method("set_density_function", static_cast<void (T::*) (const char*)> (&T::set_density_function),
							"", "", "add a density function")
			.add_method("set_density_function", static_cast<void (T::*) (SmartPtr<CplUserData<number,dim> >)>
					(&T::set_density_function), "", "", "add a density function");
		reg.add_class_to_group(name, "TwoSidedMembraneTransportFV1", tag);
	}

	// implementation of two-sided ER membrane calcium transport systems
	{
		typedef TwoSidedERCalciumTransportFV1<TDomain> TE0;
		typedef TwoSidedMembraneTransportFV1<TDomain> TEBase;
		string name = string("TwoSidedERCalciumTransportFV1").append(suffix);
		reg.add_class_<TE0, TEBase >(name, grp);
		reg.add_class_to_group(name, "TwoSidedERCalciumTransportFV1", tag);

		typedef TwoSidedIP3RFV1<TDomain> TE1;
		name = string("TwoSidedIP3RFV1").append(suffix);
		reg.add_class_<TE1, TE0>(name, grp)
			.template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "TwoSidedIP3RFV1", tag);

		typedef TwoSidedRyRFV1<TDomain> TE2;
		name = string("TwoSidedRyRFV1").append(suffix);
		reg.add_class_<TE2, TE0>(name, grp)
			.template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "TwoSidedRyRFV1", tag);

		typedef TwoSidedSERCAFV1<TDomain> TE3;
		name = string("TwoSidedSERCAFV1").append(suffix);
		reg.add_class_<TE3, TE0>(name, grp)
			.template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "TwoSidedSERCAFV1", tag);

		typedef TwoSidedERCalciumLeakFV1<TDomain> TE4;
		name = string("TwoSidedERCalciumLeakFV1").append(suffix);
		reg.add_class_<TE4, TE0>(name, grp)
			.template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "TwoSidedERCalciumLeakFV1", tag);
	}
}

/**
 * Function called for the registration of Dimension dependent parts.
 * All Functions and Classes depending on the Dimension
 * are to be placed here when registering. The method is called for all
 * available Dimension types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
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
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TAlgebra>
static void Algebra(Registry& reg, string grp)
{
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();

}

/**
 * Function called for the registration of Domain and Algebra independent parts.
 * All Functions and Classes not depending on Domain and Algebra
 * are to be placed here when registering.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
static void Common(Registry& reg, string grp)
{

}

}; // end Functionality

// end group plugin_neuro_collection
/// \}

} // end namespace neuro_collection


/**
 * This function is called when the plugin is loaded.
 */
extern "C" void
InitUGPlugin_neuro_collection(Registry* reg, string grp)
{
	grp.append("neuro_collection");
	typedef neuro_collection::Functionality Functionality;

	try{
		RegisterCommon<Functionality>(*reg,grp);
		RegisterDimensionDependent<Functionality>(*reg,grp);
		RegisterDomainDependent<Functionality>(*reg,grp);
		RegisterAlgebraDependent<Functionality>(*reg,grp);
		RegisterDomainAlgebraDependent<Functionality>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// namespace ug
