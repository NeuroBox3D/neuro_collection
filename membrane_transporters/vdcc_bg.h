/*
 * vdcc_bg.h
 *
 *  Created on: 05.02.2013
 *      Author: mbreit
 */

#ifndef __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__VDCC_BG_H__
#define __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__VDCC_BG_H__

#include "common/common.h"
#include "membrane_transporter_interface.h"
#include "lib_disc/domain.h"
#include "lib_grid/lg_base.h"
#include "lib_disc/spatial_disc/elem_disc/inner_boundary/inner_boundary.h"
#include "../plugins/experimental/membrane_potential_mapping/vm2ug.h"

#include <locale>	// for control over the decimal separator (point instead of comma, please!)


namespace ug{
namespace neuro_collection{

///@addtogroup plugin_neuro_collection
///@{

/// Class interface for Borg Graham type VGCCs of the plasma membrane.
/** This class is an interface for the Borg-Graham-type voltage-gated calcium channels
 *	in the plasma membrane (see chapter 8 of "Interpretations of data and mechanisms
 *	for hippocampal pyramidal cell models", Borg-Graham (1998) in Cerebral Cortex,
 *	Vol. 13: Models of Cortical Circuits.
 *
 *	The unknowns of the discretization (the so-called "gates" of the channel) are
 *	not added to the system unknowns, but held separately as vertex attachments, since
 *	they are governed by ODEs and are therefore merely updated by an implicit time
 *	schema once per time step. This will have to be done by the user (in the lua script)
 *	via the ITimeDiscretization method prepare_step_elem() before the execution of any
 *	time step assembling.
 *
 *	The class implements simple versions of N-, L- and T type channels but may very
 *	well be generalized to more complex models thereof. The type of channel can be set
 *	before calling the init() method, N-type channel parameters are taken as default.
 *
 *	This class does not handle the procuration of the values for the membrane
 *	potential. This must be dealt with in a specialization of this class.
 *
 *	Any class specializing this interface _must_ implement the virtual method:
 *	- void update_potential(number newTime);
 *	they _can_ re-implement the virtual method
 *	- void init(number time).
 *
 *
 *	The units required for this discretization are:
 * 		V_m	: mV		membrane voltage (VDCC_BG_VM2UG & VDCC_BG_VM2UG_NEURON)
 * 		V_m	: V			membrane voltage (VDCC_BG & VDCC_BG_UserData)
 * 		t	: s			time
 *  	f	: mol*s^-1	ionic flux
 *
 *  Remarks:
 *  	- Internally all membrane potentials are attached to elements in [V] for
 *  	  the element discretization!
 *
 *  	- VDCC_BG & VDCC_BG_UserData use [ms] and [mV] in gating and flux calculations!
 *  	  This is due to the use of the original gating parameter sets by Borg-Graham
 *  	  with tau_0 in [ms] and V_12 in [mV].
 *  	  Note, that the update potential method already takes care of this, when calling
 *  	  the corresponding gating & flux calculation methods.
 *
**/

template<typename TDomain>
class VDCC_BG : public IMembraneTransporter
{
	public:
		/// channel types N, L and T
		enum {BG_Ntype, BG_Ltype, BG_Ttype};
		enum{_CCYT_=0, _CEXT_};

		static const int dim =TDomain::dim;	//!< world dimension

		// some type definitions
		typedef VDCC_BG<TDomain> this_type;
		typedef typename GeomObjBaseTypeByDim<dim>::base_obj_type elem_t;
		typedef typename elem_t::side side_t;

	protected:
		const number R;			///< universal gas constant
		const number T;			///< temperature
		const number F;			///< Faraday constant

		/// holds the parameters of a channel type
		struct GatingParams
		{
			GatingParams(number _z, number _v, number _t) : z(_z), V_12(_v), tau_0(_t){};
			number z;
			number V_12;
			number tau_0;
		};

	public:
	// inherited from IMembraneTransport
		/**
		 * @brief Constructor for the Borg-Graham type channel interface
		 *
		 * This constructor not only needs information on the functions involved, but also on the subsets
		 * and the approximation space, as it needs to create side attachments for each of the sides in
		 * the subsets involved.
		 *
		 * @param fcts		functions vector
		 * @param subsets	subsets vector
		 * @param approx	underlying approximation space
		 */
		VDCC_BG
		(
			const std::vector<std::string>& fcts,
			const std::vector<std::string>& subsets,
			SmartPtr<ApproximationSpace<TDomain> > approx
		);

		/**
		 * @brief Constructor for the Borg-Graham type channel interface
		 *
		 * This constructor not only needs information on the functions involved, but also on the subsets
		 * and the approximation space, as it needs to create side attachments for each of the sides in
		 * the subsets involved.
		 *
		 * @param fcts		functions as comma-separated c-string
		 * @param subsets	subsets as comma-separated c-string
		 * @param approx	underlying approximation space
		 */
		VDCC_BG
		(
			const char* fcts,
			const char* subsets,
			SmartPtr<ApproximationSpace<TDomain> > approx
		);

		/// @copydoc IMembraneTransporter::IMembraneTransporter()
		virtual ~VDCC_BG();

		/// @copydoc IMembraneTransporter::prep_timestep_elem()
		virtual void prep_timestep_elem
		(
			const number time,
			const LocalVector& u,
			GridObject* elem
		);

		/// @copydoc IMembraneTransporter::calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const
		virtual void calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const;

/*
		/// @copydoc IMembraneTransporter::calc_flux(const std::vector<number>&, number&,  flux) const
		virtual number calc_flux(const std::vector<number>& u, size_t index);
*/
		/// @copydoc IMembraneTransporter::calc_flux_deriv()
		virtual void calc_flux_deriv(const std::vector<number>& u, GridObject* e, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const;

		/// @copydoc IMembraneTransporter::n_dependencies()
		virtual const size_t n_dependencies() const;

		/// @copydoc IMembraneTransporter::n_fluxes()
		virtual size_t n_fluxes() const;

		/// @copydoc IMembraneTransporter::flux_from_to()
		virtual const std::pair<size_t,size_t> flux_from_to(size_t flux_i) const;

		/// @copydoc IMembraneTransporter::name()
		virtual const std::string name() const;

		/// @copydoc IMembraneTransporter::check_supplied_functions()
		virtual void check_supplied_functions() const;

		/// @copydoc IMembraneTransporter::print_units()
		virtual void print_units() const;

	// own methods
		/// sets the channel type
		template<int TType> void set_channel_type();

		/// initializes the defined channel type
		/** During the initialization, the necessary attachments are attached to the vertices
		 *	and their values calculated by the equilibrium state for the start membrane potential.
		**/

		/**
		 * @brief Sets the permeability of this channel
		 * @param perm    permeability values
		 */
		void set_permeability(const number perm);

		virtual void init(number time);

		/// updates the potential values in the corresponding attachments to new time.
		/**
		 * This method needs to be called before update_gating() if potential is non-constant.
		 * @param newTime new point in time
		 */
		virtual void update_potential(side_t* elem) = 0;

		/// updates the gating parameters
		/**
		 * This method needs to be called before ionic_current().
		 * @param newTime new point in time
		 */
		void update_gating(side_t* elem);

		/// updates internal time if necessary
		virtual void update_time(number newTime);

	protected:
		/// calculates the equilibrium state of a gating "particle"
		/** The calculation is done with respect to the given gating parameters set (which represents
		 *	one gating "particle") and the given membrane potential (to be specified in [mV]!).
		**/
		number calc_gating_start(GatingParams& gp, number Vm);

		/// calculates the next state of a gating "particle"
		/** The calculation is done with respect to the given gating parameters set (which represents
		 *	one gating "particle"), the current value of this gating particle as well as the given membrane
		 *	potential (to be specified in [mV]!). The new value represents the "particle" state
		 *	at current time + dt (to be specified in [ms]!).
		**/
		void calc_gating_step(GatingParams& gp, number Vm, number dt, number& currVal);

	private:
		void after_construction();

	protected:
		/// whether this channel disposes of an inactivating gate
		bool has_hGate() const {return this->m_channelType == BG_Ntype || this->m_channelType == BG_Ttype;}

	protected:
		SmartPtr<TDomain> m_dom;							//!< underlying domain
		SmartPtr<Grid> m_mg;								//!< underlying multigrid
		SmartPtr<DoFDistribution> m_dd;						//!< underlying surface dof distribution
		ConstSmartPtr<MGSubsetHandler> m_sh;				//!< underlying subset handler
		typename TDomain::position_accessor_type& m_aaPos;	//!< underlying position accessor

		std::vector<std::string> m_vSubset;					//!< subsets this channel exists on

		ADouble m_MGate;							//!< activating gating "particle"
		ADouble m_HGate;							//!< inactivating gating "particle"
		ADouble m_Vm;								//!< membrane voltage (in Volt)

		Grid::AttachmentAccessor<side_t, ADouble> m_aaMGate;	//!< accessor for activating gate
		Grid::AttachmentAccessor<side_t, ADouble> m_aaHGate;	//!< accessor for inactivating gate
		Grid::AttachmentAccessor<side_t, ADouble> m_aaVm;		//!< accessor for membrane potential

		GatingParams m_gpMGate;						//!< gating parameter set for activating gate
		GatingParams m_gpHGate;						//!< gating parameter set for inactivating gate

		number m_time;								//!< current time
		number m_oldTime;							//!< time step before current time

		number m_perm;								//!< channel permeability [m^3/s] (= diff coeff * cross section / membrane thickness)
		int m_mp, m_hp;								//!< powers for gating parameters

		int m_channelType;							//!< channel type

		bool m_initiated;							//!< indicates whether channel has been initialized by init()
};


/// Borg Graham type VGCCs with Vm2uG membrane potential supply.
/** This class is a specialization of the Borg-Graham interface.
 *	It supplies the channel with the necessary membrane potential values by a Vm2uG object,
 *	that has to be fed file names for the files containing the V_m values produced by Neuron.
**/

template<typename TDomain>
class VDCC_BG_VM2UG : public VDCC_BG<TDomain>
{
	protected:
		using VDCC_BG<TDomain>::R;			//!< universal gas constant
		using VDCC_BG<TDomain>::T;			//!< temperature (310K)
		using VDCC_BG<TDomain>::F;			//!< Faraday constant
		using VDCC_BG<TDomain>::has_hGate;
		typedef typename GeomObjBaseTypeByDim<TDomain::dim>::base_obj_type elem_t;
		typedef typename elem_t::side side_t;

	public:
		typedef Vm2uG<std::string> vmProvType;

	public:
		/**
		 * @brief constructor with vectors
		 *
		 * @param fcts			functions as vector of string
		 * @param subsets		subsets as vector of string
		 * @param approx		approximation space
		 * @param baseName		base file name for V_m file input (default: "timesteps/timestep_")
		 * @param timeFmt		format of time in file names (default: "%.4f")
		 * @param ext			file extension of input files (default: ".dat")
		 * @param posCanChange	whether coordinates in input files can change (default: false)
		 */
		VDCC_BG_VM2UG
		(
			const std::vector<std::string>& fcts,
			const std::vector<std::string>& subsets,
			SmartPtr<ApproximationSpace<TDomain> > approx,
			const std::string baseName = "timesteps/timestep_",
			const char* timeFmt = "%.4f",
			const std::string ext = ".dat",
			const bool posCanChange = false
		);

		/**
		 * @brief constructor with c-strings
		 *
		 * @param fcts			functions as comma-separated c-string
		 * @param subsets		subsets as comma-separated c-string
		 * @param approx		approximation space
		 * @param baseName		base file name for V_m file input (default: "timesteps/timestep_")
		 * @param timeFmt		format of time in file names (default: "%.4f")
		 * @param ext			file extension of input files (default: ".dat")
		 * @param posCanChange	whether coordinates in input files can change (default: false)
		 */
		VDCC_BG_VM2UG
		(
			const char* fcts,
			const char* subsets,
			SmartPtr<ApproximationSpace<TDomain> > approx,
			const std::string baseName = "timesteps/timestep_",
			const char* timeFmt = "%.4f",
			const std::string ext = ".dat",
			const bool posCanChange = false
		);

		/// destructor
		virtual ~VDCC_BG_VM2UG();

		/// @copydoc VDCC_BG<TDomain>::init()
		virtual void init(number time);

		/// @copydoc VDCC_BG<TDomain>::update_potential()
		virtual void update_potential(side_t* elem);

		/// @copydoc VDCC_BG<TDomain>::update_time()
		virtual void update_time(number newTime);

		/// @copydoc IMembraneTransporter::print_units()
		virtual void print_units() const;

		/**
		 * @brief setting the times for which files are present
		 * Using this method, it is possible to tell the VDCC implementation at which
		 * points in time there exists a file with the corresponding membrane potential
		 * values. Points in time must be separated by a constant interval and may be
		 * shifted by an offset.
		 *
		 * @param fileInterval	constant interval separating two points in time
		 * @param fileOffset	offset (the first point in time, default: 0.0)
		 */
		void set_file_times(const number fileInterval, const number fileOffset = 0.0)
		{
			m_fileInterval = fileInterval;
			m_fileOffset = fileOffset;
		}

	private:
		vmProvType m_vmProvider;		//!< the Vm2uG object
		std::string m_tFmt;				//!< time format for the membrane potential files
		number m_fileInterval;			//!< intervals in which voltage files are available
		number m_fileOffset;			//!< offset of time intervals for which voltage files are available

		std::string m_timeAsString;
};

#ifdef MPMNEURON
/// Borg Graham type VGCCs with Vm2uG membrane potential supply by NEURON.
/** This class is a specialization of the Borg-Graham interface.
 *	It supplies the channel with the necessary membrane potential values by a Vm2uG object,
 *	that has to be fed file names for the files containing the V_m values produced by Neuron.
**/
template<typename TDomain>
class VDCC_BG_VM2UG_NEURON : public VDCC_BG<TDomain>
{
private:
	SmartPtr<Transformator> m_NrnInterpreter;
	SmartPtr<Vm2uG<std::string> > m_vmProvider;

	protected:
		using VDCC_BG<TDomain>::R;		//!< universal gas constant
		using VDCC_BG<TDomain>::T;		//!< temperature (310K)
		using VDCC_BG<TDomain>::F;		//!< Faraday constant

		typedef typename GeomObjBaseTypeByDim<TDomain::dim>::base_obj_type elem_t;
		typedef typename elem_t::side side_t;

	public:
		typedef Vm2uG<std::string> vmProvType;

	public:
		/**
		 * @brief constructor with vectors
		 *
		 * @param fcts			functions as vector of string
		 * @param subsets		subsets as vector of string
		 * @param approx		approximation space
		 * @param transformator	transformator object (default: Transformator)
		 * @param baseName		base file name for V_m file input (default: "timesteps/timestep_")
		 * @param timeFmt		format of time in file names (default: "%.4f")
		 * @param ext			file extension of input files (default: ".dat")
		 * @param posCanChange	whether coordinates in input files can change (default: false)
		 */
		VDCC_BG_VM2UG_NEURON
		(
			const std::vector<std::string>& fcts,
			const std::vector<std::string>& subsets,
			SmartPtr<ApproximationSpace<TDomain> > approx,
			SmartPtr<Transformator> transformator = make_sp(new Transformator),
			const std::string baseName = "timesteps/timestep_",
			const char* timeFmt = "%.4f",
			const std::string ext = ".dat",
			const bool posCanChange = false
		);

		/**
		 * @brief constructor with c-strings
		 *
		 * @param fcts			functions as comma-separated c-string
		 * @param subsets		subsets as comma-separated c-string
		 * @param approx		approximation space
		 * @param transformator	transformator object (default: Transformator)
		 * @param baseName		base file name for V_m file input (default: "timesteps/timestep_")
		 * @param timeFmt		format of time in file names (default: "%.4f")
		 * @param ext			file extension of input files (default: ".dat")
		 * @param posCanChange	whether coordinates in input files can change (default: false)
		 */
		VDCC_BG_VM2UG_NEURON
		(
			const char* fcts,
			const char* subsets,
			SmartPtr<ApproximationSpace<TDomain> > approx,
			SmartPtr<Transformator> transformator = make_sp(new Transformator),
			const std::string baseName = "timesteps/timestep_",
			const char* timeFmt = "%.4f",
			const std::string ext = ".dat",
			const bool posCanChange = false
		);

		/// destructor
		virtual ~VDCC_BG_VM2UG_NEURON();

		// inherited from BorgGraham
		virtual void init(number time);

		// update internal time if necessary
		virtual void update_time(number newTime);

		// update membrane potential
		virtual void update_potential(side_t* elem);

		/// @copydoc IMembraneTransporter::print_units()
		virtual void print_units() const;

		// set the transformator
		inline void set_transformator(SmartPtr<Transformator> transformator) {
			this->m_NrnInterpreter = transformator;
		}

		// set the mapper
		inline void set_mapper(SmartPtr<Vm2uG<std::string> > mapper) {
			this->m_vmProvider = mapper;
		}

	private:
		/// whether this channel disposes of an inactivating gate
		bool has_hGate() {return this->m_channelType == VDCC_BG<TDomain>::BG_Ntype
								|| this->m_channelType == VDCC_BG<TDomain>::BG_Ttype;}

	private:
		std::string m_tFmt;				//!< time format for the membrane potential files
		number m_vmTime;
};
#endif



/// Borg Graham type VGCCs with UserData membrane potential supply.
/** This class is a specialization of the Borg-Graham interface.
 *	It supplies the channel with the necessary membrane potential values by a UserData object,
 *	i.e. constant UserData or UserData provided by a lua function.
**/
template<typename TDomain>
class VDCC_BG_UserData : public VDCC_BG<TDomain>
{
	public:
		static const int dim = TDomain::dim;	//!< world dimension

	protected:
		typedef typename GeomObjBaseTypeByDim<TDomain::dim>::base_obj_type elem_t;
		typedef typename elem_t::side side_t;

		using VDCC_BG<TDomain>::R;			//!< universal gas constant
		using VDCC_BG<TDomain>::T;			//!< temperature (310K)
		using VDCC_BG<TDomain>::F;			//!< Faraday constant
		using VDCC_BG<TDomain>::has_hGate;	//!< Faraday constant


	public:
		/**
		 * @brief constructor with vectors
		 *
		 * @param fcts		functions as vector of string
		 * @param subsets	subsets as vector of string
		 * @param approx	approximation space
		 */
		VDCC_BG_UserData
		(
			const std::vector<std::string>& fcts,
			const std::vector<std::string>& subsets,
			SmartPtr<ApproximationSpace<TDomain> > approx
		);

		/**
		 * @brief constructor with c-strings
		 *
		 * @param fcts		functions as comma-separated c-string
		 * @param subsets	subsets as comma-separated c-string
		 * @param approx	approximation space
		 */
		VDCC_BG_UserData
		(
			const char* fcts,
			const char* subsets,
			SmartPtr<ApproximationSpace<TDomain> > approx
		);

		/// destructor
		virtual ~VDCC_BG_UserData();

		/// adding potential information for pumps/channels in membrane
		void set_potential_function(const char* name);
		void set_potential_function(const number value);
		void set_potential_function(SmartPtr<CplUserData<number, dim> > spPotFct);

		/// @copydoc VDCC_BG<TDomain>::update_potential()
		virtual void update_potential(side_t* elem);

	private:
		SmartPtr<CplUserData<number,dim> > m_spPotential;		//!< the UserData for potential
		bool m_bIsConstData;
};

///@}


} // namespace neuro_collection
} // namespace ug

#include "vdcc_bg_impl.h"

#endif // __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__VDCC_BG_H__
