/*
 * vdcc_bg.h
 *
 *  Created on: 05.02.2013
 *      Author: mbreit
 */

#ifndef __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__VDCC_BG_H__
#define __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__VDCC_BG_H__

#include "../membrane_transporter_interface.h"

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
 *  	- Internally, all membrane potentials are attached to elements in [V] for
 *  	  the element discretization!
 *
 *  	- VDCC_BG & VDCC_BG_UserData use [ms] and [mV] in gating and flux calculations!
 *  	  This is due to the use of the original gating parameter sets by Borg-Graham
 *  	  with tau_0 in [ms] and V_12 in [mV].
 *  	  Note that the update potential method already takes care of this, when calling
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

		static const int dim = TDomain::dim;	//!< world dimension

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

		/// @copydoc IMembraneTransporter::prep_timestep()
		virtual void prep_timestep
		(
			number future_time, const number time, VectorProxyBase* upb
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
		virtual size_t n_dependencies() const;

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

        /**
         * @brief Sets the permeability of this channel
         * @param perm    permeability values
         */
		void set_permeability(const number perm);

        /// initializes the defined channel type
        /** During the initialization, the necessary attachments are attached to the vertices
         *  and their values calculated by the equilibrium state for the start membrane potential.
        **/
		virtual void init(number time);

		/// updates the potential values in the corresponding attachments to new time.
		/**
		 * This method needs to be called before update_gating() if potential is non-constant.
		 * @param newTime new point in time
		 */
		virtual void update_potential(side_t* elem) = 0;

		/// updates the gating parameters
		/**
		 * This method needs to be called before calc_flux().
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
		/// whether this channel has an inactivating gate
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
		number m_initTime;							//!< time of initialization
		number m_oldTime;							//!< time step before current time

		number m_perm;								//!< channel permeability [m^3/s] (= diff coeff * cross section / membrane thickness)
		int m_mp, m_hp;								//!< powers for gating parameters

		int m_channelType;							//!< channel type

		bool m_initiated;							//!< indicates whether channel has been initialized by init()
};


///@}


} // namespace neuro_collection
} // namespace ug

#include "vdcc_bg_impl.h"

#endif // __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__VDCC_BG_H__
