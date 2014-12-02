/*
 * one_sided_borg_graham_fv1.h
 *
 *  Created on: 05.02.2013
 *      Author: mbreit
 */

#ifndef __H__UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__ONE_SIDED_BORG_GRAHAM_FV1__
#define __H__UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__ONE_SIDED_BORG_GRAHAM_FV1__


#include "lib_grid/lg_base.h"
#include <vector>
#include <stdio.h>
#include "../plugins/experimental/membrane_potential_mapping/vm2ug.h"
#include "one_sided_membrane_transport_fv1.h"
#include <locale>	// for control over the decimal separator (point instead of comma, please!)
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/hfv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"

namespace ug
{
namespace neuro_collection
{

/// Class interface for Borg Graham type VGCCs of the plasma membrane.
/** This class is an interface the Borg-Graham-type voltage-gated calcium channels
 *	in the plasma membrane (see chapter 8 of "Interpretations of data and mechanisms
 *	for hippocampal pyramidal cell models", Borg-Graham (1998) in Cerebral Cortex,
 *	Vol. 13: Models of Cortical Circuits.
 *
 *	The unknowns of the discretization (the so-called "particles" of the channel) are
 *	not added to the system unknowns, but held separately as vertex attachments, since
 *	they are governed by ODEs ans are therefore merely updated by an explicit time
 *	schema once per time step. This will have to be done by the user (in the lua script)
 *	via the method update_gating() before the execution of every time step assembling.
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
 *	they _can_ reimplement the virtual method
 *	- void init(number time).
 *
 *
 *	The units required for this discretization are:
 * 		V_m	: mV		membrane voltage
 * 		t	: s			time
 *  	f	: mol*s^-1	ionic flux
 *
**/
template<typename TDomain>
class OneSidedBorgGrahamFV1 : public OneSidedMembraneTransportFV1<TDomain>
{
	public:
		/// channel types N, L and T
		enum {BG_Ntype, BG_Ltype, BG_Ttype};

	protected:
		using OneSidedMembraneTransportFV1<TDomain>::R;		//!< universal gas constant
		using OneSidedMembraneTransportFV1<TDomain>::T;		//!< temperature (310K)
		using OneSidedMembraneTransportFV1<TDomain>::F;		//!< Faraday constant
		using OneSidedMembraneTransportFV1<TDomain>::dim;	//!< world dimension

		const number RHO_BG;			//!< default channel density in the membrane

	private:
		typedef OneSidedBorgGrahamFV1<TDomain> this_type;
		typedef typename DependentNeumannBoundaryFV1<TDomain>::NFluxCond NFluxCond;
		typedef typename DependentNeumannBoundaryFV1<TDomain>::NFluxDerivCond NFluxDerivCond;
		typedef typename GeomObjBaseTypeByDim<dim>::base_obj_type elem_t;
		typedef typename elem_t::side side_t;

	private:
		/// holds the paramters of a channel type
		struct GatingParams
		{
			GatingParams(number _z, number _v, number _t) : z(_z), V_12(_v), tau_0(_t){};
			number z;
			number V_12;
			number tau_0;
		};

	public:
		/// constructor
		OneSidedBorgGrahamFV1(const char* functions, const char* subsets, ApproximationSpace<TDomain>& approx);

		/// destructor
		virtual	~OneSidedBorgGrahamFV1() {};

		/// sets the channel type
		template<int TType> void set_channel_type();

		/// initializes the defined channel type
		/** During the initialization, the necessary attachments are attached to the vertices
		 *	and their values calculated by the equilibrium state for the start membrane potential.
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
		 * This method needs to be called before ionic_current().
		 * @param newTime new point in time
		 */
		void update_gating(side_t* elem);

		/// provides the ionic current (mol*s^-1) at a given vertex
		number ionic_current(side_t* e);

		/// updates internal time if necessary
		virtual void update_time(number newTime) {m_oldTime = m_time; m_time = newTime;};

		// inherited from IElemDisc
		virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);

		template<typename TElem>
		void prep_timestep_elem(const number time, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		// inherited from FV1MyNeumannBoundaryElemDisc
		virtual bool fluxDensityFct(const std::vector<LocalVector::value_type>& u, GridObject* e, const MathVector<dim>& coords, int si, NFluxCond& fc);
		virtual bool fluxDensityDerivFct(const std::vector<LocalVector::value_type>& u, GridObject* e, const MathVector<dim>& coords, int si, NFluxDerivCond& fdc);

	protected:
		/// calculates the equilibrium state of a gating "particle"
		/** The calculation is done with respect to the given gating paramters set (which represents
		 *	one gating "particle") and the given membrane potential.
		**/
		number calc_gating_start(GatingParams& gp, number Vm);

		/// calculates the next state of a gating "particle"
		/** The calculation is done with respect to the given gating paramters set (which represents
		 *	one gating "particle"), the current value of this gating particle as well as the given membrane
		 *	potential. The new value represents the "particle" state at current time + dt.
		**/
		void calc_gating_step(GatingParams& gp, number Vm, number dt, number& currVal);


	protected:
		/// whether this channel disposes of an inactivating gate
		bool has_hGate() {return this->m_channelType == OneSidedBorgGrahamFV1<TDomain>::BG_Ntype
								|| this->m_channelType == OneSidedBorgGrahamFV1<TDomain>::BG_Ttype;}

	protected:
		SmartPtr<TDomain> m_dom;					//!< underlying domain
		SmartPtr<Grid> m_mg;						//!< underlying multigrid
		SmartPtr<DoFDistribution> m_dd;				//!< underlying surface dof distribution
		ConstSmartPtr<MGSubsetHandler> m_sh;		//!< underlying subset handler
		typename TDomain::position_accessor_type& m_aaPos;	//!< underlying position accessor

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

		number m_lambda;							//!< channel conductivity (C/(V*s))
		int m_mp, m_hp;								//!< powers for gating parameters

		int m_channelType;							//!< channel type

		bool m_initiated;							//!< indicates whether channel has been initialized by init()

	private:
		struct RegisterFV1
		{
				RegisterFV1(this_type* pThis) : m_pThis(pThis){}
				this_type* m_pThis;
				template< typename TElem > void operator()(TElem&)
				{
					if (m_pThis->m_bNonRegularGrid)
						m_pThis->register_fv1_func<TElem, HFV1ManifoldGeometry<TElem, dim> >();
					else
						m_pThis->register_fv1_func<TElem, FV1ManifoldGeometry<TElem, dim> >();

				}
		};

		void register_all_fv1_funcs();

		template <typename TElem, typename TFVGeom>
		void register_fv1_func();

	private:
		bool m_bNonRegularGrid;
};



/// Borg Graham type VGCCs with Vm2uG membrane potential supply.
/** This class is a specialization of the Borg-Graham interface.
 *	It supplies the channel with the necessary membrane potential values by a Vm2uG object,
 *	that has to be fed file names for the files containing the V_m values produced by Neuron.
**/
template<typename TDomain>
class OneSidedBorgGrahamFV1WithVM2UG : public OneSidedBorgGrahamFV1<TDomain>
{
	protected:
		using OneSidedBorgGrahamFV1<TDomain>::R;		//!< universal gas constant
		using OneSidedBorgGrahamFV1<TDomain>::T;		//!< temperature (310K)
		using OneSidedBorgGrahamFV1<TDomain>::F;		//!< Faraday constant
		using OneSidedBorgGrahamFV1<TDomain>::has_hGate;//!< Faraday constant
		typedef typename GeomObjBaseTypeByDim<TDomain::dim>::base_obj_type elem_t;
		typedef typename elem_t::side side_t;

	public:
		typedef Vm2uG<std::string> vmProvType;

	public:
		/// constructor
		OneSidedBorgGrahamFV1WithVM2UG(const char* functions,
							const char* subsets,
							ApproximationSpace<TDomain>& approx,
							const std::string baseName = "timesteps/timestep_",
							const char* timeFmt = "%.4f",
							const std::string ext = ".dat",
							const bool posCanChange = false)
			: OneSidedBorgGrahamFV1<TDomain>(functions, subsets, approx),
			  m_vmProvider(baseName, ext, !posCanChange), m_tFmt(timeFmt),
			  m_fileInterval(0.0), m_fileOffset(0.0) {};

		/// destructor
		virtual ~OneSidedBorgGrahamFV1WithVM2UG() {};

		// inherited from BorgGraham
		virtual void init(number time);

		// update membrane potential
		virtual void update_potential(side_t* elem);

		// update internal time if necessary
		virtual void update_time(number newTime);

		// set times for which files with potential values are available
		void set_file_times(const number fileInterval, const number fileOffset = 0.0)
		{
			m_fileInterval = fileInterval;
			m_fileOffset = fileOffset;
		}

	private:
		vmProvType m_vmProvider;		//!< the Vm2uG object
		std::string m_tFmt;				//!< time format for the membrane potential files
		number m_fileInterval;			//!< intervals in which voltage files are available
		number m_fileOffset;				//!< offset of time intervals for which voltage files are available

		std::string m_timeAsString;
};

#ifdef MPMNEURON
/// Borg Graham type VGCCs with Vm2uG membrane potential supply by NEURON.
/** This class is a specialization of the Borg-Graham interface.
 *	It supplies the channel with the necessary membrane potential values by a Vm2uG object,
 *	that has to be fed file names for the files containing the V_m values produced by Neuron.
**/
template<typename TDomain>
class OneSidedBorgGrahamFV1WithVM2UGNEURON : public OneSidedBorgGrahamFV1<TDomain>
{
private:
	SmartPtr<Transformator> m_NrnInterpreter;
	std::string timeAsString = "";
	protected:
		using OneSidedBorgGrahamFV1<TDomain>::R;		//!< universal gas constant
		using OneSidedBorgGrahamFV1<TDomain>::T;		//!< temperature (310K)
		using OneSidedBorgGrahamFV1<TDomain>::F;		//!< Faraday constant

		typedef typename GeomObjBaseTypeByDim<TDomain::dim>::base_obj_type elem_t;
		typedef typename elem_t::side side_t;

	public:
		typedef Vm2uG<std::string> vmProvType;

	public:
		/// constructor
		OneSidedBorgGrahamFV1WithVM2UGNEURON(const char* functions,
							const char* subsets,
							ApproximationSpace<TDomain>& approx,
							SmartPtr<Transformator> transformator = Transformator(),
							const std::string baseName = "timesteps/timestep_",
							const char* timeFmt = "%.4f",
							const std::string ext = ".dat",
							const bool posCanChange = false)
			: OneSidedBorgGrahamFV1<TDomain>(functions, subsets, approx),
			  m_NrnInterpreter(transformator), m_tFmt(timeFmt), m_vmTime(0.0) {
		};

		/// destructor
		virtual ~OneSidedBorgGrahamFV1WithVM2UGNEURON() {};

		// inherited from BorgGraham
		virtual void init(number time);

		// update internal time if necessary
		virtual void update_time(number newTime);

		// update membrane potential
		virtual void update_potential(side_t* elem);

		// set the transformator
		inline void set_transformator(SmartPtr<Transformator> transformator) {
			this->m_NrnInterpreter = transformator;
		}

	private:
		/// whether this channel disposes of an inactivating gate
		bool has_hGate() {return this->m_channelType == OneSidedBorgGrahamFV1<TDomain>::BG_Ntype
								|| this->m_channelType == OneSidedBorgGrahamFV1<TDomain>::BG_Ttype;}

	private:
		vmProvType m_vmProvider;		//!< the Vm2uG object
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
class OneSidedBorgGrahamFV1WithUserData : public OneSidedBorgGrahamFV1<TDomain>
{
	protected:
		using OneSidedBorgGrahamFV1<TDomain>::R;		//!< universal gas constant
		using OneSidedBorgGrahamFV1<TDomain>::T;		//!< temperature (310K)
		using OneSidedBorgGrahamFV1<TDomain>::F;		//!< Faraday constant
		using OneSidedBorgGrahamFV1<TDomain>::has_hGate;//!< Faraday constant
		typedef typename GeomObjBaseTypeByDim<TDomain::dim>::base_obj_type elem_t;
		typedef typename elem_t::side side_t;

		///	world dimension
		using OneSidedBorgGrahamFV1<TDomain>::dim;	//!< world dimension

	public:
		/// constructor
		OneSidedBorgGrahamFV1WithUserData(const char* functions,
							const char* subsets,
							ApproximationSpace<TDomain>& approx)
			: OneSidedBorgGrahamFV1<TDomain>(functions, subsets, approx), m_bIsConstData(false) {};

		/// destructor
		virtual ~OneSidedBorgGrahamFV1WithUserData() {};

		/// adding potential information for pumps/channels in membrane
		void set_potential_function(const char* name);
		void set_potential_function(const number value);
		void set_potential_function(SmartPtr<CplUserData<number, dim> > spPotFct);

		// update membrane potential
		virtual void update_potential(side_t* elem);

	private:
		SmartPtr<CplUserData<number,dim> > m_spPotential;		//!< the UserData for potential
		bool m_bIsConstData;
};


} // namespace neuro_collection
} // namespace ug



#include "one_sided_borg_graham_fv1_impl.h"

#endif /* __H__UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__ONE_SIDED_BORG_GRAHAM_FV1__ */
