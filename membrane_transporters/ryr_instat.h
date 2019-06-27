/*
 *	Alternative discretization for the RyR calcium channel in the ER membrane with time delay
 *
 *  Created on: Nov 12, 2015
 *      Author: marcuskessler, mbreit
 */

#ifndef UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__RYR_INSTAT_H
#define UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__RYR_INSTAT_H

#include "membrane_transporter_interface.h"
#include "lib_disc/spatial_disc/elem_disc/inner_boundary/inner_boundary.h"


namespace ug {
namespace neuro_collection {


///@addtogroup plugin_neuro_collection
///@{


/// Discretization for the RyR calcium channel in the ER membrane
/**
 * This class implements the IMembraneTransporter interface to provide flux densities
 * and their derivatives for the Keizer & Levine (1996) RyR model (instationary).
 *
 * Units used in the implementation of this channel:
 * [Ca_cyt]  mM (= mol/m^3)
 * [Ca_er]   mM (= mol/m^3)
 *
 * Ca flux   mol/s
 */

template<typename TDomain>
class RyRinstat : public IMembraneTransporter
{
	public:
		enum{_CCYT_=0, _CER_};

		static const int dim = TDomain::dim;	//!< world dimension

		typedef typename GeomObjBaseTypeByDim<dim>::base_obj_type elem_t;
		typedef typename elem_t::side side_t;

	protected:
		const number R;			// universal gas constant
		const number T;			// temperature
		const number F;			// Faraday constant

		const number KAplus;		// calcium binding (C1 --> O1)
		const number KBplus;		// calcium binding (O1 --> O2)
		const number KCplus;		// O1 --> C2
		const number KAminus;		// C1 <-- O1
		const number KBminus;		// O1 <-- O2
		const number KCminus;		// O1 <-- C2
		const number MU_RYR;	// RyR channel conductance

		const number REF_CA_ER;	// reference endoplasmatic Ca2+ concentration (for conductances)

	public:
		/// @copydoc IMembraneTransporter::IMembraneTransporter(const std::vector<std::string)
		RyRinstat(const std::vector<std::string>& fcts, const std::vector<std::string>& subsets, SmartPtr<ApproximationSpace<TDomain> > approx);

		/// @copydoc IMembraneTransporter::IMembraneTransporter()
		RyRinstat(const char* fcts, const char* subsets, SmartPtr<ApproximationSpace<TDomain> > approx);

		/// @copydoc IMembraneTransporter::IMembraneTransporter()
		virtual ~RyRinstat();

		/// @copydoc IMembraneTransporter::prepare_timestep()
		virtual void prepare_timestep(number future_time, const number time, VectorProxyBase* upb);

		/// @copydoc IMembraneTransporter::calc_flux()
		virtual void calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const;

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

	protected:
		// constructing directives to be called from every constructor
		void construct(const std::vector<std::string>& subsets, SmartPtr<ApproximationSpace<TDomain> > approx);

		// init gating variables to equilibrium
		void init(number time, VectorProxyBase* upb);

		// calculate open probability at grid object
		template <typename TBaseElem>
		number open_prob(GridObject* o) const;

		number open_prob_for_grid_object(GridObject* o) const;

	protected:
		SmartPtr<TDomain> m_dom;					//!< underlying domain
		SmartPtr<MultiGrid> m_mg;					//!< underlying multigrid
		SmartPtr<DoFDistribution> m_dd;				//!< underlying surface dof distribution

		std::vector<size_t> m_vSubset;				//!< subset indices this mechanism works on

		ADouble m_aO2;								//!< proportion of channels currently in state O2
		ADouble m_aC1;								//!< proportion of channels currently in state C1
		ADouble m_aC2;								//!< proportion of channels currently in state C2
		ANumber m_aOavg; 							//!< average open probability between two consecutive time steps
		ANumber m_aCaOld; 							//!< old calcium value

		Grid::AttachmentAccessor<Vertex, ADouble> m_aaO2;		//!< accessor for channel states
		Grid::AttachmentAccessor<Vertex, ADouble> m_aaC1;		//!< accessor for channel states
		Grid::AttachmentAccessor<Vertex, ADouble> m_aaC2;		//!< accessor for channel states
		Grid::AttachmentAccessor<Vertex, ANumber> m_aaOavg;		//!< accessor for avg. open prob.
		Grid::AttachmentAccessor<Vertex, ANumber> m_aaCaOld;	//!< accessor for old calcium value

		number m_time;								//!< current time
		number m_initTime;							//!< time of initialization
		bool m_initiated;							//!< indicates whether channel has been initialized by init()
};

///@}


} // namespace neuro_collection
} // namespace ug

#endif // UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__RYR_INSTAT_H

