/*
 * vdcc_bg.cpp
 *
 *  Created on: 05.02.2013
 *      Author: mbreit
 */

#include "vdcc_bg.h"
#include "lib_disc/spatial_disc/elem_disc/inner_boundary/inner_boundary.h"	// InnerBoundaryConstants

namespace ug{
namespace neuro_collection{


template<typename TDomain>
VDCC_BG<TDomain>::VDCC_BG
(
	const std::vector<std::string>& fcts,
	const std::vector<std::string>& subsets,
	SmartPtr<ApproximationSpace<TDomain> > approx
)
: IMembraneTransporter(fcts),
  R(8.314), T(310.0), F(96485.0),
  m_dom(approx->domain()), m_mg(m_dom->grid()), m_dd(approx->dof_distribution(GridLevel())),
  m_sh(m_dom->subset_handler()), m_aaPos(m_dom->position_accessor()), m_vSubset(subsets),
  m_gpMGate(3.4, -21.0, 1.5), m_gpHGate(-2.0, -40.0, 75.0), m_time(0.0), m_oldTime(0.0),
  m_perm(3.8e-19), m_mp(2), m_hp(1), m_channelType(BG_Ntype), m_initiated(false)
{
	after_construction();
}

template<typename TDomain>
VDCC_BG<TDomain>::VDCC_BG
(
	const char* fcts,
	const char* subsets,
	SmartPtr<ApproximationSpace<TDomain> > approx
)
: IMembraneTransporter(fcts),
  R(8.314), T(310.0), F(96485.0),
  m_dom(approx->domain()), m_mg(m_dom->grid()), m_dd(approx->dof_distribution(GridLevel())),
  m_sh(m_dom->subset_handler()), m_aaPos(m_dom->position_accessor()), m_vSubset(TokenizeString(subsets)),
  m_gpMGate(3.4, -21.0, 1.5), m_gpHGate(-2.0, -40.0, 75.0), m_time(0.0), m_oldTime(0.0),
  m_perm(3.8e-19), m_mp(2), m_hp(1), m_channelType(BG_Ntype), m_initiated(false)
{
	after_construction();
}

template<typename TDomain>
void VDCC_BG<TDomain>::after_construction()
{
	// process subsets

	//	remove white space
	for (size_t i = 0; i < m_vSubset.size(); ++i)
		RemoveWhitespaceFromString(m_vSubset[i]);

	//	if no subset passed, clear subsets
	if (m_vSubset.size() == 1 && m_vSubset[0].empty())
		m_vSubset.clear();

	//	if subsets passed with separator, but not all tokens filled, throw error
	for (size_t i = 0; i < m_vSubset.size(); ++i)
	{
		if (m_vSubset[i].empty())
		{
			UG_THROW("Error while setting subsets in " << name() << ": passed "
					 "subset string lacks a subset specification at position "
					 << i << "(of " << m_vSubset.size()-1 << ")");
		}
	}

	// attach attachments
	if (m_mg->template has_attachment<side_t>(this->m_MGate))
		UG_THROW("Attachment necessary for Borg-Graham channel dynamics "
				 "could not be made, since it already exists.");
	m_mg->template attach_to<side_t>(this->m_MGate);

	if (has_hGate())
	{
		if (m_mg->template has_attachment<side_t>(this->m_HGate))
			UG_THROW("Attachment necessary for Borg-Graham channel dynamics "
					 "could not be made, since it already exists.");
		m_mg->template attach_to<side_t>(this->m_HGate);
	}

	if (m_mg->template has_attachment<side_t>(this->m_Vm))
		UG_THROW("Attachment necessary for Borg-Graham channel dynamics "
				 "could not be made, since it already exists.");
	m_mg->template attach_to<side_t>(this->m_Vm);


	// create attachment accessors
	m_aaMGate = Grid::AttachmentAccessor<side_t, ADouble>(*m_mg, m_MGate);
	if (has_hGate()) m_aaHGate = Grid::AttachmentAccessor<side_t, ADouble>(*m_mg, m_HGate);
	m_aaVm = Grid::AttachmentAccessor<side_t, ADouble>(*m_mg, m_Vm);
}


template<typename TDomain>
VDCC_BG<TDomain>::~VDCC_BG()
{
	m_mg->template detach_from<side_t>(this->m_MGate);
	m_mg->template detach_from<side_t>(this->m_HGate);
	m_mg->template detach_from<side_t>(this->m_Vm);
}


template<typename TDomain>
number VDCC_BG<TDomain>::calc_gating_start(GatingParams& gp, number Vm)
{
	return 1.0 / (1.0 + exp(-gp.z * (Vm - gp.V_12) * 1e-3*F/(R*T)));
}


template<typename TDomain>
void VDCC_BG<TDomain>::calc_gating_step(GatingParams& gp, number Vm, number dt, number& currVal)
{
	const number gp_inf = calc_gating_start(gp,Vm);

	// forward step: implicit
	if (dt >= 0)
	{
		// For calculating the next gating step it is recommended to use a time step
		// size no larger than 1e-5s = 1e-2ms in order to meet a sufficient accuracy.
		if (dt > 1e-2)
		{
			number vdcc_dt = 1e-2;
			number t0 = 0.0;

			// loop intermediate time steps until the current intermediate time point
			// is the final time point dt
			while (t0 < dt)
			{
				// compute next time point
				number t = t0 + vdcc_dt;

				// check if out of bounds, if yes:
				// set to final time point and adjust step size accordingly
				if (t > dt)
				{
					t = dt;
					vdcc_dt = dt - t0;
				}

				// compute next gating value
				currVal = (gp.tau_0*currVal + vdcc_dt*gp_inf) / (gp.tau_0 + vdcc_dt);

				// save new time as current time
				t0 = t;
			}
		}
		// sufficiently small time step size already specified
		else
			currVal = (gp.tau_0*currVal + dt*gp_inf) / (gp.tau_0 + dt);
	}

	// backward step: explicit
	else
	{
		// For calculating the next gating step it is recommended to use a time step
		// size no larger than 1e-5s = 1e-2ms in order to meet a sufficient accuracy.
		if (dt < -1e-2)
		{
			number vdcc_dt = -1e-2;
			number t0 = 0.0;

			// loop intermediate time steps until the current intermediate time point
			// is the final time point dt
			while (t0 > dt)
			{
				// compute next time point
				number t = t0 + vdcc_dt;

				// check if out of bounds, if yes:
				// set to final time point and adjust step size accordingly
				if (t < dt)
				{
					t = dt;
					vdcc_dt = dt - t0;
				}

				// compute next gating value
				currVal += vdcc_dt/gp.tau_0 * (gp_inf - currVal);

				// save new time as current time
				t0 = t;
			}
		}
		// sufficiently small time step size already specified
		currVal += dt/gp.tau_0 * (gp_inf - currVal);
	}
}


template<typename TDomain>
void VDCC_BG<TDomain>::calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const
{
	side_t* elem = dynamic_cast<side_t*>(e);
	if (!elem) UG_THROW("OneSidedBorgGrahamFV1 fluxDensityFunction called with the wrong type of element.");

	number gating = pow(m_aaMGate[elem], m_mp);
	if (has_hGate()) gating *= pow(m_aaHGate[elem], m_hp);

	// flux derived from Goldman-Hodgkin-Katz equation,
	number maxFlux;
	number vm = m_aaVm[elem];
	number caCyt = u[_CCYT_];		// cytosolic Ca2+ concentration
	number caExt = u[_CEXT_];		// extracellular Ca2+ concentration

	// near V_m == 0: approximate by first order Taylor to avoid relative errors and div-by-0
	if (fabs(vm) < 1e-8) maxFlux = m_perm * ((caExt - caCyt) - F/(R*T) * (caExt + caCyt)*vm);
	else maxFlux = -m_perm * 2*F/(R*T) * vm * (caExt - caCyt*exp(2*F/(R*T)*vm)) / (1.0 - exp(2*F/(R*T)*vm));

	flux[0] = gating * maxFlux;

	//UG_LOGN(std::setprecision(std::numeric_limits<long double>::digits10 + 1) << "VDCC flux:" << flux[0]);

	//UG_COND_THROW(flux[0] != flux[0],
	//	"VDCC NaN: gating = " << gating << ", maxFlux = " << maxFlux);
}


template<typename TDomain>
void VDCC_BG<TDomain>::calc_flux_deriv(const std::vector<number>& u, GridObject* e, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const
{
	side_t* elem = dynamic_cast<side_t*>(e);
	if (!elem) UG_THROW("OneSidedBorgGrahamFV1 fluxDensityFunction called with the wrong type of element.");

	number gating = pow(m_aaMGate[elem], m_mp);
	if (has_hGate()) gating *= pow(m_aaHGate[elem], m_hp);

	number dMaxFlux_dCyt, dMaxFlux_dExt;
	number vm = m_aaVm[elem];

	// near V_m == 0: approximate by first order Taylor to avoid relative errors and div-by-0
	if (fabs(vm) < 1e-8)
	{
		dMaxFlux_dCyt =  m_perm * (-1.0 - F/(R*T) * vm);
		dMaxFlux_dExt =  m_perm * (1.0 - F/(R*T) * vm);
	}
	else
	{
		dMaxFlux_dCyt = m_perm * 2*F/(R*T) * vm / (exp(-2*F/(R*T)*vm) - 1.0);
		dMaxFlux_dExt = m_perm * 2*F/(R*T) * vm / (exp(2*F/(R*T)*vm) - 1.0);
	}

	size_t i = 0;
	if (!has_constant_value(_CCYT_))
	{
		flux_derivs[0][i].first = local_fct_index(_CCYT_);
		flux_derivs[0][i].second = gating * dMaxFlux_dCyt;
		i++;
	}
	if (!has_constant_value(_CEXT_))
	{
		flux_derivs[0][i].first = local_fct_index(_CEXT_);
		flux_derivs[0][i].second = gating * dMaxFlux_dExt;
	}
}


template<typename TDomain>
size_t VDCC_BG<TDomain>::n_dependencies() const
{
	size_t n = 2;
	if (has_constant_value(_CCYT_)) n--;
	if (has_constant_value(_CEXT_)) n--;

	return n;
}


template<typename TDomain>
size_t VDCC_BG<TDomain>::n_fluxes() const
{
	return 1;
}


template<typename TDomain>
const std::pair<size_t,size_t> VDCC_BG<TDomain>::flux_from_to(size_t flux_i) const
{
	size_t from, to;
	if (allows_flux(_CCYT_)) to = local_fct_index(_CCYT_); else to = InnerBoundaryConstants::_IGNORE_;
	if (allows_flux(_CEXT_)) from = local_fct_index(_CEXT_); else from = InnerBoundaryConstants::_IGNORE_;

	return std::pair<size_t, size_t>(from, to);
}


template<typename TDomain>
const std::string VDCC_BG<TDomain>::name() const
{
	return std::string("VDCC_BG");
}


template<typename TDomain>
void VDCC_BG<TDomain>::check_supplied_functions() const
{
	// Check that not both, inner and outer calcium concentrations are not supplied;
	// in that case, calculation of a flux would be of no consequence.
	if (!allows_flux(_CCYT_) && !allows_flux(_CEXT_))
	{
		UG_THROW("Supplying neither cytosolic nor extracellular calcium concentrations is not allowed.\n"
				"This would mean that the flux calculation would be of no consequence\n"
				"and this channel would not do anything.");
	}
}


template<typename TDomain>
void VDCC_BG<TDomain>::print_units() const
{
	std::string nm = name();
	size_t n = nm.size();
	UG_LOG(std::endl);
	UG_LOG("+------------------------------------------------------------------------------+"<< std::endl);
	UG_LOG("|  Units used in the implementation of " << nm << std::string(n>=40?0:40-n, ' ') << "|" << std::endl);
	UG_LOG("|------------------------------------------------------------------------------|"<< std::endl);
	UG_LOG("|    Input                                                                     |"<< std::endl);
	UG_LOG("|      [Ca_cyt]  mM (= mol/m^3)                                                |"<< std::endl);
	UG_LOG("|      [Ca_ext]  mM (= mol/m^3)                                                |"<< std::endl);
	UG_LOG("|      V_m       V                                                             |"<< std::endl);
	UG_LOG("|                                                                              |"<< std::endl);
	UG_LOG("|    Output                                                                    |"<< std::endl);
	UG_LOG("|      Ca flux   mol/s                                                         |"<< std::endl);
	UG_LOG("+------------------------------------------------------------------------------+"<< std::endl);
	UG_LOG(std::endl);
}


template<typename TDomain>
void VDCC_BG<TDomain>::set_permeability(const number perm)
{
	m_perm = perm;
}


template<typename TDomain>
void VDCC_BG<TDomain>::init(number time)
{
	m_time = time;
	m_initTime = time;

	typedef typename DoFDistribution::traits<side_t>::const_iterator itType;
	SubsetGroup ssGrp;
	try { ssGrp = SubsetGroup(m_dom->subset_handler(), this->m_vSubset);}
	UG_CATCH_THROW("Subset group creation failed.");

	for (std::size_t si = 0; si < ssGrp.size(); si++)
	{
		itType iterBegin = m_dd->template begin<side_t>(ssGrp[si]);
		itType iterEnd = m_dd->template end<side_t>(ssGrp[si]);

		for (itType iter = iterBegin; iter != iterEnd; ++iter)
		{
			// get potential data
			update_potential(*iter);
			number vm = m_aaVm[*iter];

			// calculate corresponding start condition for gates
			m_aaMGate[*iter] = calc_gating_start(m_gpMGate, 1e3*vm);
			if (has_hGate()) m_aaHGate[*iter] = calc_gating_start(m_gpHGate, 1e3*vm);
		}
	}

	this->m_initiated = true;
}


template<typename TDomain>
void VDCC_BG<TDomain>::update_gating(side_t* elem)
{
	if (!this->m_initiated)
		UG_THROW("Borg-Graham not initialized.\n"
				  << "Do not forget to do so before any updates by calling init(initTime).");

	// set new gating particle values
	number dt = 1e3*(m_time - m_oldTime);	// calculating in ms
	calc_gating_step(m_gpMGate, 1e3*m_aaVm[elem], dt, m_aaMGate[elem]);
	if (has_hGate())
		calc_gating_step(m_gpHGate, 1e3*m_aaVm[elem], dt, m_aaHGate[elem]);
}

template<typename TDomain>
void VDCC_BG<TDomain>::update_time(const number newTime)
{
	if (newTime != m_time)
	{
		m_oldTime = m_time;
		m_time = newTime;
	}
};

template<typename TDomain>
void VDCC_BG<TDomain>::prep_timestep
(
    number future_time, const number time, VectorProxyBase* upb
)
{
    // initiate if this has not already been done (or init again; stationary case)
    if (!m_initiated || future_time == m_initTime)
        init(time);

	update_time(future_time);
    bool backwardsStep = m_time < m_oldTime;

	SubsetGroup ssGrp;
    try { ssGrp = SubsetGroup(m_dom->subset_handler(), this->m_vSubset);}
    UG_CATCH_THROW("Subset group creation failed.");

    typedef typename DoFDistribution::traits<side_t>::const_iterator it_type;


    // TODO: Think about updating only on the base level and then propagating upwards.
    //       Typically, the potential does not need very fine resolution.
    //       This would save a lot of work for very fine surface levels.
    for (std::size_t si = 0; si < ssGrp.size(); si++)
    {
        // loop sides and update potential and then gatings
        it_type it = m_dd->begin<side_t>(ssGrp[si]);
        it_type it_end = m_dd->end<side_t>(ssGrp[si]);
        if (backwardsStep)
        {
			for (; it != it_end; ++it)
			{
				update_gating(*it);
				update_potential(*it);
			}
        }
        else
        {
			for (; it != it_end; ++it)
			{
				update_potential(*it);
				update_gating(*it);
			}
        }
    }
}



// explicit template specializations
#ifdef UG_DIM_1
	template class VDCC_BG<Domain1d>;
#endif
#ifdef UG_DIM_2
	template class VDCC_BG<Domain2d>;
#endif
#ifdef UG_DIM_3
	template class VDCC_BG<Domain3d>;
#endif

} // neuro_collection
} // namespace ug
