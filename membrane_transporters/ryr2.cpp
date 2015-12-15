/*
 *	Discretization for the RyR calcium channel in the ER membrane
 *
 *  Created on: Nov 12, 2015
 *      Author: marcuskessler
 */

#include "ryr2.h"


namespace ug{
namespace neuro_collection{

template<typename TDomain>
RyR2<TDomain>::
RyR2
(
	const std::vector<std::string>& fcts,
	SmartPtr<ApproximationSpace<TDomain> > approx
)
: IMembraneTransporter(fcts),
R(8.314), T(310.0), F(96485.0),
KAplus(1500.0e12), KBplus(1500.0e12), KCplus(1.75), MU_RYR(5.0e-11),
KAminus(28.8), KBminus(385.9), KCminus(0.1), REF_CA_ER(2.5e-1),
m_time(0.0), m_oldTime(0.0), m_initiated(false)
{}

template<typename TDomain>
RyR2<TDomain>::
RyR2
(
	const char* fcts,
	SmartPtr<ApproximationSpace<TDomain> > approx
)
: IMembraneTransporter(fcts),
R(8.314), T(310.0), F(96485.0),
KAplus(1500.0e12), KBplus(1500.0e12), KCplus(1.75), MU_RYR(5.0e-11),
KAminus(28.8), KBminus(385.9), KCminus(0.1), REF_CA_ER(2.5e-1),
m_time(0.0), m_oldTime(0.0), m_initiated(false)
{
	// save underlying multigrid
	m_mg = approx->domain()->grid();

	// save underlying surface dof distribution
	m_dd = approx->dof_distribution(GridLevel(), false);

	if (m_mg->template has_attachment<side_t>(this->m_aO1))
	UG_THROW("Attachment necessary for RyR channel dynamics "
			 "could not be made, since it already exists.");
	m_mg->template attach_to<side_t>(this->m_aO1);

	if (m_mg->template has_attachment<side_t>(this->m_aO2))
		UG_THROW("Attachment necessary for RyR channel dynamics "
				 "could not be made, since it already exists.");
	m_mg->template attach_to<side_t>(this->m_aO2);

	if (m_mg->template has_attachment<side_t>(this->m_aC1))
		UG_THROW("Attachment necessary for RyR channel dynamics "
				 "could not be made, since it already exists.");
	m_mg->template attach_to<side_t>(this->m_aC1);

	if (m_mg->template has_attachment<side_t>(this->m_aC2))
		UG_THROW("Attachment necessary for RyR channel dynamics "
				 "could not be made, since it already exists.");
	m_mg->template attach_to<side_t>(this->m_aC2);

	m_aaO1 = Grid::AttachmentAccessor<side_t, ADouble>(*m_mg, m_aO1);
	m_aaO2 = Grid::AttachmentAccessor<side_t, ADouble>(*m_mg, m_aO2);
	m_aaC1 = Grid::AttachmentAccessor<side_t, ADouble>(*m_mg, m_aC1);
	m_aaC2 = Grid::AttachmentAccessor<side_t, ADouble>(*m_mg, m_aC2);
}


template<typename TDomain>
RyR2<TDomain>::~RyR2()
{
	m_mg->template detach_from<side_t>(this->m_aO1);
	m_mg->template detach_from<side_t>(this->m_aO2);
	m_mg->template detach_from<side_t>(this->m_aC1);
	m_mg->template detach_from<side_t>(this->m_aC2);
};


template <typename TDomain>
void RyR2<TDomain>::prep_timestep(const number time, VectorProxyBase* upb)
{
	// get solution u with which to prepare time step (this code only accepts CPUAlgebra type)
	typedef CPUAlgebra::vector_type v_type;
	typedef VectorProxy<v_type> vp_type;
	vp_type* up = dynamic_cast<vp_type*>(upb);
	UG_COND_THROW(!up, "Wrong algebra type!");
	const v_type& u = up->m_v;

	// for DoF index storage
	std::vector<DoFIndex> dofIndex;

	// update time
	if (time != m_time)
	{
		m_oldTime = m_time;
		m_time = time;
	}
	number dt = m_time - m_oldTime;

	// loop sides and update potential and then gatings
	typedef typename DoFDistribution::traits<side_t>::const_iterator it_type;
	it_type it = m_dd->begin<side_t>();
	it_type it_end = m_dd->end<side_t>();

	for (; it != it_end; ++it)
	{
		number& pO2 = m_aaO2[*it];
		number& pO1 = m_aaO1[*it];
		number& pC2 = m_aaC2[*it];
		number& pC1 = m_aaC1[*it];

		// get ca_cyt and ca_er;
		// we suppose our approx space to be 1st order Lagrange (linear, DoFs in the vertices)
		// and interpolate value at the center of the element
		number ca_cyt = 0.0;
		m_dd->dof_indices(*it, _CCYT_, dofIndex, false, true);
		UG_ASSERT(dofIndex.size() > 0, "No DoF found for function " << _CCYT_ << "on element.");
		for (size_t i = 0; i < dofIndex.size(); ++i) ca_cyt += DoFRef(u, dofIndex[i]);
		ca_cyt /= dofIndex.size();

		number ca_er = 0.0;
		m_dd->dof_indices(*it, _CER_, dofIndex, false, true);
		UG_ASSERT(dofIndex.size() > 0, "No DoF found for function " << _CER_ << "on element.");
		for (size_t i = 0; i < dofIndex.size(); ++i) ca_er += DoFRef(u, dofIndex[i]);
		ca_er /= dofIndex.size();

		// implicit Euler!
		number pO1num = (1 + dt * KBminus) * ((-pC2 + 1 + dt * KCminus) * (1 + dt * KAplus * pow(ca_cyt, 3)) - pC1 * (1 + dt * KCminus)) - pO2 * (1 + dt * KCminus) * (1 + dt * KAplus * pow(ca_cyt, 3));
		number pO1den = (1 + dt * KBminus) * ((1 + dt * KCminus) * (dt * KAminus) + (dt * KCplus + 1 + dt * KCminus) * (1 + dt * KAplus * pow(ca_cyt, 3))) + (dt * KBplus * pow(ca_cyt, 4)) * (1 + dt * KCminus) * (1 + dt * KAplus * pow(ca_cyt, 3));
		pO1 = pO1num / pO1den;
		pO2 = (pO2 + dt * KBplus * pow(ca_cyt, 4) * pO1) / (1 + dt * KBminus);
		pC1 = (pC1 + dt * KAminus * pO1) / (1 + dt * KAplus * pow(ca_cyt, 3));
		pC2 = (pC2 + dt * KCplus * pO1) / (1 + dt * KCminus);
	}
}


//TODO:  add init method; init in equilibrium state!
template<typename TDomain>
void RyR2<TDomain>::init(number time, VectorProxyBase* upb)
{
	this->m_time = time;

	int test = 1;
	// get solution u with which to prepare time step (this code only accepts CPUAlgebra type)
	typedef CPUAlgebra::vector_type v_type;
	typedef VectorProxy<v_type> vp_type;
	vp_type* up = dynamic_cast<vp_type*>(upb);
	UG_COND_THROW(!up, "Wrong algebra type!");
	const v_type& u = up->m_v;

	// for DoF index storage
	std::vector<DoFIndex> dofIndex;

	typedef typename DoFDistribution::traits<side_t>::const_iterator it_type;

	it_type it = m_dd->begin<side_t>();
	it_type it_end = m_dd->end<side_t>();

	for (; it != it_end; ++it)
	{
		number& pO2 = m_aaO2[*it];
		number& pO1 = m_aaO1[*it];
		number& pC2 = m_aaC2[*it];
		number& pC1 = m_aaC1[*it];

		// get ca_cyt and ca_er;
		// we suppose our approx space to be 1st order Lagrange (linear, DoFs in the vertices)
		// and interpolate value at the center of the element
		number ca_cyt = 0.0;
		m_dd->dof_indices(*it, _CCYT_, dofIndex, false, true);
		UG_ASSERT(dofIndex.size() > 0, "No DoF found for function " << _CCYT_ << "on element.");
		for (size_t i = 0; i < dofIndex.size(); ++i) ca_cyt += DoFRef(u, dofIndex[i]);
		ca_cyt /= dofIndex.size();

		number ca_er = 0.0;
		m_dd->dof_indices(*it, _CER_, dofIndex, false, true);
		UG_ASSERT(dofIndex.size() > 0, "No DoF found for function " << _CER_ << "on element.");
		for (size_t i = 0; i < dofIndex.size(); ++i) ca_er += DoFRef(u, dofIndex[i]);
		ca_er /= dofIndex.size();

		//calculate equilibrium
		pO1 = (1.0 / (1.0 + KCplus/KCminus + (1/KAplus/KAminus*pow(ca_cyt,3)) + KBplus/KBminus*pow(ca_cyt,4)));
		pO2 = (KBplus/KBminus*pow(ca_cyt,4) / (1.0 + KCplus/KCminus + (1/KAplus/KAminus*pow(ca_cyt,3)) + KBplus/KBminus*pow(ca_cyt,4)));
		pC1 = ((1/KAplus/KAminus*pow(ca_cyt,3)) / (1.0 + KCplus/KCminus + (1/KAplus/KAminus*pow(ca_cyt,3)) + KBplus/KBminus*pow(ca_cyt,4)));
		pC2 = (KCplus/KCminus / (1.0 + KCplus/KCminus + (1/KAplus/KAminus*pow(ca_cyt,3)) + KBplus/KBminus*pow(ca_cyt,4)));
	}

	this->m_initiated = true;
}

template<typename TDomain>
void RyR2<TDomain>::calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const
{
	// cast to side_type
	side_t* elem = dynamic_cast<side_t*>(e);
	if (!elem) UG_THROW("RyR2 fluxDensityFunction called with the wrong type of element.");

	number caCyt = u[_CCYT_];	// cytosolic Ca2+ concentration
	number caER = u[_CER_];		// ER Ca2+ concentration

	// membrane current corresponding to diffusion pressure
	number current = R*T/(4*F*F) * MU_RYR/REF_CA_ER * (caER - caCyt);

	// open probability
	number pOpen = m_aaO1[elem] + m_aaO2[elem];
	//number pOpen =  (1.0 + KB*pow(caCyt,3)) / (1.0 + KC + 1.0/(KA*pow(caCyt,4)) + KB*pow(caCyt,3));

	flux[0] = pOpen * current;
}


template<typename TDomain>
void RyR2<TDomain>::calc_flux_deriv(const std::vector<number>& u, GridObject* e, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const
{
	// we do not need to calculate any derivatives
	// as the whole channel dynamic is only considered in an explicit fashion
	size_t i = 0;
	if (!has_constant_value(_CCYT_))
	{
		flux_derivs[0][i].first = local_fct_index(_CCYT_);
		flux_derivs[0][i].second = 0.0;
		i++;
	}
	if (!has_constant_value(_CER_))
	{
		flux_derivs[0][i].first = local_fct_index(_CER_);
		flux_derivs[0][i].second = 0.0;
		i++;
	}
}


template<typename TDomain>
size_t RyR2<TDomain>::n_dependencies() const
{
	size_t n = 2;
	if (has_constant_value(_CCYT_))
		n--;
	if (has_constant_value(_CER_))
		n--;

	return n;
}


template<typename TDomain>
size_t RyR2<TDomain>::n_fluxes() const
{
	return 1;
};


template<typename TDomain>
const std::pair<size_t,size_t> RyR2<TDomain>::flux_from_to(size_t flux_i) const
{
    size_t from, to;
    if (allows_flux(_CCYT_)) to = local_fct_index(_CCYT_); else to = InnerBoundaryConstants::_IGNORE_;
    if (allows_flux(_CER_)) from = local_fct_index(_CER_); else from = InnerBoundaryConstants::_IGNORE_;

    return std::pair<size_t, size_t>(from, to);
}


template<typename TDomain>
const std::string RyR2<TDomain>::name() const
{
	return std::string("RyR2");
};


template<typename TDomain>
void RyR2<TDomain>::check_supplied_functions() const
{
	// Check that not both, inner and outer calcium concentrations are not supplied;
	// in that case, calculation of a flux would be of no consequence.
	if (!allows_flux(_CCYT_) && !allows_flux(_CER_))
	{
		UG_THROW("Supplying neither cytosolic nor endoplasmic calcium concentrations is not allowed.\n"
				"This would mean that the flux calculation would be of no consequence\n"
				"and this channel would not do anything.");
	}
}


template<typename TDomain>
void RyR2<TDomain>::print_units() const
{
	std::string nm = name();
	size_t n = nm.size();
	UG_LOG(std::endl);
	UG_LOG("+------------------------------------------------------------------------------+"<< std::endl);
	UG_LOG("|  Units used in the implementation of " << nm << std::string(n>=40?0:40-n, ' ') << "|" << std::endl);
	UG_LOG("|------------------------------------------------------------------------------|"<< std::endl);
	UG_LOG("|    Input                                                                     |"<< std::endl);
	UG_LOG("|      [Ca_cyt]  mM (= mol/m^3)                                                |"<< std::endl);
	UG_LOG("|      [Ca_er]   mM (= mol/m^3)                                                |"<< std::endl);
	UG_LOG("|                                                                              |"<< std::endl);
	UG_LOG("|    Output                                                                    |"<< std::endl);
	UG_LOG("|      Ca flux   mol/s                                                         |"<< std::endl);
	UG_LOG("+------------------------------------------------------------------------------+"<< std::endl);
	UG_LOG(std::endl);
}


// explicit template specializations
#ifdef UG_DIM_1
	template class RyR2<Domain1d>;
#endif
#ifdef UG_DIM_2
	template class RyR2<Domain2d>;
#endif
#ifdef UG_DIM_3
	template class RyR2<Domain3d>;
#endif


} // namespace neuro_collection
} // namespace ug



