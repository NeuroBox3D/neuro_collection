/*
 * vdcc_bg_impl.h
 *
 *  Created on: 05.02.2013
 *      Author: mbreit
 */

#include "vdcc_bg.h"

namespace ug{
namespace neuro_collection{


template<typename TDomain>
template<int TType>
void VDCC_BG<TDomain>::set_channel_type()
{
	if (m_initiated)
	{
		UG_THROW("Borg-Graham channel type can not be changed after initialization.\n"
				 "Call set_channel_type() BEFORE init().");
	}

	m_channelType = TType;
	check_supplied_functions();

	switch (TType)
	{
		// - all gating params according to table 5, page 98 of the Borg-Graham article
		//   "Interpretations of data and mechanisms for hippocampal pyramidal cell models"
		// - conductance value lambda of 5.5 pS from Gollasch et al. - "Single calcium
		//   channel currents of arterial smooth muscle at physiological calcium concentrations"
		//	 and used to derive permeability values p_T, p_N and p_L inspired by
		//   Fisher et al. - "Properties and distribution of single
		//   voltage-gated calcium channels in adult hippocampal neurons" in an offset way:
		//
		//	 p_T = 1.9e-19 m^3/s; p_N = 2*p_T; p_L = 3*p_T
		//
		//	 where p = lambda * R * T / (valence^2 * F^2 * Ca_o)
		//
		case BG_Ntype:
			m_gpMGate = GatingParams(3.4, -21.0, 1.5);
			m_gpHGate = GatingParams(-2.0, -40.0, 75.0);
			m_perm = 3.8e-19;
			m_mp = 2;
			m_hp = 1;
			break;
		case BG_Ltype:
			m_gpMGate = GatingParams(4.6, -1.0, 1.5);
			m_gpHGate = GatingParams(0.0, 0.0, 0.0);
			m_perm = 5.7e-19;
			m_mp = 2;
			m_hp = 0;
			break;
		case BG_Ttype:
			m_gpMGate = GatingParams(3.0, -36.0, 1.5);
			m_gpHGate = GatingParams(-5.2, -68.0, 10.0);
			m_perm = 1.9e-19;
			m_mp = 2;
			m_hp = 1;
			break;
		default:
			UG_THROW("Type of Borg-Graham channel does not match any of the pre-implemented.");
			break;
	}
}



template <typename TDomain>
template <typename TVector>
void VDCC_BG<TDomain>::calculate_steady_state(SmartPtr<TVector> u, number vm) const
{
	typedef typename DoFDistribution::traits<side_t>::const_iterator it_type;
	ConstSmartPtr<ApproximationSpace<TDomain> > approxSpace = this->approx_space();
	UG_COND_THROW(!approxSpace.valid(), "Approximation space not present in Borg-Graham VDCC discretization."
		<< std::endl << " Did you forget to add it to the domain discretization?");
	ConstSmartPtr<DoFDistribution> dd = approxSpace->dof_distribution(GridLevel(), false);

	// get global fct index for our functions
	FunctionGroup fctGrp(dd->dof_distribution_info());
	fctGrp.add(IMembraneTransporter::m_vFct);
	const size_t ind_m = fctGrp.unique_id(_M_ - m_localIndicesOffset);
	size_t ind_h = 0;
	if (has_hGate())
		ind_h = fctGrp.unique_id(_H_ - m_localIndicesOffset);

	// loop dof distro vertices of plasma membrane subsets
	std::vector<DoFIndex> dofIndexM;
	std::vector<DoFIndex> dofIndexH;

	SubsetGroup ssg(approxSpace->domain()->subset_handler(), this->m_vSubset);
	size_t si_sz = ssg.size();
	for (size_t si = 0; si < si_sz; ++si)
	{
		it_type it = dd->begin<side_t>(ssg[si]);
		it_type it_end = dd->end<side_t>(ssg[si]);
		for (; it != it_end; ++it)
		{
			// get DoFs for all involved unknowns
			dd->dof_indices(*it, ind_m, dofIndexM, true, true);
			if (has_hGate())
				dd->dof_indices(*it, ind_h, dofIndexH, true, true);

			UG_COND_THROW(has_hGate() && dofIndexM.size() != dofIndexH.size(),
				"Not the same number of DoFs on the same element for state variables m ("
				<< dofIndexM.size() << ") and h (" << dofIndexH.size() << ") "
				"for " << ElementDebugInfo(*approxSpace->domain()->grid(), *it));

			// write equilibrium state probs to functions
			const size_t sz = dofIndexM.size();
			for (size_t i = 0; i < sz; ++i)
			{
				DoFRef(*u, dofIndexM[i]) = calc_gating_start(m_gpMGate, 1e3*vm);
				if (has_hGate())
					DoFRef(*u, dofIndexH[i]) = calc_gating_start(m_gpHGate, 1e3*vm);
			}
		}
	}
}


} // neuro_collection
} // namespace ug
