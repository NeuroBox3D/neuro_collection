/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2013-02-05
 *
 * This file is part of NeuroBox, which is based on UG4.
 *
 * NeuroBox and UG4 are free software: You can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3
 * (as published by the Free Software Foundation) with the following additional
 * attribution requirements (according to LGPL/GPL v3 §7):
 *
 * (1) The following notice must be displayed in the appropriate legal notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 *
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 *
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating PDE based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * "Stepniewski, M., Breit, M., Hoffer, M. and Queisser, G.
 *   NeuroBox: computational mathematics in multiscale neuroscience.
 *   Computing and visualization in science (2019).
 * "Breit, M. et al. Anatomically detailed and large-scale simulations studying
 *   synapse loss and synchrony using NeuroBox. Front. Neuroanat. 10 (2016), 8"
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
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

	const bool had_hGate = has_hGate();

	m_channelType = TType;

	// check availability of gating functions
	if (!m_bUseGatingAttachments)
	{
		UG_COND_THROW(has_hGate() && !is_supplied(_H_),
			"You provided the m gating parameter as function, but not the h gating parameter.\n"
			"This is not allowed. Either provide both or none.");
	}
	check_supplied_functions();

	if (m_bUseGatingAttachments)
	{
		if (has_hGate() && !had_hGate)
		{
			if (!m_mg->template has_attachment<vm_grid_object>(this->m_HGate))
				m_mg->template attach_to<vm_grid_object>(this->m_HGate);
			m_aaHGate = Grid::AttachmentAccessor<vm_grid_object, ADouble>(*m_mg, m_HGate);
		}
		else if (!has_hGate() && had_hGate)
		{
			if (m_mg->template has_attachment<vm_grid_object>(this->m_HGate))
				m_mg->detach_from<vm_grid_object>(this->m_HGate);
		}
	}

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
template <typename TGridFunction>
void VDCC_BG<TDomain>::calculate_steady_state(SmartPtr<TGridFunction> u, number vm) const
{
	typedef typename DoFDistribution::traits<Vertex>::const_iterator it_type;
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

	const bool bUseP1Interpolation = (u->local_finite_element_id(ind_m).type() == LFEID::LAGRANGE
		&& u->local_finite_element_id(ind_m).order() == 1)
		&& (!has_hGate() || (u->local_finite_element_id(ind_h).type() == LFEID::LAGRANGE
			&& u->local_finite_element_id(ind_h).order() == 1));
	UG_COND_THROW(!bUseP1Interpolation, "Only implemented for P1 functions at the moment.");

	// loop dof distro vertices of plasma membrane subsets
	std::vector<DoFIndex> dofIndexM;
	std::vector<DoFIndex> dofIndexH;

	SubsetGroup ssg(approxSpace->domain()->subset_handler(), this->m_vSubset);
	size_t si_sz = ssg.size();
	for (size_t si = 0; si < si_sz; ++si)
	{
		it_type it = dd->begin<Vertex>(ssg[si]);
		it_type it_end = dd->end<Vertex>(ssg[si]);
		for (; it != it_end; ++it)
		{
			// get DoFs for all involved unknowns
			dd->inner_dof_indices(*it, ind_m, dofIndexM, true);
			if (has_hGate())
				dd->inner_dof_indices(*it, ind_h, dofIndexH, true);

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
