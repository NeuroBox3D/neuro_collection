/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2017-08-16
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

#include "lib_grid/algorithms/debug_util.h"                 // for ElementDebugInfo


namespace ug {
namespace neuro_collection {

template <typename TDomain>
template <typename TGridFunction>
void RyRLinearized<TDomain>::calculate_steady_state(SmartPtr<TGridFunction> u) const
{
	typedef typename DoFDistribution::traits<Vertex>::const_iterator it_type;
	ConstSmartPtr<ApproximationSpace<TDomain> > approxSpace = this->approx_space();
	UG_COND_THROW(!approxSpace.valid(), "Approximation space not present in implicit RyR discretization."
		<< std::endl << " Did you forget to add it to the domain discretization?");
	ConstSmartPtr<DoFDistribution> dd = approxSpace->dof_distribution(GridLevel(), false);

	// get global fct index for our functions
	FunctionGroup fctGrp(dd->dof_distribution_info());
	fctGrp.add(IMembraneTransporter::m_vFct);
	const size_t ind_ccyt = fctGrp.unique_id(_CCYT_);
	const size_t ind_o2 = fctGrp.unique_id(_O2_);
	const size_t ind_c1 = fctGrp.unique_id(_C1_);
	const size_t ind_c2 = fctGrp.unique_id(_C2_);

	const bool bUseP1Interpolation =
		u->local_finite_element_id(ind_o2).type() == LFEID::LAGRANGE
		&& u->local_finite_element_id(ind_o2).order() == 1
		&& u->local_finite_element_id(ind_c1).type() == LFEID::LAGRANGE
		&& u->local_finite_element_id(ind_c1).order() == 1
		&& u->local_finite_element_id(ind_c2).type() == LFEID::LAGRANGE
		&& u->local_finite_element_id(ind_c2).order() == 1;
	UG_COND_THROW(!bUseP1Interpolation, "Only implemented for P1 functions at the moment, "
		"but either o2, c1 or c2 is not.");

	// loop dof distro vertices of ER membrane subsets
	std::vector<DoFIndex> dofIndexCa;
	std::vector<DoFIndex> dofIndexO2;
	std::vector<DoFIndex> dofIndexC1;
	std::vector<DoFIndex> dofIndexC2;

	SubsetGroup ssg(this->approx_space()->domain()->subset_handler(), this->m_vSubset);
	size_t si_sz = ssg.size();
	for (size_t si = 0; si < si_sz; ++si)
	{
		it_type it = dd->begin<Vertex>(ssg[si]);
		it_type it_end = dd->end<Vertex>(ssg[si]);
		for (; it != it_end; ++it)
		{
			// get DoFs for all involved unknowns
			dd->inner_dof_indices(*it, ind_o2, dofIndexO2, true);
			dd->inner_dof_indices(*it, ind_c1, dofIndexC1, true);
			dd->inner_dof_indices(*it, ind_c2, dofIndexC2, true);

			number ca_cyt = 0.0;
			if (!this->has_constant_value(_CCYT_, ca_cyt))
			{
				dd->inner_dof_indices(*it, ind_ccyt, dofIndexCa, true);
				UG_COND_THROW(dofIndexCa.size() != dofIndexO2.size(),
					"Not the same number of DoFs on the same element for cytosolic calcium ("
					<< dofIndexCa.size() << ") and channel state O2 (" << dofIndexO2.size() << ") "
					"for " << ElementDebugInfo(*this->approx_space()->domain()->grid(), *it));
			} // else the constant value has been written to ca_cyt by has_constant_value()

			UG_COND_THROW(dofIndexC1.size() != dofIndexO2.size(),
				"Not the same number of DoFs on the same element for channel state C1 ("
				<< dofIndexC1.size() << ") and channel state O2 (" << dofIndexO2.size() << ") "
				"for " << ElementDebugInfo(*this->approx_space()->domain()->grid(), *it));

			UG_COND_THROW(dofIndexC2.size() != dofIndexO2.size(),
				"Not the same number of DoFs on the same element for channel state C2 ("
				<< dofIndexC2.size() << ") and channel state O2 (" << dofIndexO2.size() << ") "
				"for " << ElementDebugInfo(*this->approx_space()->domain()->grid(), *it));

			// write equilibrium state probs to functions
			std::size_t sz = dofIndexO2.size();
			for (size_t i = 0; i < sz; ++i)
			{
				if (!this->has_constant_value(_CCYT_))
					ca_cyt = DoFRef(*u, dofIndexCa[i]);

				ca_cyt *= this->scale_input(_CCYT_);

				// calculate equilibrium
				number KA = KAplus/KAminus * ca_cyt*ca_cyt*ca_cyt*ca_cyt;
				number KB = KBplus/KBminus * ca_cyt*ca_cyt*ca_cyt;
				number KC = KCplus/KCminus;

				number denom_inv = 1.0 / (1.0 + KC + 1.0/KA + KB);

				DoFRef(*u, dofIndexO2[i]) = KB * denom_inv / this->scale_input(_O2_);
				DoFRef(*u, dofIndexC1[i]) = denom_inv / KA / this->scale_input(_C1_);
				DoFRef(*u, dofIndexC2[i]) = KC * denom_inv / this->scale_input(_C2_);
			}
		}
	}
}



} // namespace neuro_collection
} // namespace ug



