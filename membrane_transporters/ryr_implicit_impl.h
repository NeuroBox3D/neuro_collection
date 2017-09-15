/*
 * ryr_implicit_impl.h
 * Fully implicit RyR implementation.
 *
 * Date:   2017-08-16
 * Author: mbreit
 */

#include "lib_grid/algorithms/debug_util.h"                 // for ElementDebugInfo


namespace ug {
namespace neuro_collection {

template <typename TDomain>
template <typename TVector>
void RyRImplicit<TDomain>::calculate_steady_state(SmartPtr<TVector> u) const
{
	typedef typename DoFDistribution::traits<side_t>::const_iterator it_type;
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


	// loop dof distro vertices of ER membrane subsets
	std::vector<DoFIndex> dofIndexCa;
	std::vector<DoFIndex> dofIndexO2;
	std::vector<DoFIndex> dofIndexC1;
	std::vector<DoFIndex> dofIndexC2;

	SubsetGroup ssg(this->approx_space()->domain()->subset_handler(), this->m_vSubset);
	size_t si_sz = ssg.size();
	for (size_t si = 0; si < si_sz; ++si)
	{
		it_type it = dd->begin<side_t>(ssg[si]);
		it_type it_end = dd->end<side_t>(ssg[si]);
		for (; it != it_end; ++it)
		{
			// get DoFs for all involved unknowns
			dd->dof_indices(*it, ind_o2, dofIndexO2, true, true);
			dd->dof_indices(*it, ind_c1, dofIndexC1, true, true);
			dd->dof_indices(*it, ind_c2, dofIndexC2, true, true);

			number ca_cyt = 0.0;
			if (!this->has_constant_value(_CCYT_, ca_cyt))
			{
				dd->dof_indices(*it, ind_ccyt, dofIndexCa, true, true);
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





template <typename TDomain>
template <typename TVector>
void RyRImplicit_1drotsym<TDomain>::calculate_steady_state(SmartPtr<TVector> u) const
{
	typedef typename DoFDistribution::traits<Edge>::const_iterator it_type;
	ConstSmartPtr<ApproximationSpace<TDomain> > approxSpace = this->approx_space();
	UG_COND_THROW(!approxSpace.valid(), "Approximation space not present in implicit RyR discretization."
		<< std::endl << " Did you forget to add it to the domain discretization?");
	ConstSmartPtr<DoFDistribution> dd = approxSpace->dof_distribution(GridLevel(), false);

	// get global fct index for our functions
	FunctionGroup fctGrp(dd->dof_distribution_info());
	fctGrp.add(IElemDisc<TDomain>::m_vFct);
	const size_t ind_ccyt = fctGrp.unique_id(_CCYT_);
	const size_t ind_o2 = fctGrp.unique_id(_O2_);
	const size_t ind_c1 = fctGrp.unique_id(_C1_);
	const size_t ind_c2 = fctGrp.unique_id(_C2_);


	// loop dof distro vertices of ER membrane subsets
	std::vector<DoFIndex> dofIndexCa;
	std::vector<DoFIndex> dofIndexO2;
	std::vector<DoFIndex> dofIndexC1;
	std::vector<DoFIndex> dofIndexC2;

	SubsetGroup ssg(this->approx_space()->domain()->subset_handler(), this->m_vSubset);
	size_t si_sz = ssg.size();
	for (size_t si = 0; si < si_sz; ++si)
	{
		it_type it = dd->begin<Edge>(ssg[si]);
		it_type it_end = dd->end<Edge>(ssg[si]);
		for (; it != it_end; ++it)
		{
			// get DoFs for all involved unknowns (in vertices) -- ONLY P1 Lagrange allowed!
			dd->dof_indices(*it, ind_ccyt, dofIndexCa, true, true);
			UG_ASSERT(dofIndexCa.size() == 2, "Not exactly 2 DoFs found for function " << ind_ccyt
				<< " in edge " << ElementDebugInfo(*this->approx_space()->domain()->grid(), *it));

			dd->dof_indices(*it, ind_o2, dofIndexO2, true, true);
			UG_ASSERT(dofIndexO2.size() == 2, "Not exactly 2 DoFs found for function " << ind_o2
				<< " in edge " << ElementDebugInfo(*this->approx_space()->domain()->grid(), *it));

			dd->dof_indices(*it, ind_c1, dofIndexC1, true, true);
			UG_ASSERT(dofIndexC1.size() == 2, "Not exactly 2 DoFs found for function " << ind_c1
				<< " in edge " << ElementDebugInfo(*this->approx_space()->domain()->grid(), *it));

			dd->dof_indices(*it, ind_c2, dofIndexC2, true, true);
			UG_ASSERT(dofIndexC2.size() == 2, "Not exactly 2 DoFs found for function " << ind_c2
				<< " in edge " << ElementDebugInfo(*this->approx_space()->domain()->grid(), *it));

			// write equilibrium state probs to functions
			for (size_t i = 0; i < 2; ++i)
			{
				number ca_cyt = DoFRef(*u, dofIndexCa[i]) * m_scale_cc;

				// calculate equilibrium
				number KA = KAplus/KAminus * ca_cyt*ca_cyt*ca_cyt*ca_cyt;
				number KB = KBplus/KBminus * ca_cyt*ca_cyt*ca_cyt;
				number KC = KCplus/KCminus;

				number denom_inv = 1.0 / (1.0 + KC + 1.0/KA + KB);

				DoFRef(*u, dofIndexO2[i]) = KB * denom_inv;
				DoFRef(*u, dofIndexC1[i]) = denom_inv / KA;
				DoFRef(*u, dofIndexC2[i]) = KC * denom_inv;
			}
		}
	}
}


} // namespace neuro_collection
} // namespace ug



