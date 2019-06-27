/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2019-02-06
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

#include "ryr_discrete.h"

#include "common/error.h"  // for UG_THROW
#include "common/util/string_util.h"  // for TokenizeTrimString
#include "lib_disc/common/function_group.h"  // for FunctionGroup
#include "lib_disc/common/multi_index.h"  // for DoFRef
#include "lib_disc/spatial_disc/ass_tuner.h"  // for CT_CONSTRAINTS
#include "lib_grid/algorithms/debug_util.h"  // for ElementDebugInfo
#include "lib_grid/tools/subset_group.h"  // for SubsetGroup



namespace ug {
namespace neuro_collection {



template <typename TDomain, typename TAlgebra>
RyRDiscrete<TDomain, TAlgebra>::RyRDiscrete(const char* functions, const char* subsets)
: R(8.314), T(310.0), F(96485.0),
  KAplus(1500.0e12), KBplus(1500.0e9), KCplus(1.75),
  KAminus(28.8), KBminus(385.9), KCminus(0.1),
  MU_RYR(5.0e-11), REF_CA_ER(2.5e-1),
  m_vSubsetNames(TokenizeTrimString(subsets)),
  m_vFunctionNames(TokenizeTrimString(functions)),
  m_cutoffOpenProb(0.002)
{}


template <typename TDomain, typename TAlgebra>
RyRDiscrete<TDomain, TAlgebra>::
RyRDiscrete(const std::vector<std::string>& functions, const std::vector<std::string>& subsets)
: R(8.314), T(310.0), F(96485.0),
  KAplus(1500.0e12), KBplus(1500.0e9), KCplus(1.75),
  KAminus(28.8), KBminus(385.9), KCminus(0.1),
  MU_RYR(5.0e-11), REF_CA_ER(2.5e-1),
  m_vSubsetNames(subsets),
  m_vFunctionNames(functions),
  m_cutoffOpenProb(0.002)
{}


template <typename TDomain, typename TAlgebra>
RyRDiscrete<TDomain, TAlgebra>::~RyRDiscrete()
{}



template <typename TDomain, typename TAlgebra>
void RyRDiscrete<TDomain, TAlgebra>::adjust_defect
(
	vector_type& d,
	const vector_type& u,
	ConstSmartPtr<DoFDistribution> dd,
	int type,
	number time,
	ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
	const std::vector<number>* vScaleMass,
	const std::vector<number>* vScaleStiff
)
{
	std::vector<DoFIndex> vDI(5);

#ifdef UG_PARALLEL
	const DistributedGridManager& dgm = *dd->multi_grid()->distributed_grid_manager();
#endif

	// get number of stages of time stepping scheme (or 1 in stationary case)
	size_t nTimes = 1;
	if (vSol.valid() && vScaleMass && vScaleStiff)
	{
		nTimes = vSol->size();

		// it may happen (LIMEX) that the sizes of vScaleMass and vScaleStiff are less than that of vSol
		// in that case simply forget about additional solutions given
		if (vScaleMass->size() < nTimes || vScaleStiff->size() < nTimes)
			nTimes = std::min(vScaleMass->size(), vScaleStiff->size());
	}

	// loop subsets
	typedef typename DoFDistribution::traits<Vertex>::const_iterator vrt_it_type;
	const size_t nSI = m_vSI.size();
	const size_t nFct = m_vFctMap.size();

	for (size_t i = 0; i < nSI; ++i)
	{
		int si = m_vSI[i];

		// loop all vertices in the subset
		vrt_it_type it = dd->begin<Vertex>(si);
		vrt_it_type itEnd = dd->end<Vertex>(si);
		for (; it != itEnd; ++it)
		{
			Vertex* vrt = *it;

#ifdef UG_PARALLEL
			// do not process horizontal slaves in the parallel case to prevent double counting
			if (dgm.get_status(vrt) & ES_H_SLAVE)
				continue;
#endif

			// get all DoF indices that concern us here
			vDI.clear();
			for (size_t j = 0; j < nFct; ++j)
			{
				size_t fct = m_vFctMap[j];

				// we do not need hanging vertices as all transport mechanism
				// vertices must be in the coarse grid
				dd->dof_indices(vrt, fct, vDI, false, false);
			}
			UG_COND_THROW(vDI.size() != 5, "Not exactly 5 DoF indices for vertex.");


			// add stiffness defect for cytosolic calcium / ER calcium, gating params
			for (size_t k = 0; k < nTimes; ++k)
			{
				if (vScaleStiff && !(*vScaleStiff)[k])
					continue;

				number caCyt;
				number caER;
				number o2;
				number c1;
				number c2;

				if (vSol.valid())
				{
					caCyt = 1e3 * DoFRef(*vSol->solution(k), vDI[_CCYT_]);   // scale from M to mM
					caER = 1e3 * DoFRef(*vSol->solution(k), vDI[_CER_]);   // scale from M to mM
					o2 = DoFRef(*vSol->solution(k), vDI[_O2_]);
					c1 = DoFRef(*vSol->solution(k), vDI[_C1_]);
					c2 = DoFRef(*vSol->solution(k), vDI[_C2_]);
				}
				else
				{
					caCyt = 1e3 * DoFRef(u, vDI[_CCYT_]);   // scale from M to mM
					caER = 1e3 * DoFRef(u, vDI[_CER_]);   // scale from M to mM
					o2 = DoFRef(u, vDI[_O2_]);
					c1 = DoFRef(u, vDI[_C1_]);
					c2 = DoFRef(u, vDI[_C2_]);
				}
				const number o1 = 1.0 - (o2 + c1 + c2);

				number pOpen = 1.0 - (c1 + c2);

				// really close channels below cutoff (0.002 corresponds to equil. cy_cyt of ~8e-5mM)
				// using a cubic spline between x0=cutoff and x1=2*cutoff
				// satisfying f'(x0) = 0; f(x0) = 0, f'(x1) = 1, f(x1) = x1
				if (pOpen < 2*m_cutoffOpenProb)
				{
					if (pOpen > m_cutoffOpenProb)
					{
						const number x0 = m_cutoffOpenProb;
						const number x1 = 2*m_cutoffOpenProb;
						const number sa = - (x1-x0) / ((x1+x0)*(x1+x0)*(x1+x0));
						const number sb = 2.0*(x1*x1 + x1*x0 + x0*x0) / ((x1+x0)*(x1+x0)*(x1+x0));
						const number sc = -3.0*x0*x0*sa - 2*x0*sb;
						const number sd = -x0*x0*x0*sa - x0*x0*sb - x0*sc;
						pOpen = sa*pOpen*pOpen*pOpen + sb*pOpen*pOpen + sc*pOpen + sd;
					}
					else
						pOpen = 0.0;
				}

				number current = pOpen * R*T/(4*F*F) * MU_RYR/REF_CA_ER * (caER - caCyt);

				number dt = 1.0;
				if (vScaleStiff)
					dt = (*vScaleStiff)[k];

				DoFRef(d, vDI[_CCYT_]) -= current * dt * 1e15;  // scale from mol to mol (um/dm)^3
				DoFRef(d, vDI[_CER_]) += current * dt * 1e15;  // scale from mol to mol (um/dm)^3

				DoFRef(d, vDI[_O2_]) -= dt * (KBplus * caCyt*caCyt*caCyt * o1 - KBminus * o2);
				DoFRef(d, vDI[_C1_]) -= dt * (KAminus * o1 - KAplus * caCyt*caCyt*caCyt*caCyt * c1);
				DoFRef(d, vDI[_C2_]) -= dt * (KCplus * o1 - KCminus * c2);
			}

			// add mass defects for gating parameters
			if (vSol.valid() && vScaleMass)
			{
				for (size_t k = 0; k < nTimes; ++k)
				{
					DoFRef(d, vDI[_O2_]) += (*vScaleMass)[k] * DoFRef(*vSol->solution(k), vDI[_O2_]);
					DoFRef(d, vDI[_C1_]) += (*vScaleMass)[k] * DoFRef(*vSol->solution(k), vDI[_C1_]);
					DoFRef(d, vDI[_C2_]) += (*vScaleMass)[k] * DoFRef(*vSol->solution(k), vDI[_C2_]);
				}
			}
		}
	}
}


template <typename TDomain, typename TAlgebra>
void RyRDiscrete<TDomain, TAlgebra>::adjust_jacobian
(
	matrix_type& J,
	const vector_type& u,
	ConstSmartPtr<DoFDistribution> dd,
	int type,
	number time,
	ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
	const number s_a0
)
{
	std::vector<DoFIndex> vDI(5);

#ifdef UG_PARALLEL
	const DistributedGridManager& dgm = *dd->multi_grid()->distributed_grid_manager();
#endif

	// loop subsets
	typedef typename DoFDistribution::traits<Vertex>::const_iterator vrt_it_type;
	const size_t nSI = m_vSI.size();
	const size_t nFct = m_vFctMap.size();

	for (size_t i = 0; i < nSI; ++i)
	{
		int si = m_vSI[i];

		// loop all vertices in the subset
		vrt_it_type it = dd->begin<Vertex>(si);
		vrt_it_type itEnd = dd->end<Vertex>(si);
		for (; it != itEnd; ++it)
		{
			Vertex* vrt = *it;

#ifdef UG_PARALLEL
			// do not process horizontal slaves in the parallel case to prevent double counting
			if (dgm.get_status(vrt) & ES_H_SLAVE)
				continue;
#endif

			// get all DoF indices that concern us here
			vDI.clear();
			for (size_t j = 0; j < nFct; ++j)
			{
				size_t fct = m_vFctMap[j];

				// we do not need hanging vertices as all transport mechanism
				// vertices must be in the coarse grid
				dd->dof_indices(vrt, fct, vDI, false, false);
			}
			UG_COND_THROW(vDI.size() != 5, "Not exactly 5 DoF indices for vertex.");


			// add stiffness entries for cytosolic calcium / ER calcium, gating params
			if (s_a0)
			{
				const number caCyt = 1e3 * DoFRef(u, vDI[_CCYT_]);   // scale from M to mM
				const number caER = 1e3 * DoFRef(u, vDI[_CER_]);   // scale from M to mM
				const number o2 = DoFRef(u, vDI[_O2_]);
				const number c1 = DoFRef(u, vDI[_C1_]);
				const number c2 = DoFRef(u, vDI[_C2_]);
				const number o1 = 1.0 - (o2 + c1 + c2);

				number pOpen = 1.0 - (c1 + c2);
				number dPOdC12 = -1.0;

				// really close channels below cutoff (0.002 corresponds to equil. cy_cyt of ~8e-5mM)
				// using a cubic spline between x0=cutoff and x1=2*cutoff
				// satisfying f'(x0) = 0; f(x0) = 0, f'(x1) = 1, f(x1) = x1
				if (pOpen < 2*m_cutoffOpenProb)
				{
					if (pOpen > m_cutoffOpenProb)
					{
						const number x0 = m_cutoffOpenProb;
						const number x1 = 2*m_cutoffOpenProb;
						const number sa = - (x1-x0) / ((x1+x0)*(x1+x0)*(x1+x0));
						const number sb = 2.0*(x1*x1 + x1*x0 + x0*x0) / ((x1+x0)*(x1+x0)*(x1+x0));
						const number sc = -3.0*x0*x0*sa - 2*x0*sb;
						const number sd = -x0*x0*x0*sa - x0*x0*sb - x0*sc;

						dPOdC12 = -(3.0*sa*pOpen*pOpen + 2*sb*pOpen + sc);
						pOpen = sa*pOpen*pOpen*pOpen + sb*pOpen*pOpen + sc*pOpen + sd;
					}
					else
					{
						pOpen = 0.0;
						dPOdC12 = 0.0;
					}
				}

				DoFRef(J, vDI[_CCYT_], vDI[_CCYT_]) += pOpen * R*T/(4*F*F) * MU_RYR/REF_CA_ER * s_a0 * 1e18;
				DoFRef(J, vDI[_CCYT_], vDI[_CER_]) -= pOpen * R*T/(4*F*F) * MU_RYR/REF_CA_ER * s_a0 * 1e18;
				DoFRef(J, vDI[_CCYT_], vDI[_C1_]) -= dPOdC12 * R*T/(4*F*F) * MU_RYR/REF_CA_ER * (caER - caCyt) * s_a0 * 1e15;
				DoFRef(J, vDI[_CCYT_], vDI[_C2_]) -= dPOdC12 * R*T/(4*F*F) * MU_RYR/REF_CA_ER * (caER - caCyt) * s_a0 * 1e15;

				DoFRef(J, vDI[_CER_], vDI[_CCYT_]) -= pOpen * R*T/(4*F*F) * MU_RYR/REF_CA_ER * s_a0 * 1e18;
				DoFRef(J, vDI[_CER_], vDI[_CER_]) += pOpen * R*T/(4*F*F) * MU_RYR/REF_CA_ER * s_a0 * 1e18;
				DoFRef(J, vDI[_CER_], vDI[_C1_]) += dPOdC12 * R*T/(4*F*F) * MU_RYR/REF_CA_ER * (caER - caCyt) * s_a0 * 1e15;
				DoFRef(J, vDI[_CER_], vDI[_C2_]) += dPOdC12 * R*T/(4*F*F) * MU_RYR/REF_CA_ER * (caER - caCyt) * s_a0 * 1e15;

				DoFRef(J, vDI[_O2_], vDI[_CCYT_]) -= s_a0 * KBplus * 3.0*caCyt*caCyt * o1 * 1e3;
				DoFRef(J, vDI[_O2_], vDI[_O2_]) += s_a0 * (KBminus + KBplus * caCyt*caCyt*caCyt);
				DoFRef(J, vDI[_O2_], vDI[_C1_]) += s_a0 * KBplus * caCyt*caCyt*caCyt;
				DoFRef(J, vDI[_O2_], vDI[_C2_]) += s_a0 * KBplus * caCyt*caCyt*caCyt;

				DoFRef(J, vDI[_C1_], vDI[_CCYT_]) += s_a0 * KAplus * 4.0*caCyt*caCyt*caCyt * c1 * 1e3;
				DoFRef(J, vDI[_C1_], vDI[_O2_]) += s_a0 * KAminus;
				DoFRef(J, vDI[_C1_], vDI[_C1_]) += s_a0 * (KAminus + KAplus * caCyt*caCyt*caCyt*caCyt);
				DoFRef(J, vDI[_C1_], vDI[_C2_]) += s_a0 * KAminus;

				DoFRef(J, vDI[_C2_], vDI[_O2_]) += s_a0 * KCplus;
				DoFRef(J, vDI[_C2_], vDI[_C1_]) += s_a0 * KCplus;
				DoFRef(J, vDI[_C2_], vDI[_C2_]) += s_a0 * (KCminus + KCplus);
			}

			// add mass defects for gating parameters
			DoFRef(J, vDI[_O2_], vDI[_O2_]) += 1.0;
			DoFRef(J, vDI[_C1_], vDI[_C1_]) += 1.0;
			DoFRef(J, vDI[_C2_], vDI[_C2_]) += 1.0;
		}
	}
}


template <typename TDomain, typename TAlgebra>
void RyRDiscrete<TDomain, TAlgebra>::adjust_solution
(
	vector_type& u,
	ConstSmartPtr<DoFDistribution> dd,
	int type,
	number time
)
{
	// nothing to do
}


template <typename TDomain, typename TAlgebra>
void RyRDiscrete<TDomain, TAlgebra>::adjust_linear
(
	matrix_type& mat,
	vector_type& rhs,
	ConstSmartPtr<DoFDistribution> dd,
	int type,
	number time
)
{
	UG_THROW("RyRDiscrete::adjust_linear() not implemented.");
}


template <typename TDomain, typename TAlgebra>
void RyRDiscrete<TDomain, TAlgebra>::adjust_rhs
(
	vector_type& rhs,
	const vector_type& u,
	ConstSmartPtr<DoFDistribution> dd,
	int type,
	number time
)
{
	UG_THROW("RyRDiscrete::adjust_rhs() not implemented.");
}


template <typename TDomain, typename TAlgebra>
int RyRDiscrete<TDomain, TAlgebra>::type() const
{
	return CT_CONSTRAINTS;
}


template <typename TDomain, typename TAlgebra>
void RyRDiscrete<TDomain, TAlgebra>::
set_approximation_space(SmartPtr<ApproximationSpace<TDomain> > approxSpace)
{
	// call base class method
	IDomainConstraint<TDomain, TAlgebra>::set_approximation_space(approxSpace);

	// convert subset names to indices
	SubsetGroup sg(this->m_spApproxSpace->domain()->subset_handler());
	const size_t nSs = m_vSubsetNames.size();
	m_vSI.resize(nSs);
	try {sg.add(m_vSubsetNames);}
	UG_CATCH_THROW("Subsets could not be added to subset group.");
	for (size_t i = 0; i < nSs; ++i)
		m_vSI[i] = sg[i];

	// convert function names to indices
	FunctionGroup fg(this->m_spApproxSpace->function_pattern());
	const size_t nFct = m_vFunctionNames.size();
	UG_COND_THROW(nFct != 5, "RyRDiscrete needs exactly 5 functions in constructor.");
	m_vFctMap.resize(nFct);
	try {fg.add(m_vFunctionNames);}
	UG_CATCH_THROW("Functions could not be added to function group.");
	for (size_t i = 0; i < nFct; ++i)
		m_vFctMap[i] = fg[i];
}


template <typename TDomain, typename TAlgebra>
void RyRDiscrete<TDomain, TAlgebra>::calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const
{
	const number caCyt = u[_CCYT_];   // scale from M to mM
	const number caER = u[_CER_];   // scale from M to mM
	const number c1 = u[_C1_];
	const number c2 = u[_C2_];
	number pOpen = 1.0 - (c1 + c2);

	// really close channels below cutoff (0.002 corresponds to equil. cy_cyt of ~8e-5mM)
	// using a cubic spline between x0=cutoff and x1=2*cutoff
	// satisfying f'(x0) = 0; f(x0) = 0, f'(x1) = 1, f(x1) = x1
	if (pOpen < 2*m_cutoffOpenProb)
	{
		if (pOpen > m_cutoffOpenProb)
		{
			const number x0 = m_cutoffOpenProb;
			const number x1 = 2*m_cutoffOpenProb;
			const number sa = -1.0 / ((x1-x0)*(x1+x0));
			const number sb = 2.0 / (x1-x0);
			const number sc = -3.0*x0*x0*sa - 2*x0*sb;
			const number sd = -x0*x0*x0*sa - x0*x0*sb - x0*sc;
			pOpen = sa*pOpen*pOpen*pOpen + sb*pOpen*pOpen + sc*pOpen + sd;
		}
		else
			pOpen = 0.0;
	}

	flux[0] = pOpen * R*T/(4*F*F) * MU_RYR/REF_CA_ER * (caER - caCyt);
}


template <typename TDomain, typename TAlgebra>
number RyRDiscrete<TDomain, TAlgebra>::scale_input(const size_t i) const
{
	switch (i)
	{
		case _CCYT_:
		case _CER_: return 1e3;
		case _O2_:
		case _C1_:
		case _C2_: return 1.0;
		default: return 1.0;
	}
}


template <typename TDomain, typename TAlgebra>
void RyRDiscrete<TDomain, TAlgebra>::calculate_steady_state(SmartPtr<vector_type> u) const
{
	typedef typename DoFDistribution::traits<Vertex>::const_iterator it_type;
	UG_COND_THROW(!this->m_spApproxSpace.valid(), "Approximation space not present in discrete RyR discretization."
		<< std::endl << " Did you forget to add it to the domain discretization?");
	ConstSmartPtr<DoFDistribution> dd = this->m_spApproxSpace->dof_distribution(GridLevel(), false);

	// get global fct index for our functions
	FunctionGroup fctGrp(dd->dof_distribution_info());
	fctGrp.add(m_vFunctionNames);
	const size_t ind_ccyt = fctGrp.unique_id(_CCYT_);
	const size_t ind_o2 = fctGrp.unique_id(_O2_);
	const size_t ind_c1 = fctGrp.unique_id(_C1_);
	const size_t ind_c2 = fctGrp.unique_id(_C2_);


	// loop dof distro vertices of ER membrane subsets
	std::vector<DoFIndex> dofIndexCa;
	std::vector<DoFIndex> dofIndexO2;
	std::vector<DoFIndex> dofIndexC1;
	std::vector<DoFIndex> dofIndexC2;

	SubsetGroup ssg(this->m_spApproxSpace->domain()->subset_handler(), m_vSubsetNames);
	size_t si_sz = ssg.size();
	for (size_t si = 0; si < si_sz; ++si)
	{
		it_type it = dd->begin<Vertex>(ssg[si]);
		it_type it_end = dd->end<Vertex>(ssg[si]);
		for (; it != it_end; ++it)
		{
			// get DoFs for all involved unknowns
			dd->dof_indices(*it, ind_o2, dofIndexO2, true, true);
			dd->dof_indices(*it, ind_c1, dofIndexC1, true, true);
			dd->dof_indices(*it, ind_c2, dofIndexC2, true, true);
			dd->dof_indices(*it, ind_ccyt, dofIndexCa, true, true);

			UG_COND_THROW(dofIndexO2.size() != 1,
				"Not exactly 1 DoF on vertex for o2 (" << dofIndexCa.size() << ") for "
				<< ElementDebugInfo(*this->m_spApproxSpace->domain()->grid(), *it));
			UG_COND_THROW(dofIndexCa.size() != dofIndexO2.size(),
				"Not the same number of DoFs on the same element for cytosolic calcium ("
				<< dofIndexCa.size() << ") and channel state O2 (" << dofIndexO2.size() << ") "
				"for " << ElementDebugInfo(*this->m_spApproxSpace->domain()->grid(), *it));
			UG_COND_THROW(dofIndexC1.size() != dofIndexO2.size(),
				"Not the same number of DoFs on the same element for channel state C1 ("
				<< dofIndexC1.size() << ") and channel state O2 (" << dofIndexO2.size() << ") "
				"for " << ElementDebugInfo(*this->m_spApproxSpace->domain()->grid(), *it));
			UG_COND_THROW(dofIndexC2.size() != dofIndexO2.size(),
				"Not the same number of DoFs on the same element for channel state C2 ("
				<< dofIndexC2.size() << ") and channel state O2 (" << dofIndexO2.size() << ") "
				"for " << ElementDebugInfo(*this->m_spApproxSpace->domain()->grid(), *it));

			// write equilibrium state probs to functions
			number ca_cyt = 1e3 * DoFRef(*u, dofIndexCa[0]);  // scale from M to mM

			// calculate equilibrium
			number KA = KAplus/KAminus * ca_cyt*ca_cyt*ca_cyt*ca_cyt;
			number KB = KBplus/KBminus * ca_cyt*ca_cyt*ca_cyt;
			number KC = KCplus/KCminus;

			number denom_inv = 1.0 / (1.0 + KC + 1.0/KA + KB);

			DoFRef(*u, dofIndexO2[0]) = KB * denom_inv;
			DoFRef(*u, dofIndexC1[0]) = denom_inv / KA;
			DoFRef(*u, dofIndexC2[0]) = KC * denom_inv;
		}
	}
}


template <typename TDomain, typename TAlgebra>
void RyRDiscrete<TDomain, TAlgebra>::set_cutoff_open_probability(number cutoffProb)
{
	m_cutoffOpenProb = cutoffProb;
}


// explicit template specializations
#ifdef UG_CPU_1
	#ifdef UG_DIM_1
		template class RyRDiscrete<Domain1d, CPUAlgebra>;
	#endif
	#ifdef UG_DIM_2
		template class RyRDiscrete<Domain2d, CPUAlgebra>;
	#endif
	#ifdef UG_DIM_3
		template class RyRDiscrete<Domain3d, CPUAlgebra>;
	#endif
#endif

} // end namespace neuro_collection
} // end namespace ug

