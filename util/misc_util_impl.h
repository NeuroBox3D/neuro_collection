/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2019-07-05
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

#include "lib_disc/common/multi_index.h"  // for DoFIndex
#include "lib_disc/domain_traits.h"  // for domain_traits


namespace ug {
namespace neuro_collection {


// //////////////////////////////////////////////////////////////////////////////
// 	Check values are within bounds
// //////////////////////////////////////////////////////////////////////////////
template <typename TGridFunction, typename TBaseElem>
unsigned long MarkOutOfRangeElems
(
	SmartPtr<IRefiner> refiner,
	ConstSmartPtr<TGridFunction> u,
	size_t cmp,
	number lowerBnd,
	number upperBnd
)
{
	typedef typename TGridFunction::template traits<TBaseElem>::const_iterator elem_it;
	typedef typename domain_traits<TGridFunction::domain_type::dim>::element_type elem_type;
	typedef typename Grid::traits<elem_type>::secure_container elem_list;

	Grid& grid = *refiner->grid();

	unsigned long nMarked = 0;
	elem_list el;

	// loop the elements in the subset
	std::vector<DoFIndex> vDI;
	elem_it it = u->template begin<TBaseElem>();
	elem_it itEnd = u->template end<TBaseElem>();
	for (; it != itEnd; ++it)
	{
		TBaseElem* elem = *it;

		// loop indices at this element
		const size_t nInd = u->inner_dof_indices(elem, cmp, vDI, true);
		for (size_t i = 0; i < nInd; ++i)
		{
			const number& val = DoFRef(*u, vDI[i]);
			if (val < lowerBnd || val > upperBnd)
			{
				// mark neighbors for refinement
				grid.associated_elements(el, elem);
				const size_t elSz = el.size();
				for (size_t e = 0; e < elSz; ++e)
				{
					if (refiner->get_mark(el[e]) == RM_NONE)
					{
						refiner->mark(el[e], RM_FULL);
						++nMarked;
					}
				}
			}
		}
	}

	return nMarked;
}

template <typename TGridFunction>
void MarkOutOfRangeElems
(
	SmartPtr<IRefiner> refiner,
	ConstSmartPtr<TGridFunction> u,
	size_t cmp,
	number lowerBnd,
	number upperBnd
)
{
	unsigned long nMarked = 0;
	if (u->max_fct_dofs(cmp, 0))
		nMarked += MarkOutOfRangeElems<TGridFunction, Vertex>(refiner, u, cmp, lowerBnd, upperBnd);
	if (u->max_fct_dofs(cmp, 1))
		nMarked += MarkOutOfRangeElems<TGridFunction, Edge>(refiner, u, cmp, lowerBnd, upperBnd);
	if (u->max_fct_dofs(cmp, 2))
		nMarked += MarkOutOfRangeElems<TGridFunction, Face>(refiner, u, cmp, lowerBnd, upperBnd);
	if (u->max_fct_dofs(cmp, 3))
		nMarked += MarkOutOfRangeElems<TGridFunction, Volume>(refiner, u, cmp, lowerBnd, upperBnd);

#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator pc;
		nMarked = pc.allreduce(nMarked, PCL_RO_SUM);
	}
#endif


	if (nMarked)
		UG_LOGN("  +++ Marked for refinement: " << nMarked << " elements.");
}

#if 0
template <typename TGridFunction>
void MarkOutOfRangeElems
(
	SmartPtr<IRefiner> refiner,
	ConstSmartPtr<TGridFunction> u,
	const char* fctNames,
	number lowerBnd,
	number upperBnd
)
{
	std::vector<std::string> vFctNames;
	TokenizeTrimString (fctNames, vFctNames);
	FunctionGroup fctGroup (u->function_pattern());
	const size_t nFct = vFctNames.size();
	for (size_t f = 0; f < nFct; ++f)
	{
		try {fctGroup.add(vFctNames[f]);}
		UG_CATCH_THROW("Could not add function " << vFctNames[f] << " to function group.");
	}

	const size_t fctGrpSz = fctGroup.size();
	for (size_t f = 0; f < fctGrpSz; ++f)
		MarkOutOfRangeElems(refiner, u, fctGroup[f], lowerBnd, upperBnd);
}
#endif

} // namespace ug
} // namespace neuro_collection

