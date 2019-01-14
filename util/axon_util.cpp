/*
 * axon_util.cpp
 *
 *  Created on: 2018-10-29
 *      Author: mbreit
 */

#include "axon_util.h"

#include <cstddef>                                                  // for size_t, NULL
#include <stack>

#include "common/error.h"                                           // for UG_CATCH_THROW, UG_THROW, UG_CO...
#include "common/types.h"                                           // for number
#include "lib_disc/common/geometry_util.h"                          // for SideNormal
#include "lib_disc/dof_manager/dof_distribution.h"                  // for DoFDistribution, DoFDistributio...
#include "lib_disc/domain.h"                                        // for Domain2d, Domain3d...
#include "lib_disc/domain_traits.h"                                 // for domain_traits
#include "lib_disc/domain_util.h"                                   // for CollectCornerCoordinates
#include "lib_disc/function_spaces/approximation_space.h"           // for ApproximationSpace
#include "lib_grid/grid/neighborhood_util.h"                        // for GetConnectedNeighbors
#include "lib_grid/grid_objects/grid_dim_traits.h"                  // for grid_dim_traits
#include "lib_grid/tools/grid_level.h"                              // for GridLevel
#include "lib_grid/tools/subset_group.h"                            // for SubsetGroup

namespace ug {
namespace neuro_collection {


template <typename TDomain>
void unmark_ranvier_areas
(
	SmartPtr<IRefiner> refiner,
	SmartPtr<ApproximationSpace<TDomain> > approx,
	const std::string& ranvierSubsets,
	bool doUnmark
)
{
	typedef typename domain_traits<TDomain::dim>::element_type elem_type;
	typedef typename domain_traits<TDomain::dim>::side_type side_type;
	typedef typename DoFDistribution::traits<elem_type>::const_iterator const_elem_iterator;
	typedef typename DoFDistribution::traits<side_type>::const_iterator const_side_iterator;

	// get surface dof distribution
	const typename TDomain::position_accessor_type aaPos = approx->domain()->position_accessor();
	ConstSmartPtr<DoFDistribution> dd = approx->dof_distribution(GridLevel(), false);

	// get grid
	Grid& grid = *approx->domain()->grid();

	SubsetGroup sg(approx->domain()->subset_handler());
	try {sg.add(TokenizeString(ranvierSubsets));}
	UG_CATCH_THROW("Ranvier subsets could not be identified in SubsetHandler.");

	std::vector<MathVector<TDomain::dim> > corners;

	// calculate the ranvier centers
	const size_t nSs = sg.size();
	std::vector<number> ranvierCentersX(nSs);
	std::vector<int> ranvierCenterFound(nSs, 0);
	for (size_t ss = 0; ss < nSs; ++ss)
	{
		ranvierCentersX[ss] = 0.0;

		// iterate ranvier elements
		const int ssDim = sg.dim(ss);
		if (ssDim == TDomain::dim)
		{
			const_elem_iterator iter = dd->template begin<elem_type>(sg[ss]);
			const_elem_iterator iterEnd = dd->template end<elem_type>(sg[ss]);
			for (; iter != iterEnd; ++iter)
			{
				ranvierCentersX[ss] += CalculateGridObjectCenter(*iter, aaPos)[0];
				++ranvierCenterFound[ss];
			}
		}
		else if (ssDim == TDomain::dim - 1)
		{
			const_side_iterator iter = dd->template begin<side_type>(sg[ss]);
			const_side_iterator iterEnd = dd->template end<side_type>(sg[ss]);
			for (; iter != iterEnd; ++iter)
			{
				ranvierCentersX[ss] += CalculateGridObjectCenter(*iter, aaPos)[0];
				++ranvierCenterFound[ss];
			}
		}
		else
			UG_THROW("Ranvier subsets must either be full-dim elements or sides.")
	}

	// communicate ranvier centers
#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator pc;
		std::vector<number> ranvierCentersXSend(ranvierCentersX);
		pc.allreduce(ranvierCentersXSend, ranvierCentersX, PCL_RO_SUM);
		std::vector<int> ranvierCenterFoundSend(ranvierCenterFound);
		pc.allreduce(ranvierCenterFoundSend, ranvierCenterFound, PCL_RO_SUM);
	}
#endif

	for (size_t ss = 0; ss < nSs; ++ss)
	{
		UG_COND_THROW(ranvierCenterFound[ss] == 0,
			"Ranvier node '" << sg.name(ss) << "' not found on any proc.");

		ranvierCentersX[ss] = ranvierCentersX[ss] / ranvierCenterFound[ss];
	}


	MathVector<TDomain::dim> normal;

	const_elem_iterator iter = dd->template begin<elem_type>();
	const_elem_iterator iterEnd = dd->template end<elem_type>();
	for (; iter != iterEnd; ++iter)
	{
		elem_type* elem = *iter;

		// calculate object center
		number objectCenter = CalculateGridObjectCenter(elem, aaPos)[0];
		number objectSize = VecDistance(aaPos[elem->vertex(0)], aaPos[elem->vertex(1)]);

		// only work on this elem if center is at one of the ranvier coords
		bool isRanvierElem = false;
		for (size_t ss = 0; ss < nSs; ++ss)
		{
			if (fabs(objectCenter - ranvierCentersX[ss]) < 1e-3 * objectSize)
			{
				isRanvierElem = true;
				break;
			}
		}
		if (!isRanvierElem)
			continue;

		// unmark element
		if (doUnmark)
			refiner->mark(elem, RM_NONE);

		// loop sides to remove their marks as well
		ReferenceObjectID roid = elem->reference_object_id();
		CollectCornerCoordinates(corners, *elem, *approx->domain(), true);

		// FIXME: This is not enough: In parallel cases, it might happen that
		//        a ranvier area element has neighbors, but on other procs.
		// we will mark the elem RM_CLOSURE if all neighbors will be refined (and at least one exists)
		bool allAxialNeighborsWillBeRefined = true;
		int nAxialNeighbors = 0;

		typename Grid::traits<side_type>::secure_container sl;
		grid.associated_elements_sorted(sl, elem);
		const size_t slSz = sl.size();
		for (size_t s = 0; s < slSz; ++s)
		{
			side_type* side = sl[s];
			if (doUnmark)
				refiner->mark(side, RM_NONE);

			// if we have an axial side, check whether neighbor will be refined
			// if not: do not refine this element either
			SideNormal<TDomain::dim>(roid, normal, s, &corners[0]);
			VecNormalize(normal, normal);
			if (fabs(normal[0]) > 0.05)
			{
				elem_type* nb = GetConnectedNeighbor(grid, side, elem);
				if (nb)
				{
					++nAxialNeighbors;
					if (
						refiner->get_mark(nb) == RM_NONE
#ifdef UG_PARALLEL
					&& !grid.distributed_grid_manager()->contains_status(nb, ES_V_MASTER)
#endif
					)
					{
						allAxialNeighborsWillBeRefined = false;
					}
				}
			}
		}

		if (nAxialNeighbors && allAxialNeighborsWillBeRefined)
			refiner->mark(elem, RM_CLOSURE);
	}
}


// explicit template instantiations
#ifdef UG_DIM_1
	template void unmark_ranvier_areas<Domain1d>(SmartPtr<IRefiner>, SmartPtr<ApproximationSpace<Domain1d> >, const std::string&, bool);
#endif
#ifdef UG_DIM_2
	template void unmark_ranvier_areas<Domain2d>(SmartPtr<IRefiner>, SmartPtr<ApproximationSpace<Domain2d> >, const std::string&, bool);
#endif
#ifdef UG_DIM_3
	template void unmark_ranvier_areas<Domain3d>(SmartPtr<IRefiner>, SmartPtr<ApproximationSpace<Domain3d> >, const std::string&, bool);
#endif


} // namspace neuro_collection
} // namespace ug

