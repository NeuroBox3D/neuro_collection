/*
 * misc_util.cpp
 *
 *  Created on: 2018-01-30
 *      Author: mbreit
 */

#include "misc_util.h"

#include "lib_disc/domain.h"                                        // for Domain1d, Doma...
#include "lib_disc/domain_traits.h"                                 // for domain_traits
#include "lib_disc/function_spaces/approximation_space.h"           // for ApproximationSpace
#include "lib_grid/algorithms/geom_obj_util/anisotropy_util.h"      // for is_anisotropic
#include "lib_grid/tools/grid_level.h"                              // for GridLevel


namespace ug {
namespace neuro_collection {


template <typename TDomain>
void mark_global(SmartPtr<IRefiner> refiner, SmartPtr<ApproximationSpace<TDomain> > approx)
{
	typedef typename domain_traits<TDomain::dim>::element_type elem_type;
	typedef typename DoFDistribution::traits<elem_type>::const_iterator const_iterator;

	// get surface dof distribution
	ConstSmartPtr<DoFDistribution> dd = approx->dof_distribution(GridLevel(), false);

	const_iterator iter = dd->template begin<elem_type>();
	const_iterator iterEnd = dd->template end<elem_type>();

	// loop elements for marking
	for (; iter != iterEnd; ++iter)
		refiner->mark(*iter, RM_REFINE);
}



template <typename TDomain>
void mark_anisotropic
(
	SmartPtr<IRefiner> refiner,
	SmartPtr<ApproximationSpace<TDomain> > approx,
	number thresholdRatio
)
{
	typedef typename domain_traits<TDomain::dim>::element_type elem_type;
	typedef typename domain_traits<TDomain::dim>::side_type side_type;
	typedef typename DoFDistribution::traits<elem_type>::const_iterator const_iterator;

	Grid& grid = *refiner->grid();
	typename TDomain::position_accessor_type aaPos = approx->domain()->position_accessor();

	// get surface dof distribution and prepare loop over all surface elements
	ConstSmartPtr<DoFDistribution> dd = approx->dof_distribution(GridLevel(), false);
	const_iterator iter = dd->template begin<elem_type>();
	const_iterator iterEnd = dd->template end<elem_type>();

	// loop elements for marking
	std::vector<Edge*> longEdges;
	for (; iter != iterEnd; ++iter)
	{
		longEdges.clear();
		AnisotropyState state = long_edges_of_anisotropic_elem(*iter, grid, aaPos, thresholdRatio, longEdges);
		if (state != ISOTROPIC)
		{
			// mark elem
			refiner->mark(*iter, RM_CLOSURE);

			// mark long edges
			const size_t nEdges = longEdges.size();
			for (size_t e = 0; e < nEdges; ++e)
				refiner->mark(longEdges[e], RM_FULL);

			// mark all sides
			typename Grid::traits<side_type>::secure_container sl;
			grid.associated_elements(sl, *iter);
			const size_t slSz = sl.size();
			for (size_t s = 0; s < slSz; ++s)
				if (refiner->get_mark(sl[s]) != RM_FULL)
					refiner->mark(sl[s], RM_CLOSURE);
		}
	}
}


template <typename TDomain>
void mark_anisotropic_onlyX
(
	SmartPtr<IRefiner> refiner,
	SmartPtr<ApproximationSpace<TDomain> > approx,
	number thresholdRatio
)
{
	typedef typename domain_traits<TDomain::dim>::element_type elem_type;
	typedef typename domain_traits<TDomain::dim>::side_type side_type;
	typedef typename DoFDistribution::traits<elem_type>::const_iterator const_iterator;

	Grid& grid = *refiner->grid();
	typename TDomain::position_accessor_type aaPos = approx->domain()->position_accessor();

	// get surface dof distribution and prepare loop over all surface elements
	ConstSmartPtr<DoFDistribution> dd = approx->dof_distribution(GridLevel(), false);
	const_iterator iter = dd->template begin<elem_type>();
	const_iterator iterEnd = dd->template end<elem_type>();

	// loop elements for marking
	std::vector<Edge*> longEdges;
	for (; iter != iterEnd; ++iter)
	{
		longEdges.clear();
		AnisotropyState state = long_edges_of_anisotropic_elem(*iter, grid, aaPos, thresholdRatio, longEdges);
		if (state == ISOTROPIC)
			continue;

		UG_COND_THROW(!longEdges.size(), "Element is anisotropic, but no long edges present.");

		// check whether edges point in x direction
		Edge* longEdge = longEdges[0];
		MathVector<TDomain::dim> dir;
		VecSubtract(dir, aaPos[longEdge->vertex(1)], aaPos[longEdge->vertex(0)]);
		VecNormalize(dir, dir);
		if (fabs(dir[0]) > 0.9)
		{
			// mark elem
			refiner->mark(*iter, RM_CLOSURE);

			// mark long edges
			const size_t nEdges = longEdges.size();
			for (size_t e = 0; e < nEdges; ++e)
				refiner->mark(longEdges[e], RM_FULL);

			// mark all sides
			typename Grid::traits<side_type>::secure_container sl;
			grid.associated_elements(sl, *iter);
			const size_t slSz = sl.size();
			for (size_t s = 0; s < slSz; ++s)
				if (refiner->get_mark(sl[s]) != RM_FULL)
					refiner->mark(sl[s], RM_CLOSURE);
		}
	}
}



// explicit template specializations
#ifdef UG_DIM_1
	template void mark_global<Domain1d>(SmartPtr<IRefiner>, SmartPtr<ApproximationSpace<Domain1d> >);
	template void mark_anisotropic<Domain1d>(SmartPtr<IRefiner>, SmartPtr<ApproximationSpace<Domain1d> >, number);
	template void mark_anisotropic_onlyX<Domain1d>(SmartPtr<IRefiner>, SmartPtr<ApproximationSpace<Domain1d> >, number);
#endif
#ifdef UG_DIM_2
	template void mark_global<Domain2d>(SmartPtr<IRefiner>, SmartPtr<ApproximationSpace<Domain2d> >);
	template void mark_anisotropic<Domain2d>(SmartPtr<IRefiner>, SmartPtr<ApproximationSpace<Domain2d> >, number);
	template void mark_anisotropic_onlyX<Domain2d>(SmartPtr<IRefiner>, SmartPtr<ApproximationSpace<Domain2d> >, number);
#endif
#ifdef UG_DIM_3
	template void mark_global<Domain3d>(SmartPtr<IRefiner>, SmartPtr<ApproximationSpace<Domain3d> >);
	template void mark_anisotropic<Domain3d>(SmartPtr<IRefiner>, SmartPtr<ApproximationSpace<Domain3d> >, number);
	template void mark_anisotropic_onlyX<Domain3d>(SmartPtr<IRefiner>, SmartPtr<ApproximationSpace<Domain3d> >, number);
#endif


}  // namespace neuro_collection
}  // namespace ug
