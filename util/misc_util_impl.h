/*
 * misc_util_impl.h
 *
 *  Created on: 2018-01-30
 *      Author: mbreit
 */

#include <cstddef>                                                  // for size_t, NULL

#include "common/error.h"                                           // for UG_CATCH_THROW, UG_THROW, UG_CO...

#include "lib_disc/common/multi_index.h"                            // for DoFIndex, DoFRef
#include "lib_disc/dof_manager/dof_distribution.h"                  // for DoFDistribution, DoFDistributio...
#include "lib_disc/domain_traits.h"                                 // for domain_traits
#include "lib_disc/function_spaces/approximation_space.h"           // for ApproximationSpace
#include "lib_grid/tools/grid_level.h"                              // for GridLevel
#include "lib_grid/tools/surface_view.h"                            // for SurfaceView::ConstSurfaceViewEl...


namespace ug {
namespace neuro_collection {



template <typename TBaseElem, typename TGridFunction>
static void scale_dof_indices
(
	ConstSmartPtr<DoFDistribution> dd,
	SmartPtr<TGridFunction> vecOut,
	ConstSmartPtr<TGridFunction> vecIn,
	const std::vector<number>& vScale
)
{
	typename DoFDistribution::traits<TBaseElem>::const_iterator iter, iterEnd;
	std::vector<DoFIndex> vInd;

	try
	{
		// iterate all elements (including SHADOW_RIM_COPY!)
		iter = dd->template begin<TBaseElem>(SurfaceView::ALL);
		iterEnd = dd->template end<TBaseElem>(SurfaceView::ALL);
		for (; iter != iterEnd; ++iter)
		{
			for (size_t fi = 0; fi < dd->num_fct(); ++fi)
			{
				size_t nInd = dd->inner_dof_indices(*iter, fi, vInd);

				// remember multi indices
				for (size_t dof = 0; dof < nInd; ++dof)
					DoFRef(*vecOut, vInd[dof]) = vScale[fi] * DoFRef(*vecIn, vInd[dof]);
			}
		}
	}
	UG_CATCH_THROW("Error while scaling vector.")
}


template <typename TGridFunction>
void scale_dimless_vector
(
	SmartPtr<TGridFunction> scaledVecOut,
	ConstSmartPtr<TGridFunction> dimlessVecIn,
	const std::vector<number>& scalingFactors
)
{
	// check that the correct numbers of scaling factors are given
	size_t n = scalingFactors.size();
	UG_COND_THROW(n != dimlessVecIn->num_fct(), "Number of scaling factors (" << n << ") "
			"does not match number of functions given in dimless vector (" << dimlessVecIn->num_fct() << ").");

	// check that input and output vectors have the same number of components and dofs
	UG_COND_THROW(n != scaledVecOut->num_fct(), "Input and output vectors do not have "
			"the same number of functions (" << n << " vs. " << scaledVecOut->num_fct() << ").");
	for (size_t fct = 0; fct < n; ++fct)
	{
		UG_COND_THROW(dimlessVecIn->num_dofs(fct) != scaledVecOut->num_dofs(fct),
				"Input and output vectors do not have the same number of DoFs for function " << fct
				<< " (" << dimlessVecIn->num_dofs(fct) << " vs. " << scaledVecOut->num_dofs(fct) << ").");
	}

	ConstSmartPtr<DoFDistribution> dd = dimlessVecIn->dof_distribution();

	if (dd->max_dofs(VERTEX))
		scale_dof_indices<Vertex, TGridFunction>(dd, scaledVecOut, dimlessVecIn, scalingFactors);
	if (dd->max_dofs(EDGE))
		scale_dof_indices<Edge, TGridFunction>(dd, scaledVecOut, dimlessVecIn, scalingFactors);
	if (dd->max_dofs(FACE))
		scale_dof_indices<Face, TGridFunction>(dd, scaledVecOut, dimlessVecIn, scalingFactors);
	if (dd->max_dofs(VOLUME))
		scale_dof_indices<Volume, TGridFunction>(dd, scaledVecOut, dimlessVecIn, scalingFactors);
}




template <typename TDomain>
void mark_global(SmartPtr<IRefiner> refiner, SmartPtr<ApproximationSpace<TDomain> > approx)
{
	typedef typename domain_traits<TDomain::dim>::element_type elem_type;
	typedef typename DoFDistribution::traits<elem_type>::const_iterator const_iterator;

	// get surface dof distribution
	ConstSmartPtr<DoFDistribution> dd = approx->dof_distribution(GridLevel(), false);

	const_iterator iter = dd->template begin<elem_type>();
	const_iterator iterEnd = dd->template end<elem_type>();

//	loop elements for marking
	for (; iter != iterEnd; ++iter)
		refiner->mark(*iter, RM_REFINE);
}



} // namspace neuro_collection
} // namespace ug

