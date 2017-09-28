/*
 * ca_wave_util_impl.h
 *
 *  Created on: 2017-09-18
 *      Author: mbreit
 */

#include "ca_wave_util.h"

#include "common/error.h"                       // UG_THROW etc.
#include "lib_disc/common/multi_index.h"        // DoFIndex

#include <sstream>
#include <vector>


namespace ug {
namespace neuro_collection {

template <typename TGridFunction>
number maxRyRFluxDensity
(
	ConstSmartPtr<TGridFunction> u,
	const char* fctNames,
	const char* subsetNames,
	ConstSmartPtr<RyRImplicit<typename TGridFunction::domain_type> > ryr
)
{
	typedef RyRImplicit<typename TGridFunction::domain_type> RyRType;

	// get (surface) DoF distro
	ConstSmartPtr<DoFDistribution> dd = u->dof_distribution();

	// get function indices for c_c, c_e, c1, c2 via FunctionGroup
	FunctionGroup fctGrp;
	try {fctGrp = u->fct_grp_by_name(fctNames);}
	UG_CATCH_THROW("Something went wrong during creation of function group for \n"
		"function names '" << fctNames << "' and provided solution.");
	UG_COND_THROW(fctGrp.size() != 4, "Function group must have exactly 5 entries.\n"
		"Make sure you provided exactly 4 function names in the following order:\n"
		"ca_cyt, ca_er, c1, c2.");

	size_t fi_cc = fctGrp[0];
	size_t fi_ce = fctGrp[1];
	size_t fi_c1 = fctGrp[2];
	size_t fi_c2 = fctGrp[3];

	// get subset indices via SubsetGroup
	SubsetGroup ssGrp;
	try {ssGrp = u->subset_grp_by_name(subsetNames);}
	UG_CATCH_THROW("Something went wrong during creation of subset group for \n"
		"membrane subset names '" << subsetNames << "' and provided solution.");

	// loop membrane subsets
	double currMax = 0.0;
	size_t nsi = ssGrp.size();
	for (size_t s = 0; s < nsi; ++s)
	{
		int si = ssGrp[s];

		// loop vertices on that subset
		DoFDistribution::traits<Vertex>::const_iterator it, itEnd;
		it = dd->begin<Vertex>(si);
		itEnd = dd->end<Vertex>(si);
		for (; it != itEnd; ++it)
		{
			Vertex* vrt = *it;

			// get function values in that vertex
			std::vector<number> values(5, 0.0);

			std::vector<DoFIndex> ind;
			UG_COND_THROW(!dd->is_def_in_subset(fi_cc, si),
				"Function " << fi_cc << " is not defined on subset " << si << ".");
			size_t numInd = dd->dof_indices(vrt, fi_cc, ind, true, true);
			UG_ASSERT(numInd == 1, "More (or less) than one function index found on a vertex!");
			values[RyRType::_CCYT_] = ryr->scale_input(RyRType::_CCYT_) * DoFRef(*u, ind[0]);

			UG_COND_THROW(!dd->is_def_in_subset(fi_ce, si),
				"Function " << fi_ce << " is not defined on subset " << si << ".");
			numInd = dd->dof_indices(vrt, fi_ce, ind, true, true);
			UG_ASSERT(numInd == 1, "More (or less) than one function index found on a vertex!");
			values[RyRType::_CER_] = ryr->scale_input(RyRType::_CER_) * DoFRef(*u, ind[0]);

			UG_COND_THROW(!dd->is_def_in_subset(fi_c1, si),
				"Function " << fi_c1 << " is not defined on subset " << si << ".");
			numInd = dd->dof_indices(vrt, fi_c1, ind, true, true);
			UG_ASSERT(numInd == 1, "More (or less) than one function index found on a vertex!");
			values[RyRType::_C1_] = ryr->scale_input(RyRType::_C1_) * DoFRef(*u, ind[0]);

			UG_COND_THROW(!dd->is_def_in_subset(fi_c2, si),
				"Function " << fi_c2 << " is not defined on subset " << si << ".");
			numInd = dd->dof_indices(vrt, fi_c2, ind, true, true);
			UG_ASSERT(numInd == 1, "More (or less) than one function index found on a vertex!");
			values[RyRType::_C2_] = ryr->scale_input(RyRType::_C2_) * DoFRef(*u, ind[0]);

			// calculate flux density
			std::vector<number> flux(1, 0.0);
			ryr->calc_flux(values, vrt, flux);

			// compare to current max
			currMax = std::max(currMax, (double) fabs(flux[0]));
		}
	}

	// max over all processes
#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		double local = currMax;
		com.allreduce(&local, &currMax, 1, PCL_DT_DOUBLE, PCL_RO_MAX);
	}
#endif

	return (number) currMax;
}


template <typename TGridFunction>
number waveFrontX
(
	ConstSmartPtr<TGridFunction> u,
	const char* fctNames,
	const char* subsetNames,
	number thresh
)
{
	// get coordinate attachment accessor
	const typename TGridFunction::domain_type::position_accessor_type&
		aaPos = u->domain()->position_accessor();

	// get (surface) DoF distro
	ConstSmartPtr<DoFDistribution> dd = u->dof_distribution();

	// get function indices for c1, c2 via FunctionGroup
	FunctionGroup fctGrp;
	try {fctGrp = u->fct_grp_by_name(fctNames);}
	UG_CATCH_THROW("Something went wrong during creation of function group for \n"
		"function names '" << fctNames << "' and provided solution.");
	UG_COND_THROW(fctGrp.size() != 2, "Function group must have exactly 5 entries.\n"
		"Make sure you provided exactly 2 function names in the following order:\n"
		"c1, c2.");

	size_t fi_c1 = fctGrp[0];
	size_t fi_c2 = fctGrp[1];

	// get subset indices via SubsetGroup
	SubsetGroup ssGrp;
	try {ssGrp = u->subset_grp_by_name(subsetNames);}
	UG_CATCH_THROW("Something went wrong during creation of subset group for \n"
		"membrane subset names '" << subsetNames << "' and provided solution.");

	// loop membrane subsets
	double xMax = -std::numeric_limits<double>::max();
	size_t nsi = ssGrp.size();
	for (size_t s = 0; s < nsi; ++s)
	{
		int si = ssGrp[s];

		// loop vertices on that subset
		DoFDistribution::traits<Vertex>::const_iterator it, itEnd;
		it = dd->begin<Vertex>(si);
		itEnd = dd->end<Vertex>(si);
		for (; it != itEnd; ++it)
		{
			Vertex* vrt = *it;

			// only treat vertices to the right of the current max
			number xCoord = aaPos[vrt][0];
			if ((double) xCoord <= xMax)
				continue;

			// get function values in that vertex
			std::vector<DoFIndex> ind;

			UG_COND_THROW(!dd->is_def_in_subset(fi_c1, si),
				"Function " << fi_c1 << " is not defined on subset " << si << ".");
			size_t numInd = dd->dof_indices(vrt, fi_c1, ind, true, true);
			UG_ASSERT(numInd == 1, "More (or less) than one function index found on a vertex!");
			number c1 = DoFRef(*u, ind[0]);

			UG_COND_THROW(!dd->is_def_in_subset(fi_c2, si),
				"Function " << fi_c2 << " is not defined on subset " << si << ".");
			numInd = dd->dof_indices(vrt, fi_c2, ind, true, true);
			UG_ASSERT(numInd == 1, "More (or less) than one function index found on a vertex!");
			number c2 = DoFRef(*u, ind[0]);

			// calculate open probability
			number pOpen = 1.0 - (c1 + c2);

			// update xMax if necessary
			if (pOpen > thresh)
				xMax = xCoord;
		}
	}

	// max over all processes
#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		double local = xMax;
		com.allreduce(&local, &xMax, 1, PCL_DT_DOUBLE, PCL_RO_MAX);
	}
#endif

	return xMax;
}



} // namespace ug
} // namespace neuro_collection
