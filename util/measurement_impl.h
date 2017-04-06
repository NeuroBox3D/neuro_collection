/*
 * measurement.cpp
 *
 *  Created on: 31.05.2016
 *      Author: mbreit
 */

#include "measurement.h"

#include "common/error.h"	// UG_THROW etc.
#include "lib_disc/function_spaces/integrate.h"	// IntegrateSubset

#include <sstream>



namespace ug {
namespace neuro_collection {




// ////////////////////// //
// compute volume command //
// ////////////////////// //

template <typename TDomain, int dim>
void collectVol
(
	ConstSmartPtr<DoFDistribution> dofDistr,
	const typename TDomain::position_accessor_type& aaPos,
	size_t si, number& vol
)
{
	const int worldDim = TDomain::dim;

	// determine types
	typedef typename domain_traits<dim>::grid_base_object Elem;
	typedef typename DoFDistribution::template traits<Elem>::const_iterator iter_type;

	//	get element iterator for current subset
	iter_type iter = dofDistr->template begin<Elem>(si);
	iter_type iterEnd = dofDistr->template end<Elem>(si);

	// loop over all elements
	for (; iter != iterEnd; ++iter)
	{
		// get current element
		Elem* elem = *iter;

		ReferenceObjectID roid = elem->reference_object_id();

		// collect corner coords
		std::vector<typename TDomain::position_type> coco;
		CollectCornerCoordinates(coco, elem, aaPos, false);

		vol += ElementSize<worldDim>(roid, &coco[0]);
	}
}


template <typename TDomain>
number computeVolume(ConstSmartPtr<ApproximationSpace<TDomain> > approx, const int subset)
{
	const int worldDim = TDomain::dim;

	ConstSmartPtr<DoFDistribution> dofDistr = approx->dof_distribution(GridLevel());
	const typename TDomain::position_accessor_type& aaPos = approx->domain()->position_accessor();

	// initialize result variable
	number vol = 0.0;

	// further implementation depends on dim
	int dim = dofDistr->dim_subset(subset);

	// collect volume
	if (dim == worldDim)
		collectVol<TDomain, worldDim>(dofDistr, aaPos, subset, vol);
	else if (dim == worldDim-1 && worldDim >= 1)
		collectVol<TDomain, worldDim>=1 ? worldDim-1 : 1>(dofDistr, aaPos, subset, vol);
	else if (dim == worldDim-2 && worldDim >= 2)
		collectVol<TDomain, worldDim>=2 ? worldDim-2 : 1>(dofDistr, aaPos, subset, vol);
	else if (dim == worldDim-3 && worldDim >= 3)
		collectVol<TDomain, worldDim>=3 ? worldDim-3 : 1>(dofDistr, aaPos, subset, vol);
	else {UG_THROW("Unknown dim (" << dim << ") or worldDim (" << worldDim << ").");}

#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		// sum up volumes and integrals on all processes
		pcl::ProcessCommunicator com;
		number local = vol;
		com.allreduce(&local, &vol, 1, PCL_DT_DOUBLE, PCL_RO_SUM);
	}
#endif

	return vol;
}


template <typename TDomain>
void computeVolume
(
	ConstSmartPtr<ApproximationSpace<TDomain> > approx,
	const char* subsetNames
)
{
	// retrieve domain, dofDistr and dim from approxSpace
	ConstSmartPtr<TDomain> domain = approx->domain();
	ConstSmartPtr<DoFDistribution> dofDistr = approx->dof_distribution(GridLevel());

	// get subset group to be measured on
	SubsetGroup ssGrp;
	try {ssGrp = dofDistr->subset_grp_by_name(subsetNames);}
	UG_CATCH_THROW("At least one of the subsets in '" << subsetNames
					<< "' is not contained in the approximation space (or something else was wrong).");

	// loop subsets
	UG_LOG("\n");
	for (size_t si = 0; si < ssGrp.size(); si++)
	{
		// initialize measurement result variables
		number vol = computeVolume(approx, ssGrp[si]);

		// subset dim
		int dim = ssGrp.dim(si);

		if (dim == 3) {UG_LOG("volume of subset '" << ssGrp.name(si) << "':\t" << vol << "\n");}
		else if (dim == 2) {UG_LOG("area of subset '" << ssGrp.name(si) << "':\t" << vol << "\n");}
		else if (dim == 1) {UG_LOG("length of subset '" << ssGrp.name(si) << "':\t" << vol << "\n");}
		else if (dim == 0) {UG_LOG("Number of vertices in subset '" << ssGrp.name(si) << "':\t" << vol << "\n");}
		else {UG_THROW("Unknown dim " << dim << ".");}
	}
	UG_LOG("\n");

	return;
}


// //////////////////// //
// measurement commands //
// //////////////////// //
template <typename TGridFunction>
number takeMeasurement
(
	SmartPtr<TGridFunction> solution,
	const number time,
	const char* subsetNames,
	const char* functionNames,
	const char* outFileName
) {
	return takeMeasurement(solution, time, subsetNames, functionNames, outFileName, "");
}


template <typename TGridFunction>
number takeMeasurement
(
	SmartPtr<TGridFunction> solution,
	const number time,
	const char* subsetNames,
	const char* functionNames,
	const char* outFileName,
	const char* outFileExt
)
{
	typedef typename TGridFunction::domain_type domain_type;
	const int worldDim = domain_type::dim;

	// retrieve approximation space, domain, dofDistr and dim from grid function
	ConstSmartPtr<ApproximationSpace<domain_type> > approx = solution->approx_space();

	// get subset group to be measured on
	SubsetGroup ssGrp;
	try {ssGrp = solution->subset_grp_by_name(subsetNames);}
	UG_CATCH_THROW("At least one of the subsets in '" << subsetNames
					<< "' is not contained in the approximation space (or something else was wrong).");

	// get function group to be measured
	FunctionGroup fctGrp;
	try {fctGrp = solution->fct_grp_by_name(functionNames);}
	UG_CATCH_THROW("At least one of the functions in '" << functionNames
					<< "' is not contained in the approximation space (or something else was wrong).");

	number value = 0.;
	number vol = 0.;
	// loop subsets
	for (size_t si = 0; si < ssGrp.size(); si++)
	{
		vol = computeVolume(approx, ssGrp[si]);

		int dim = ssGrp.dim(si);

		// loop functions
		for (size_t fi = 0; fi < fctGrp.size(); fi++)
		{

			// special case for vertex
			if (dim == 0)
			{
				value = StdFuncIntegralOnVertex(solution, fctGrp[fi], ssGrp[si]);
			}
			else
			{
				SmartPtr<IIntegrand<number, TGridFunction::dim> > spIntegrand
					= make_sp(new StdFuncIntegrand<TGridFunction>(solution, fctGrp[fi]));

				if (dim == worldDim)
					value = IntegrateSubset<TGridFunction, worldDim>(spIntegrand, solution, ssGrp[si], 1, "best");
				else if (dim == worldDim-1 && worldDim > 1)
					value = IntegrateSubset<TGridFunction, (worldDim>1 ? worldDim-1 : 1)>(spIntegrand, solution, ssGrp[si], 1, "best");
				else if (dim == worldDim-2 && worldDim > 2)
					value = IntegrateSubset<TGridFunction, (worldDim>2 ? worldDim-2 : 1)>(spIntegrand, solution, ssGrp[si], 1, "best");
				else {UG_THROW("Unknown dim (" << dim << ") or worldDim (" << worldDim << ").");}
			}
		#ifdef UG_PARALLEL
			// sum over processes
			if (pcl::NumProcs() > 1)
			{
				pcl::ProcessCommunicator com;
				number local = value;
				com.allreduce(&local, &value, 1, PCL_DT_DOUBLE, PCL_RO_SUM);
			}

			// check if this proc is output proc
			if (GetLogAssistant().is_output_process())
			{
		#endif

			// construct outFile name
			std::ostringstream ofnss(outFileName, std::ios_base::app);
			ofnss << "_" << ssGrp.name(si) << "_" << fctGrp.name(fi) << outFileExt;

			// create if first time step, append otherwise
			std::ofstream outFile;
			if (time == 0.0) outFile.open(ofnss.str().c_str(), std::ios_base::out);
			else outFile.open(ofnss.str().c_str(), std::ios_base::app);

			// write measurement
			try {outFile << time << "\t" << value/vol << "\n";}
			UG_CATCH_THROW("Output file " << ofnss.str() << " could not be written to.");
			outFile.close();

		#ifdef UG_PARALLEL
			}
		#endif
		}
	}

	return value/vol;
}



} // namespace ug
} // namespace neuro_collection
