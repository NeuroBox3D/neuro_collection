/*
 * solution_impexp_util_impl.h
 *
 *  Created on: 31.07.2012
 *      Author: mbreit
 */

#include <cstdlib>
#include <limits>
#include <vector>
#include <sstream>
#include <fstream>


#ifdef NC_WITH_MPM
	#include "../plugins/MembranePotentialMapping/vm2ug_rework.h"
#endif


namespace ug {
namespace neuro_collection {


// forward declaration
template <typename TDomain>
number computeVolume(ConstSmartPtr<ApproximationSpace<TDomain> > approx, const int subset);


// /////////////////////// //
// export solution command //
// /////////////////////// //

template <typename TGridFunction>
void exportSolution
(
	SmartPtr<TGridFunction> solution,
	const number time,
	const char* subsetNames,
	const char* functionNames,
	const char* outFileName
)
{
	typedef typename TGridFunction::domain_type domain_type;

	// retrieve domain, dofDistr and dim from approxSpace
	ConstSmartPtr<domain_type> domain = solution->approx_space()->domain();
	ConstSmartPtr<DoFDistribution> dofDistr = solution->dof_distribution();
	typedef typename DoFDistribution::traits<Vertex>::const_iterator itType;

	// get subset group to be measured on (if none provided: take all)
	SubsetGroup ssGrp;
	if (subsetNames == NULL)
	{
		ssGrp.set_subset_handler(dofDistr->subset_handler());
		ssGrp.add_all();
	}
	else
	{
		try {ssGrp = dofDistr->subset_grp_by_name(subsetNames);}
		UG_CATCH_THROW("At least one of the subsets in '" << subsetNames
				<< "' is not contained in the approximation space (or something else was wrong).");
	}

	// get function group to be measured (if none provided: take all)
	FunctionGroup fctGrp;
	if (functionNames == NULL)
	{
		fctGrp.set_function_pattern(dofDistr->function_pattern());
		fctGrp.add_all();
	}
	else
	{
		try {fctGrp = dofDistr->fct_grp_by_name(functionNames);}
		UG_CATCH_THROW("At least one of the functions in '" << functionNames
						<< "' is not contained in the approximation space (or something else was wrong).");
	}

	// loop functions
	for (size_t fi = 0; fi < fctGrp.size(); fi++)
	{
		// construct outFile name
		std::ostringstream ofnss(outFileName, std::ios_base::app);
		ofnss << "_" << fctGrp.name(fi) << "_" << time;

		// create if first time step, append otherwise
		std::ofstream outFile;
		outFile.precision(std::numeric_limits<number>::digits10);
		outFile.open(ofnss.str().c_str(), std::ios_base::out);

		UG_LOGN(ofnss.str());
		// loop subsets
		for (size_t si = 0; si < ssGrp.size(); si++)
		{
			//	get vertex iterator for current subset
			itType iterBegin = dofDistr->template begin<Vertex>(ssGrp[si]);
			itType iterEnd = dofDistr->template end<Vertex>(ssGrp[si]);

			// loop over all vertices
			for (itType iter = iterBegin; iter != iterEnd; ++iter)
			{
				// get current vertex
				Vertex* vert = *iter;

				// get coords
				typename domain_type::position_accessor_type posAcc = domain->position_accessor();
				typename domain_type::position_type coords = posAcc[vert];

				// get multi-indices
				std::vector<DoFIndex> multInd;
				dofDistr->dof_indices(vert, fctGrp[fi], multInd);
				if (multInd.size() != 1)
					UG_THROW("Number of multi-indices returned for vertex "
							 "of function '" << fctGrp.name(fi) << "' on subset '"
							 << ssGrp.name(si) << "' not exactly one ("
							 << multInd.size() << " instead).");

				// get value
				const number val = DoFRef(*solution, multInd[0]);

				// write solution to file
				try
				{
					for (size_t i = 0; i < coords.size(); ++i)
						outFile << coords[i] << " ";
					outFile << val << std::endl;
				}
				UG_CATCH_THROW("output file" << ofnss.str() << "could not be written to.");
			}
		}

		outFile.close();
	}

	return;
}



#ifdef NC_WITH_MPM
// /////////////////////// //
// import solution command //
// /////////////////////// //

template <typename TGridFunction>
void importSolution
(
	SmartPtr<TGridFunction> solution,
	const char* subsetNames,
	const char* functionName,
	const char* inFileName
)
{
	typedef typename TGridFunction::domain_type domain_type;
	typedef typename domain_type::position_type pos_type;
	typedef typename domain_type::position_accessor_type pos_acc_type;

	ConstSmartPtr<domain_type> domain = solution->approx_space()->domain();
	ConstSmartPtr<DoFDistribution> dofDistr = solution->dof_distribution();
	const pos_acc_type& aaPos = domain->position_accessor();

	// read values from file and fill Vm2Ug structure with it
	// Vm2uG<std::string> valueProvider(inFileName, "", true);
	static const size_t dim = TGridFunction::dim;
	membrane_potential_mapping::Mapper<dim, number> valueProvider;

	//try {valueProvider.buildTree("");}
	try {valueProvider.build_tree(inFileName, " ");}
	UG_CATCH_THROW("Underlying Vm2uG object could not build its tree on given file.");

	// get subset group to be measured on
	SubsetGroup ssGrp;
	try {ssGrp = solution->subset_grp_by_name(subsetNames);}
	UG_CATCH_THROW("At least one of the subsets in '" << subsetNames
					<< "' is not contained in the approximation space (or something else was wrong).");

	// get function group to be measured
	FunctionGroup fctGrp;
	try {fctGrp = solution->fct_grp_by_name(functionName);}
	UG_CATCH_THROW("The function '" << functionName
					<< "' is not contained in the approximation space (or something else was wrong).");
	if (fctGrp.size() != 1)
	{
		UG_THROW("The functions argument must contain exactly one valid function name, "
				 "but it contains " << fctGrp.size() << ".");
	}
	size_t fct = fctGrp[0];

	// loop subsets
	for (size_t si = 0; si < ssGrp.size(); si++)
	{
		int ssi = ssGrp[si];

		// loop vertices of subset
		DoFDistribution::traits<Vertex>::const_iterator iter, iterEnd;
		iter = dofDistr->begin<Vertex>(ssi);
		iterEnd = dofDistr->end<Vertex>(ssi);

		for (; iter != iterEnd; ++iter)
		{
			Vertex* vrt = *iter;

			// get coordinates for vertex
			const pos_type& coords = aaPos[vrt];

			// get value from provider
			double val;
			try
			{
				val = valueProvider.get_data_from_nearest_neighbor(coords);
			}
			UG_CATCH_THROW("No value for the vertex could be retrieved from file.");

			// get index for vertex
			std::vector<size_t> ind;
			if (!dofDistr->is_def_in_subset(fct, ssi))
				{UG_THROW("Function " << fct << " is not defined on subset " << ssi << ".");}
			dofDistr->inner_algebra_indices_for_fct(vrt, ind, false, fct);
			UG_ASSERT(ind.size() == 1, "More (or less) than one function index found on a vertex!");

			// assign solution at vertex to corresponding index in vector
			(*solution)[ind[0]] = val;
		}
	}
}
#endif

} // namespace neuro_collection
} // namespace ug
