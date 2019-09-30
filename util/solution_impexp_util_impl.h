/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2012-07-31
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

#include <cstddef>                                                  // for size_t, NULL
#include <fstream>                                                  // for operator<<, basic_ostream, char...
#include <limits>                                                   // for numeric_limits, numeric_limits<...
#include <vector>
#include <sstream>
#include <fstream>

#include "common/assert.h"                                          // for UG_ASSERT
#include "lib_disc/common/function_group.h"                         // for FunctionGroup
#include "lib_disc/common/multi_index.h"                            // for DoFIndex, DoFRef
#include "lib_disc/dof_manager/dof_distribution.h"                  // for DoFDistribution, DoFDistributio...
#include "lib_disc/function_spaces/approximation_space.h"           // for ApproximationSpace
#include "lib_disc/function_spaces/dof_position_util.h"             // for InnerDoFPosition
#include "lib_disc/function_spaces/grid_function.h"                 // for GridFunction
#include "lib_disc/local_finite_element/local_finite_element_id.h"  // for LFEID
#include "lib_grid/grid/grid_base_objects.h"                        // for Vertex (ptr only), GridBaseObje...
#include "lib_grid/tools/subset_group.h"                            // for SubsetGroup

// configuration file for compile options
#include "nc_config.h"

#ifdef NC_WITH_MPM
	#include "../../plugins/MembranePotentialMapping/vm2ug.h"
#endif


namespace ug {
namespace neuro_collection {


// forward declaration
template <typename TDomain>
number computeVolume(ConstSmartPtr<ApproximationSpace<TDomain> > approx, const int subset);


// /////////////////////// //
// export solution command //
// /////////////////////// //

// helper function
template <typename TGridFunction, typename TBaseElem>
static void exportSolution
(
	SmartPtr<TGridFunction> solution,
	size_t si,
	size_t fi,
	std::ofstream& ofs
)
{
	typedef typename TGridFunction::domain_type domain_type;

	// retrieve domain and dofDistr from approxSpace
	ConstSmartPtr<domain_type> domain = solution->approx_space()->domain();
	ConstSmartPtr<DoFDistribution> dofDistr = solution->dof_distribution();
	const LFEID lfeid = dofDistr->lfeid(fi);

	//	get elem iterator for current subset and elem type
	typedef typename DoFDistribution::traits<TBaseElem>::const_iterator itType;
	itType iter = dofDistr->template begin<TBaseElem>(si);
	itType iterEnd = dofDistr->template end<TBaseElem>(si);

	// loop over all elems
	for (; iter != iterEnd; ++iter)
	{
		// get current vertex
		TBaseElem* elem = *iter;

		// get coords
		std::vector<typename domain_type::position_type> coords;
		InnerDoFPosition<domain_type>(coords, elem, *domain, lfeid);

		// get multi-indices
		std::vector<DoFIndex> multInd;
		dofDistr->inner_dof_indices(elem, fi, multInd);

		UG_ASSERT(coords.size() == multInd.size(), "#DoF mismatch");

		// get values of DoFs
		number val;
		size_t nDof = multInd.size();
		for (size_t dof = 0; dof < nDof; ++dof)
		{
			val = DoFRef(*solution, multInd[dof]);

			// write solution to file
			for (size_t i = 0; i < coords[dof].size(); ++i)
				ofs << coords[dof][i] << " ";
			ofs << val << std::endl;
		}
	}
}


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

	// retrieve dofDistr from solution
	ConstSmartPtr<DoFDistribution> dofDistr = solution->dof_distribution();

	// get subset group to be measured on (if none provided: take all)
	SubsetGroup ssGrp;
	if (!*subsetNames)
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
	if (!*functionNames)
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
		ofnss << "_" << time << "_" << fctGrp.name(fi);

		// create if first time step, append otherwise
		std::ofstream outFile;
		outFile.precision(std::numeric_limits<number>::digits10);
		outFile.open(ofnss.str().c_str(), std::ios_base::out);
		if (!outFile.is_open())
			UG_THROW("Output file '" << ofnss.str() << "' could not be opened.")

		try
		{
			// loop subsets
			for (size_t si = 0; si < ssGrp.size(); si++)
			{
				if (domain_type::dim-1 >= VERTEX && dofDistr->max_fct_dofs(fctGrp[fi], VERTEX, ssGrp[si]) > 0)
					exportSolution<TGridFunction, Vertex>(solution, ssGrp[si], fctGrp[fi], outFile);
				if (domain_type::dim-1 >= EDGE && dofDistr->max_fct_dofs(fctGrp[fi], EDGE, ssGrp[si]) > 0)
					exportSolution<TGridFunction, Edge>(solution, ssGrp[si], fctGrp[fi], outFile);
				if (domain_type::dim-1 >= FACE && dofDistr->max_fct_dofs(fctGrp[fi], FACE, ssGrp[si]) > 0)
					exportSolution<TGridFunction, Face>(solution, ssGrp[si], fctGrp[fi], outFile);
				if (domain_type::dim-1 >= VOLUME && dofDistr->max_fct_dofs(fctGrp[fi], VOLUME, ssGrp[si]) > 0)
					exportSolution<TGridFunction, Volume>(solution, ssGrp[si], fctGrp[fi], outFile);
			}
		}
		UG_CATCH_THROW("Output file '" << ofnss.str() << "' could not be written to.");

		outFile.close();
	}

	return;
}



// /////////////////////// //
// import solution command //
// /////////////////////// //

#ifdef NC_WITH_MPM
// helper function
template <typename TGridFunction, typename TBaseElem>
static void importSolution
(
	SmartPtr<TGridFunction> solution,
	size_t si,
	size_t fi,
	const membrane_potential_mapping::Mapper<TGridFunction::domain_type::dim, number>& mapper
)
{
	typedef typename TGridFunction::domain_type domain_type;

	// retrieve domain and dofDistr from approxSpace
	ConstSmartPtr<domain_type> domain = solution->approx_space()->domain();
	ConstSmartPtr<DoFDistribution> dofDistr = solution->dof_distribution();
	const LFEID lfeid = dofDistr->lfeid(fi);

	//	get elem iterator for current subset and elem type
	typedef typename DoFDistribution::traits<TBaseElem>::const_iterator itType;
	itType iter = dofDistr->template begin<TBaseElem>(si);
	itType iterEnd = dofDistr->template end<TBaseElem>(si);

	// loop over all elems
	for (; iter != iterEnd; ++iter)
	{
		// get current vertex
		TBaseElem* elem = *iter;

		// get coords
		std::vector<typename domain_type::position_type> coords;
		InnerDoFPosition<domain_type>(coords, elem, *domain, lfeid);

		// get multi-indices
		std::vector<DoFIndex> multInd;
		dofDistr->inner_dof_indices(elem, fi, multInd);

		UG_ASSERT(coords.size() == multInd.size(), "#DoF mismatch");

		// get values of DoFs
		number val;
		size_t nDof = multInd.size();
		for (size_t dof = 0; dof < nDof; ++dof)
		{
			// get value from provider
			try {val = mapper.get_data_from_nearest_neighbor(coords[dof]);}
			UG_CATCH_THROW("No value could be retrieved for DoF at " << coords[dof]);

			DoFRef(*solution, multInd[dof]) = val;
		}
	}
}
#endif

template <typename TGridFunction>
void importSolution
(
	SmartPtr<TGridFunction> solution,
	const char* subsetNames,
	const char* functionNames,
	const char* inFileBaseName
)
{
#ifndef NC_WITH_MPM
	UG_THROW("importSolution uses functionality from the MembranePotentialMapping plugin, "
		"but was not compiled with it.");
#else
	typedef typename TGridFunction::domain_type domain_type;
	typedef typename domain_type::position_type pos_type;

	ConstSmartPtr<domain_type> domain = solution->approx_space()->domain();
	ConstSmartPtr<DoFDistribution> dofDistr = solution->dof_distribution();

	// get subset group to be measured on
	SubsetGroup ssGrp;
	if (!*subsetNames)
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
	if (!*functionNames)
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
		// construct inFile name
		std::ostringstream ofnss(inFileBaseName, std::ios_base::app);
		ofnss << "_" << fctGrp.name(fi);

		// read values from file and fill mapper structure with it
		membrane_potential_mapping::Mapper<domain_type::dim, number> valueProvider;
		try {valueProvider.build_tree(ofnss.str());}
		UG_CATCH_THROW("Underlying mapper object could not build its tree "
					   "on given file (" << ofnss.str() << ").");

		// loop subsets
		for (size_t si = 0; si < ssGrp.size(); si++)
		{
			if (domain_type::dim-1 >= VERTEX && dofDistr->max_fct_dofs(fctGrp[fi], VERTEX, ssGrp[si]) > 0)
				importSolution<TGridFunction, Vertex>(solution, ssGrp[si], fctGrp[fi], valueProvider);
			if (domain_type::dim-1 >= EDGE && dofDistr->max_fct_dofs(fctGrp[fi], EDGE, ssGrp[si]) > 0)
				importSolution<TGridFunction, Edge>(solution, ssGrp[si], fctGrp[fi], valueProvider);
			if (domain_type::dim-1 >= FACE && dofDistr->max_fct_dofs(fctGrp[fi], FACE, ssGrp[si]) > 0)
				importSolution<TGridFunction, Face>(solution, ssGrp[si], fctGrp[fi], valueProvider);
			if (domain_type::dim-1 >= VOLUME && dofDistr->max_fct_dofs(fctGrp[fi], VOLUME, ssGrp[si]) > 0)
				importSolution<TGridFunction, Volume>(solution, ssGrp[si], fctGrp[fi], valueProvider);
		}
	}
#endif
}

} // namespace neuro_collection
} // namespace ug
