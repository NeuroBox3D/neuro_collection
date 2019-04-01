/*
 * ca_wave_util_impl.h
 *
 *  Created on: 2017-09-18
 *      Author: mbreit
 */

#include "ca_wave_util.h"

#include "common/error.h"                       // UG_THROW etc.
#include "lib_disc/common/multi_index.h"        // DoFIndex
#include "lib_grid/common_attachments.h"        // for aPosition2
#include "lib_grid/algorithms/element_side_util.h"  // for GetOpposingSide

#include <sstream>
#include <vector>
#include <algorithm>  // for std::sort


namespace ug {
namespace neuro_collection {

template <typename TGridFunction, typename TRyRImpl>
number maxRyRFluxDensity
(
	ConstSmartPtr<TGridFunction> u,
	const char* fctNames,
	const char* subsetNames,
	ConstSmartPtr<TRyRImpl> ryr
)
{
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
			values[TRyRImpl::_CCYT_] = ryr->scale_input(TRyRImpl::_CCYT_) * DoFRef(*u, ind[0]);

			UG_COND_THROW(!dd->is_def_in_subset(fi_ce, si),
				"Function " << fi_ce << " is not defined on subset " << si << ".");
			numInd = dd->dof_indices(vrt, fi_ce, ind, true, true);
			UG_ASSERT(numInd == 1, "More (or less) than one function index found on a vertex!");
			values[TRyRImpl::_CER_] = ryr->scale_input(TRyRImpl::_CER_) * DoFRef(*u, ind[0]);

			UG_COND_THROW(!dd->is_def_in_subset(fi_c1, si),
				"Function " << fi_c1 << " is not defined on subset " << si << ".");
			numInd = dd->dof_indices(vrt, fi_c1, ind, true, true);
			UG_ASSERT(numInd == 1, "More (or less) than one function index found on a vertex!");
			values[TRyRImpl::_C1_] = ryr->scale_input(TRyRImpl::_C1_) * DoFRef(*u, ind[0]);

			UG_COND_THROW(!dd->is_def_in_subset(fi_c2, si),
				"Function " << fi_c2 << " is not defined on subset " << si << ".");
			numInd = dd->dof_indices(vrt, fi_c2, ind, true, true);
			UG_ASSERT(numInd == 1, "More (or less) than one function index found on a vertex!");
			values[TRyRImpl::_C2_] = ryr->scale_input(TRyRImpl::_C2_) * DoFRef(*u, ind[0]);

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
	SmartPtr<TGridFunction> u,
	const char* fctNames,
	const char* subsetName,
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

	UG_COND_THROW(fctGrp.size() < 1 || fctGrp.size() > 2,
		"Function group must have exactly 1 or exactly 2 entries.\n"
		"Make sure you provide either exactly 1 function name (ca_cyt)\n"
		"or exactly 2 functions in the following order: c1, c2.");

	// get subset indices via SubsetGroup
	SubsetGroup ssGrp;
	try {ssGrp = u->subset_grp_by_name(subsetName);}
	UG_CATCH_THROW("Something went wrong during creation of subset group for \n"
		"membrane subset names '" << subsetName << "' and provided solution.");
	UG_COND_THROW(ssGrp.size() < 1, "Subset group must have at least 1 entry.\n"
			"Make sure you provide at least 1 subset name for ER membrane subset.");
	const size_t nSs = ssGrp.size();

	double xMax = -std::numeric_limits<double>::max();

	// RyR channel state mode:
	// wave front position is calculated depending on RyR channel open probability
	if (fctGrp.size() == 2)
	{
		size_t fi_c1 = fctGrp[0];
		size_t fi_c2 = fctGrp[1];

		Vertex* maxVrt = NULL;
		number pOpenMax = 0.0;

		for (size_t s = 0; s < nSs; ++s)
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
				{
					xMax = xCoord;
					pOpenMax = pOpen;
					maxVrt = vrt;
				}
			}
		}

		// interpolate exact threshold position
		if (maxVrt)
		{
			MultiGrid& mg = *u->domain()->grid();
			MGSubsetHandler& sh = *u->domain()->subset_handler();

			typename Grid::traits<Edge>::secure_container edges;
			mg.associated_elements(edges, maxVrt);
			size_t esz = edges.size();
			for (size_t e = 0; e < esz; ++e)
			{
				Edge* ed = edges[e];
				Vertex* other = GetOpposingSide(mg, ed, maxVrt);
				if (!ssGrp.contains(sh.get_subset_index(other)))
					continue;

				number xCoord = aaPos[other][0];
				if ((double) xCoord <= xMax)
					continue;

				// now we have the right neighbor vertex (if it exists)
				// per def of the max, it has a value below thresh
				std::vector<DoFIndex> ind;
				dd->dof_indices(other, fi_c1, ind, true, true);
				number c1 = DoFRef(*u, ind[0]);
				dd->dof_indices(other, fi_c2, ind, true, true);
				number c2 = DoFRef(*u, ind[0]);
				number pOpen = 1.0 - (c1 + c2);

				xMax = (xMax*(pOpen-thresh) - xCoord*(pOpenMax - thresh)) / (pOpen-pOpenMax);
				break;
			}
		}
	}

	// cytosolic calcium mode:
	// wave front position is calculated depending on cytosolic calcium concentration
	else if (fctGrp.size() == 1)
	{
		size_t fi_cc = fctGrp[0];

		Vertex* maxVrt = NULL;
		number ccMax = 0.0;

		for (size_t s = 0; s < nSs; ++s)
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

				UG_COND_THROW(!dd->is_def_in_subset(fi_cc, si),
					"Function " << fi_cc << " is not defined on subset " << si << ".");
				size_t numInd = dd->dof_indices(vrt, fi_cc, ind, true, true);
				UG_COND_THROW(numInd != 1, "More (or less) than one function index found on a vertex!");
				const number cc = DoFRef(*u, ind[0]);

				// update xMax if necessary
				if (cc > thresh)
				{
					xMax = xCoord;
					ccMax = cc;
					maxVrt = vrt;
				}
			}
		}

		// interpolate exact threshold position
		if (maxVrt)
		{
			MultiGrid& mg = *u->domain()->grid();
			MGSubsetHandler& sh = *u->domain()->subset_handler();

			typename Grid::traits<Edge>::secure_container edges;
			mg.associated_elements(edges, maxVrt);
			size_t esz = edges.size();
			for (size_t e = 0; e < esz; ++e)
			{
				Edge* ed = edges[e];
				Vertex* other = GetOpposingSide(mg, ed, maxVrt);
				if (!ssGrp.contains(sh.get_subset_index(other)))
					continue;

				number xCoord = aaPos[other][0];
				if ((double) xCoord <= xMax)
					continue;

				// now we have the right neighbor vertex (if it exists)
				// per def of the max, it has a value below thresh
				std::vector<DoFIndex> ind;
				dd->dof_indices(other, fi_cc, ind, true, true);
				const number cc = DoFRef(*u, ind[0]);

				xMax = (xMax*(cc-thresh) - xCoord*(ccMax - thresh)) / (cc-ccMax);
				break;
			}
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





template <typename TDomain, typename TAlgebra>
WaveProfileExporter<TDomain, TAlgebra>::WaveProfileExporter
(
	SmartPtr<ApproximationSpace<TDomain> > approxSpace,
	const char* fctNames,
	const char* subsetNames,
	const std::string& fileBaseName
)
: m_spApprox(approxSpace),
  m_vFct(TokenizeString(fctNames)),
  m_vSs(TokenizeString(subsetNames)),
  m_fileName(fileBaseName),
  m_vvvDoFSeries(m_vSs.size()),
  m_vvXPos(m_vSs.size())
{
	typedef typename TDomain::position_attachment_type pos_attach_type;
	typedef typename DoFDistribution::traits<Vertex>::const_iterator vrt_it;

	ConstSmartPtr<MGSubsetHandler> sh = m_spApprox->domain()->subset_handler();
	ConstSmartPtr<DoFDistribution> dd = m_spApprox->dof_distribution(GridLevel(), false);

	SubsetGroup ssg(sh);
	try {ssg.add(m_vSs);}
	UG_CATCH_THROW("Could not add all subsets to WaveProfileExporter.");

	FunctionGroup fg(m_spApprox->function_pattern());
	try {fg.add(m_vFct);}
	UG_CATCH_THROW("Could not add all functions to WaveProfileExporter.");

	pos_attach_type& aPos = GetDefaultPositionAttachment<pos_attach_type>();
	Grid::VertexAttachmentAccessor<pos_attach_type> aaPos(*m_spApprox->domain()->grid(), aPos);

	// iterate all subsets
	size_t nfct = fg.size();
	size_t nss = ssg.size();
	for (size_t s = 0; s < nss; ++s)
	{
		const int si = ssg[s];
		m_vvvDoFSeries[s].resize(nfct);

		// iterate all surface vertices in subset
		vrt_it it = dd->begin<Vertex>(si);
		vrt_it itEnd = dd->end<Vertex>(si);
		for (; it != itEnd; ++it)
		{
			Vertex* vrt = *it;

			// save vertex
			m_vvXPos[s].push_back(aaPos[vrt][0]);

			// iterate functions
			for (size_t f = 0; f < nfct; ++f)
			{
				const size_t fct = fg[f];

				// save DoFIndices
				std::vector<DoFIndex>& vdi = m_vvvDoFSeries[s][f];
				UG_COND_THROW(!dd->is_def_in_subset(fct, si), "Function '" << m_vFct[f]
					<< "' is not defined on subset '" << m_vSs[s] << "'.");
				dd->inner_dof_indices(vrt, fct, vdi, false);
			}
		}

		// sort DoFs from left to right
		CmpVrtPos cmp(m_vvXPos[s]);
		size_t nVrt = m_vvXPos[s].size();
		std::vector<size_t> perm(nVrt);
		for (size_t i = 0; i < nVrt; ++i)
			perm[i] = i;
		std::sort(perm.begin(), perm.end(), cmp);

		{
			std::vector<number> sortedPos(nVrt);
			for (size_t i = 0; i < nVrt; ++i)
				sortedPos[i] = m_vvXPos[s][perm[i]];
			sortedPos.swap(m_vvXPos[s]);
		}

		for (size_t f = 0; f < nfct; ++f)
		{
			std::vector<DoFIndex>& vdi = m_vvvDoFSeries[s][f];
			UG_COND_THROW(vdi.size() != nVrt, "nVrt and nDoF size mismatch");
			std::vector<DoFIndex> vdiNew(nVrt);
			for (size_t i = 0; i < nVrt; ++i)
				vdiNew[i] = vdi[perm[i]];
			vdiNew.swap(vdi);
		}
	}
}


#ifdef UG_PARALLEL
struct MyCompare
{
	MyCompare(const std::vector<number>& _v) : v(_v) {};

	bool operator()(const size_t& a, const size_t& b)
	{
		return v[a] < v[b];
	}

	private:
		const std::vector<number>& v;
};

static void writeParallelFile(char* data, size_t dataSize, char* fileName, number minX)
{
	pcl::ProcessCommunicator pc;

	MPI_Status status;
	MPI_Comm m_mpiComm = pc.get_mpi_communicator();
	MPI_File fh;

	// open file
	if (MPI_File_open(m_mpiComm, fileName, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh))
		UG_THROW("Unable to open "<< fileName << ".");


	// find out correct order of procs
	size_t np = pcl::NumProcs();
	std::vector<number> allMinX(np);
	pc.allgather(&minX, 1, PCL_DT_DOUBLE, &allMinX[0], 1, PCL_DT_DOUBLE);
	std::vector<size_t> rankOrder(np);
	for (size_t i = 0; i < np; ++i)
		rankOrder[i] = i;

	MyCompare cmp(allMinX);
	std::sort(rankOrder.begin(), rankOrder.end(), cmp);


	// calculate offsets for each proc
	unsigned long mySize = dataSize;
	std::vector<unsigned long> allSizes(np);
	pc.allgather(&mySize, 1, PCL_DT_UNSIGNED_LONG, &allSizes[0], 1, PCL_DT_UNSIGNED_LONG);

	size_t myRank = pcl::ProcRank();
	size_t offset = 0;
	for (size_t i = 0; i < np; ++i)
	{
		if (rankOrder[i] == myRank)
			break;

		offset += allSizes[rankOrder[i]];
	}


	// write data at correct offset
	MPI_File_seek(fh, offset, MPI_SEEK_SET);
	MPI_File_write(fh, data, mySize, MPI_BYTE, &status);


	// close file
	MPI_File_close(&fh);
}
#endif


static void writeSerialFile(const std::string& data, const std::string& fileName)
{
	std::ofstream ofs(fileName.c_str());
	try {ofs << data;}
	UG_CATCH_THROW("Output file '" << fileName << "' could not be written to.");
	ofs.close();
}



template <typename TDomain, typename TAlgebra>
void WaveProfileExporter<TDomain, TAlgebra>::
exportWaveProfileX(ConstSmartPtr<gf_type> u, number time)
{
	const size_t nsi = m_vvvDoFSeries.size();
	for (size_t s = 0; s < nsi; ++s)
	{
		const size_t nfct = m_vvvDoFSeries[s].size();
		for (size_t f = 0; f < nfct; ++f)
		{
			const std::vector<DoFIndex>& vdi = m_vvvDoFSeries[s][f];

			// construct file name
			std::ostringstream ossFn;
			ossFn << m_fileName << "_" << m_vSs[s] << "_" << m_vFct[f] << "_" << time << ".dat";

			// write this proc's values to buffer
			std::ostringstream ossVal;
			size_t nVrt = m_vvXPos[s].size();
			for (size_t i = 0; i < nVrt; ++i)
				ossVal << m_vvXPos[s][i] << "\t" << DoFRef(*u, vdi[i]) << "\n";

#ifdef UG_PARALLEL
			if (pcl::NumProcs() > 1)
				// some old MPI implementations need non-const char* for MPI_File_open
				// and MPI_File_write, so we do a dirty const_cast here
				// (which is of no consequence if the MPI implementation uses const char*)
				writeParallelFile(const_cast<char*>(ossVal.str().c_str()), ossVal.str().size(),
					const_cast<char*>(ossFn.str().c_str()), nVrt ? m_vvXPos[s][0] : 0.0);
			else
#endif
			writeSerialFile(ossVal.str(), ossFn.str());
		}
	}
}



} // namespace ug
} // namespace neuro_collection
