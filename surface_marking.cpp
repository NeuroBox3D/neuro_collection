/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2016-03-30
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

#include "surface_marking.h"

#include "lib_disc/dof_manager/dof_distribution.h"          // for DoFDistribution, DoFDistribut...
#include "lib_disc/domain.h"                                // for Domain1d, Domain2d, Domain3d
#include "lib_disc/domain_traits.h"                         // for domain_traits
#include "lib_disc/function_spaces/approximation_space.h"   // for ApproximationSpace
#include "lib_disc/function_spaces/error_indicator_util.h"  // for ComputeMinMaxTotal
#include "lib_grid/grid/grid.h"                             // for Grid::traits, Grid::associated_elements
#include "lib_grid/grid/grid_base_objects.h"                // for Vertex
#include "lib_grid/multi_grid.h"                            // for MultiGrid
#include "lib_grid/refinement/refiner_interface.h"          // for IRefiner, RefinementMark::RM_...
#include "lib_grid/tools/grid_level.h"                      // for GridLevel
#include "lib_grid/tools/subset_handler_interface.h"        // for ISubsetHandler


namespace ug {
namespace neuro_collection {


template <typename TDomain>
void SurfaceMarking<TDomain>::add_surface(int surf_si, int adj_vol_si)
{
	m_vpSurfaces.push_back(std::make_pair(surf_si, adj_vol_si));
}

template <typename TDomain>
void SurfaceMarking<TDomain>::remove_surface(int surf_si, int adj_vol_si)
{
	std::pair<int,int> pair = std::make_pair(surf_si, adj_vol_si);
	std::vector<std::pair<int,int> >::iterator match =
		std::find(m_vpSurfaces.begin(), m_vpSurfaces.end(), pair);

	if (match != m_vpSurfaces.end())
		m_vpSurfaces.erase(match);
}


template <typename TDomain>
void SurfaceMarking<TDomain>::add_surface(const std::string& surf_si, const std::string& adj_vol_si)
{
	UG_COND_THROW(!m_spDom.valid(), "Adding surface marking pairs via subset names only works\n"
		"when the SurfaceMarking object is given a domain at construction time.")

	ConstSmartPtr<ISubsetHandler> spSH = m_spDom->subset_handler();
	try
	{
		const int srfInd = spSH->get_subset_index(surf_si.c_str());
		const int volInd = spSH->get_subset_index(adj_vol_si.c_str());
		m_vpSurfaces.push_back(std::make_pair(srfInd, volInd));
	}
	UG_CATCH_THROW("Cannot convert both subset names to indices.");
}

template <typename TDomain>
void SurfaceMarking<TDomain>::remove_surface(const std::string& surf_si, const std::string& adj_vol_si)
{
	UG_COND_THROW(!m_spDom.valid(), "Adding surface marking pairs via subset names only works \n"
		"when the SurfaceMarking object is given a domain at construction time.")

	int srfInd = 0;
	int volInd = 0;
	ConstSmartPtr<ISubsetHandler> spSH = m_spDom->subset_handler();
	try
	{
		srfInd = spSH->get_subset_index(surf_si.c_str());
		volInd = spSH->get_subset_index(adj_vol_si.c_str());
	}
	UG_CATCH_THROW("Cannot convert both subset names to indices.");


	std::pair<int,int> pair = std::make_pair(srfInd, volInd);
	std::vector<std::pair<int,int> >::iterator match =
		std::find(m_vpSurfaces.begin(), m_vpSurfaces.end(), pair);

	if (match != m_vpSurfaces.end())
		m_vpSurfaces.erase(match);
}



template <typename TDomain>
void SurfaceMarking<TDomain>::mark
(
	typename base_type::elem_accessor_type& aaError,
	IRefiner& refiner,
	ConstSmartPtr<DoFDistribution> dd
)
{
	UG_COND_THROW(m_tol == -1.0 && m_max_level == (size_t) -1,
		"SurfaceMarking needs to be passed tolerance and maximal refinement level in constructor \n"
		"when used as an error-based refinement strategy.");

	number minElemErr;
	number maxElemErr;
	number errTotal;
	size_t numElem;

	ComputeMinMaxTotal(aaError, dd, minElemErr, maxElemErr, errTotal, numElem);
	if (errTotal <= m_tol || dd->multi_grid()->num_levels() > m_max_level)
		return;

	mark(refiner, dd);
}


template <typename TDomain>
void SurfaceMarking<TDomain>::mark_without_error
(
	SmartPtr<IRefiner> refiner,
	SmartPtr<ApproximationSpace<TDomain> > approx
)
{
	// convert subsets given by name to indices using the given approxSpace
	ConstSmartPtr<ISubsetHandler> sh = approx->domain()->subset_handler();

	// get surface dof distribution
	ConstSmartPtr<DoFDistribution> dd = approx->dof_distribution(GridLevel(), false);

	mark(*refiner, dd);
}



template <typename TDomain>
void SurfaceMarking<TDomain>::mark(IRefiner& refiner, ConstSmartPtr<DoFDistribution> dd)
{
	typedef typename MultiGrid::traits<typename domain_traits<TDomain::dim>::element_type>::secure_container elem_list_type;
	typedef typename domain_traits<TDomain::dim>::side_type side_type;
	typedef typename MultiGrid::traits<side_type>::secure_container side_list_type;
	typedef typename DoFDistribution::traits<Vertex>::const_iterator iter_type;

	Grid* grid = refiner.grid();
	ConstSmartPtr<ISubsetHandler> sh = dd->subset_handler();

	size_t sz = m_vpSurfaces.size();
	for (size_t i = 0; i < sz; ++i)
	{
		int surf_si = m_vpSurfaces[i].first;
		int vol_si = m_vpSurfaces[i].second;
		iter_type iter = dd->template begin<Vertex>(surf_si);
		iter_type iterEnd = dd->template end<Vertex>(surf_si);

		//	loop elements for marking
		for (; iter != iterEnd; ++iter)
		{
			Vertex* v = *iter;
			elem_list_type el;
			grid->associated_elements(el, v);
			size_t el_sz = el.size();
			for (size_t j = 0; j < el_sz; ++j)
			{
				int elem_si = sh->get_subset_index(el[j]);
				if (elem_si == vol_si)
				{
					// mark the element for anisotropic refinement
					refiner.mark(el[j], RM_REFINE);
				}
			}
		}
	}
/*
	if (m_vIntf.size())
	{
		HangingNodeRefinerBase<MGSelector>* hnr = dynamic_cast<HangingNodeRefinerBase<MGSelector>*>(&refiner);
		if (hnr)
		{
			SmartPtr<InterfaceRefMarkAdjuster> intfRMA = make_sp(new InterfaceRefMarkAdjuster());
			intfRMA->set_subset_handler(sh);
			intfRMA->add_interfaces(m_vIntf);
			hnr->add_ref_mark_adjuster(intfRMA);
		}
		else
		{
			UG_LOGN("Warning: SurfaceMarking is supposed to work with an instance of HangingNodeRefiner\n"
					"         and might not work properly with the given refiner type.");
		}
	}
*/
}



// explicit template specializations
#ifdef UG_DIM_1
	template class SurfaceMarking<Domain1d>;
#endif
#ifdef UG_DIM_2
	template class SurfaceMarking<Domain2d>;
#endif
#ifdef UG_DIM_3
	template class SurfaceMarking<Domain3d>;
#endif



} // namespace neuro_collection
} // namespace ug
