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

#ifndef UG__PLUGINS__NEURO_COLLECTION__SURFACE_MARKING_H
#define UG__PLUGINS__NEURO_COLLECTION__SURFACE_MARKING_H

#include <cstddef>                                                 // for size_t
#include <utility>                                                 // for pair
#include <vector>                                                  // for vector
#include "common/types.h"                                          // for number
#include "common/util/smart_pointer.h"                             // for ConstSmartPtr, SmartPtr
#include "lib_disc/function_spaces/error_elem_marking_strategy.h"  // for IElementMarkingStrategy


namespace ug {

// forward declarations
class DoFDistribution;
class IRefiner;
template <typename TDomain> class ApproximationSpace;

namespace neuro_collection {


/// mark elements adjacent to given surfaces
template <typename TDomain>
class SurfaceMarking : public IElementMarkingStrategy<TDomain>
{

public:
	typedef IElementMarkingStrategy<TDomain> base_type;

	SurfaceMarking(ConstSmartPtr<TDomain > dom)
	: m_tol(-1.0), m_max_level(-1), m_spDom(dom) {}

	SurfaceMarking(number tol, size_t max_level)
	: m_tol(tol), m_max_level(max_level), m_spDom(SPNULL) {};

	SurfaceMarking(ConstSmartPtr<TDomain > dom, number tol, size_t max_level)
	: m_tol(tol), m_max_level(max_level), m_spDom(dom) {}

	void set_tolerance(number tol) {m_tol = tol;}
	void set_max_level(size_t max_level) {m_max_level = max_level;}

	void add_surface(int surf_si, int adj_vol_si);
	void remove_surface(int surf_si, int adj_vol_si);

	void add_surface(const std::string& surf_si, const std::string& adj_vol_si);
	void remove_surface(const std::string& surf_si, const std::string& adj_vol_si);

	void mark
	(
		typename base_type::elem_accessor_type& aaError,
		IRefiner& refiner,
		ConstSmartPtr<DoFDistribution> dd
	);

	void mark_without_error(SmartPtr<IRefiner> refiner, SmartPtr<ApproximationSpace<TDomain> > approx);

protected:
	void mark(IRefiner& refiner, ConstSmartPtr<DoFDistribution> dd);

protected:
	number m_tol;
	size_t m_max_level;

	ConstSmartPtr<TDomain> m_spDom;

	/// vector holding pairs of surface and adjacent
	/// element subset indices
	std::vector<std::pair<int, int> > m_vpSurfaces;
};


} // namespace neuro_collection
} // namespace ug

#endif // UG__PLUGINS__NEURO_COLLECTION__SURFACE_MARKING_H
