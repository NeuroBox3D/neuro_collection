/*
 * surface_marking.h
 *
 *  Created on: 30.03.2016
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__NEURO_COLLECTION__SURFACE_MARKING_H_
#define UG__PLUGINS__NEURO_COLLECTION__SURFACE_MARKING_H_

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

#endif // UG__PLUGINS__NEURO_COLLECTION__SURFACE_MARKING_H_
