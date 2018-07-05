/*
 * ca_wave_util.h
 *
 *  Created on: 2017-09-18
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__NEURO_COLLECTION__UTIL__CA_WAVE_UTIL_H
#define UG__PLUGINS__NEURO_COLLECTION__UTIL__CA_WAVE_UTIL_H

#include "common/types.h"  // number
#include "common/util/smart_pointer.h"
#include "lib_disc/function_spaces/grid_function.h"

#include "../membrane_transporters/ryr_implicit.h"


namespace ug {
namespace neuro_collection {

///@addtogroup plugin_neuro_collection
///@{

/**
 * @brief Find maximal RyR flux density through a membrane subset.
 *
 * This function assumes that all of the following functions are
 * defined on the membrane subset:
 * - c_c, c_e (cytosolic and ER calcium concentrations)
 * - c1, c2 (RyR state variables from the Keizer & Levine model)
 * Their names must be provided in the fctNames argument in this order.
 *
 * The function is designed to work with the fully implicit RyR
 * implementation.
 *
 * @param u             solution grid function
 * @param fctNames      names for functions ca_cyt, ca_er, c1, c2 (in this order)
 * @param subsetNames   names of all ER membrane subsets with RyR channels
 * @param ryr           RyR channel object
 * @return              maximal flux density through RyR channel (mol/(m^2*s))
 */
template <typename TGridFunction>
number maxRyRFluxDensity
(
	ConstSmartPtr<TGridFunction> u,
	const char* fctNames,
	const char* subsetNames,
	ConstSmartPtr<RyRImplicit<typename TGridFunction::domain_type> > ryr
);


/**
 * @brief Find rightmost point of RyR activity.
 *
 * This function is useful in determining the wave front position
 * of a RyR-induced calcium wave from the left to the right.
 * The wave front is determined by a user-specified threshold on the
 * open probability of the RyR channel.
 *
 * The function is designed to work with the fully implicit RyR
 * implementation.
 *
 * @param u             solution grid function
 * @param fctNames      names for functions c1, c2 (in this order)
 * @param subsetName    name of ER membrane subset with RyR channels
 * @param thresh        threshold value of (1-(c1+c2))
 * @return              rightmost vertex where threshold value is exceeded
 */
template <typename TGridFunction>
number waveFrontX
(
	SmartPtr<TGridFunction> u,
	const char* fctNames,
	const char* subsetName,
	number thresh
);

/**
 * @brief Export solution on 1d subset, sorted by x coordinate.
 *
 * Processors must hold contiguous (in direction of x) parts of the subsets
 * being processed (otherwise, the export order will be screwed up).
 *
 * @param u              solution grid function
 * @param fctNames       name for function profiles to be exported
 * @param subsetNames    names of subsets to export profiles for
 * @param fileBaseName   base name of output files
 * @param time           point in time
 */
template <typename TDomain, typename TAlgebra>
class WaveProfileExporter
{
	public:
		typedef GridFunction<TDomain, TAlgebra> gf_type;

	public:
		WaveProfileExporter
		(
			SmartPtr<ApproximationSpace<TDomain> > approxSpace,
			const char* fctNames,
			const char* subsetNames,
			const std::string& fileBaseName
		);

		void exportWaveProfileX(ConstSmartPtr<gf_type> u, number time);

	private:
		struct CmpVrtPos
		{
			public:
				CmpVrtPos(const std::vector<number>& vPos)
				: pvPos(vPos) {}

				bool operator()(const size_t& a, const size_t& b)
				{
					return pvPos[a] < pvPos[b];
				}

			private:
				const std::vector<number>& pvPos;
		};

	private:
		SmartPtr<ApproximationSpace<TDomain> > m_spApprox;
		std::vector<std::string> m_vFct;
		std::vector<std::string> m_vSs;
		const std::string m_fileName;

		std::vector<std::vector<std::vector<DoFIndex> > > m_vvvDoFSeries;
		std::vector<std::vector<number> > m_vvXPos;
};


///@}

} // namespace ug
} // namespace neuro_collection

#include "ca_wave_util_impl.h"

#endif // UG__PLUGINS__NEURO_COLLECTION__UTIL__CA_WAVE_UTIL_H
