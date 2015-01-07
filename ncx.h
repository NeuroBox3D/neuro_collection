/*
 *	Discretization for the NCX pump in the PM
 *
 *  Created on: 20.12.2011
 *      Author: mbreit
 */

#ifndef __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__NCX_H__
#define __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__NCX_H__

#include <string>


namespace ug{
namespace neuro_collection{


/// Discretization for the NCX pump in the PM
/**
 * This class implements the InnerBoundaryElemDisc to provide flux densities
 * and their derivatives for the Graupner (2003) model of NCX pumps.
 *
 * \tparam	TDomain		domain
 */

class NCX : public IMembraneTransporter
{
	protected:
		const number KD_N;			// mol*dm^-3
		const number IMAX_N;		// mol*s^-1
		const number RHO_N;//1.0e02	// um^-2

	public:
		// constructor
		NCX() : IMembraneTransporter(),
		KD_N(1.8e-6), IMAX_N(2.5e-21), RHO_N(15.0) {};

		// destructor
		virtual ~NCX() {};

		// flux output functions
		virtual void flux(const std::vector<number>& u, std::vector<number>& flux)
		{
			number caCyt = u[0];	// cytosolic Ca2+ concentration

			number gatingFactor = caCyt / (KD_N + caCyt);

			flux[0] = - gatingFactor * IMAX_N;

			// dimensional correction: concentrations are mol/dm^3, but length unit is um
			flux[0] *= 1e15;
		}


		virtual void flux_derivative(const std::vector<number>& u, std::vector<std::vector<number> >& flux_derivs)
		{
			// get values of the unknowns in associated node
			number caCyt = u[0];

			number dGating_dCyt = KD_N / std::pow(KD_N + caCyt, 2);
			flux_derivs[0][0] = - dGating_dCyt * IMAX_N;

			// dimensional correction: concentrations are mol/dm^3, but length unit is um
			flux_derivs[0][0] *= 1e15;
		}


		// return number of unknowns this transport mechanism depends on
		virtual size_t n_dependencies()
		{
			return 1;
		}

		// return number of fluxes calculated by this machanism
		virtual size_t n_fluxes()
		{
			return 1;
		};

		// from where to where do the fluxes occur
		virtual std::pair<size_t,size_t> flux_from_to(size_t flux_i)
		{
			return std::make_pair<size_t, size_t>(-1,0);
		}

		virtual std::string name()
		{
			return std::string("NCX");
		};
};


} // namespace neuro_collection
} // namespace ug

#endif // __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__NCX_H__

