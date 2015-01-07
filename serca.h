/*
 *	Discretization for the SERCA calcium pump in the ER membrane
 *
 *  Created on: 20.12.2011
 *      Author: mbreit
 */

#ifndef __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__SERCA_H__
#define __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__SERCA_H__

#include <string>


namespace ug{
namespace neuro_collection{


/// Discretization for the SERCA calcium pump in the ER membrane
/**
 * This class implements the InnerBoundaryElemDisc to provide flux densities
 * and their derivatives for the Sneyd et al. (2003) SERCA model.
 *
 * \tparam	TDomain		domain
 */

class SERCA : public IMembraneTransporter
{
	protected:
		const number VS;			// maxmimal transport power per pump (divided by [ca_ER])
		const number KS;			// concentration at which halfmaximal pumping occurs
		const number RHO_SERCA;	// average pump density (fit to zero-flux at resting)

	public:
		// constructor
		SERCA() : IMembraneTransporter(),
		VS(6.5e-27), KS(1.8e-7), RHO_SERCA(1973.0) {};

		// destructor
		virtual ~SERCA() {};

		// flux output functions
		virtual void flux(const std::vector<number>& u, std::vector<number>& flux)
		{
			// get values of the unknowns in associated node
			number caCyt = u[0];	// cytosolic Ca2+ concentration
			number caER = u[1];		// ER Ca2+ concentration

			flux[0] = - VS*caCyt / ((KS + caCyt) * caER);

			// dimensional correction: concentrations are mol/dm^3, but length unit is um
			flux[0] *= 1e15;
		}


		virtual void flux_derivative(const std::vector<number>& u, std::vector<std::vector<number> >& flux_derivs)
		{
			// get values of the unknowns in associated node
			number caCyt = u[0];
			number caER = u[1];

			// SERCA flux derivs
			flux_derivs[0][0] = - (VS*KS) / (pow(KS+caCyt,2)*caER);
			flux_derivs[0][1] = VS*caCyt / ((KS+caCyt)*caER*caER);

			// dimensional correction: concentrations are mol/dm^3, but length unit is um
			for (size_t i = 0; i < u.size(); i++)
				flux_derivs[0][i] *= 1e15;
		}


		// return number of unknowns this transport mechanism depends on
		virtual size_t n_dependencies()
		{
			return 2;
		}

		// return number of fluxes calculated by this machanism
		virtual size_t n_fluxes()
		{
			return 1;
		};

		// from where to where do the fluxes occur
		virtual std::pair<size_t,size_t> flux_from_to(size_t flux_i)
		{
			return std::make_pair<size_t, size_t>(1,0);
		}

		virtual std::string name()
		{
			return std::string("SERCA");
		};
};


} // namespace neuro_collection
} // namespace ug

#endif // __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__SERCA_H__

