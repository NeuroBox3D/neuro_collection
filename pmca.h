/*
 *	 Discretization for the PMCA pump in the PM
 *
 *  Created on: 20.12.2011
 *      Author: mbreit
 */

#ifndef __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__PMCA_H__
#define __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__PMCA_H__

#include <string>


namespace ug{
namespace neuro_collection{


/// Discretization for the PMCA pump in the PM
/**
 * This class implements the InnerBoundaryElemDisc to provide flux densities
 * and their derivatives for the Graupner (2003) model of PMCA pumps (saturated CaM).
 *
 * \tparam	TDomain		domain
 */

class PMCA : public IMembraneTransporter
{
	protected:
		const number KD_P;					// mol*dm^-3 (Elwess et al.)
		//const number KD_P = 3.4e-07;		// mol*dm^-3 (Graupner)
		const number IMAX_P;				// mol*s^-1
		const number RHO_P;		//3.0e03;	// um^-2

	public:
		// constructor
		PMCA() : IMembraneTransporter(),
		KD_P(6.0e-8), IMAX_P(1.7e-23), RHO_P(500.0) {};

		// destructor
		virtual ~PMCA() {};

		// flux output functions
		virtual void flux(const std::vector<number>& u, std::vector<number>& flux)
		{
			number caCyt = u[0];	// cytosolic Ca2+ concentration

			number gatingFactor = caCyt*caCyt / (KD_P*KD_P + caCyt*caCyt);

			flux[0] = - gatingFactor * IMAX_P;

			// dimensional correction: concentrations are mol/dm^3, but length unit is um
			flux[0] *= 1e15;
		}


		virtual void flux_derivative(const std::vector<number>& u, std::vector<std::vector<number> >& flux_derivs)
		{
			// get values of the unknowns in associated node
			number caCyt = u[0];

			number dGating_dCyt = 2*KD_P*KD_P*caCyt / std::pow(KD_P*KD_P + caCyt*caCyt, 2);
			flux_derivs[0][0] = - dGating_dCyt * IMAX_P;

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
			return std::string("PMCA");
		};
};


} // namespace neuro_collection
} // namespace ug

#endif // __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__PMCA_H__

