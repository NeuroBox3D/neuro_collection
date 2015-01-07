/*
 *	Discretization for the RyR calcium channel in the ER membrane
 *
 *  Created on: 20.12.2011
 *      Author: mbreit
 */

#ifndef __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__RYR_H__
#define __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__RYR_H__

#include <string>


namespace ug{
namespace neuro_collection{


/// Discretization for the RyR calcium channel in the ER membrane
/**
 * This class implements the InnerBoundaryElemDisc to provide flux densities
 * and their derivatives for the Keizer & Levine (1996) RyR model.
 *
 * \tparam	TDomain		domain
 */

class RyR : public IMembraneTransporter
{
	protected:
		const number R;			// universal gas constant
		const number T;			// temperature
		const number F;			// Faraday constant

		const number KA;		// calcium binding (C1 <--> O1)
		const number KB;		// calcium binding (O1 <--> O2)
		const number KC;		// O1 <--> C2
		const number RHO_RYR;	// average channel density for RyR channel
		const number MU_RYR;	// RyR channel conductance

		const number REF_CA_ER;	// reference endoplasmatic Ca2+ concentration (for conductances)

	public:
		// constructor
		RyR() : IMembraneTransporter(),
		R(8.314), T(310.0), F(96485.0),
		KA(5.21e25), KB(3.89e18), KC(17.5), RHO_RYR(0.86), MU_RYR(5.0e-11),
		REF_CA_ER(2.5e-4) {};

		// destructor
		virtual ~RyR() {};

		// flux output functions
		virtual void flux(const std::vector<number>& u, std::vector<number>& flux)
		{
			number caCyt = u[0];	// cytosolic Ca2+ concentration
			number caER = u[1];		// ER Ca2+ concentration

			// RyR flux calculation
			number current = R*T/(4*F*F) * MU_RYR/REF_CA_ER * (caER - caCyt);
			number pOpen =  (1.0 + KB*pow(caCyt,3)) / (1.0 + KC + 1.0/(KA*pow(caCyt,4)) + KB*pow(caCyt,3));

			flux[0] = pOpen * current;

			// dimensional correction: concentrations are mol/dm^3, but length unit is um
			flux[0] *= 1e15;
		}


		virtual void flux_derivative(const std::vector<number>& u, std::vector<std::vector<number> >& flux_derivs)
		{
			// get values of the unknowns in associated node
			number caCyt = u[0];
			number caER = u[1];

			// membrane potential
			number current = R*T/(4*F*F) * MU_RYR/REF_CA_ER * (caER - caCyt);
			number dI_dCyt = - R*T/(4*F*F) * MU_RYR/REF_CA_ER;
			number dI_dER = -dI_dCyt;

			// RyR flux derivs
			number schlonz1 = 1.0 + KB*pow(caCyt,3);
			number schlonz2 = 1.0 + KC + 1.0/(KA*pow(caCyt,4)) + KB*pow(caCyt,3);
			number pOpen = schlonz1 / schlonz2;

			number dOpen_dCyt = (3.0*KB*caCyt*caCyt + schlonz1/schlonz2*(4.0/(KA*pow(caCyt,5)) - 3.0*KB*caCyt*caCyt)) / schlonz2;

			flux_derivs[0][0] = dOpen_dCyt * current + pOpen * dI_dCyt;
			flux_derivs[0][1] = pOpen * dI_dER;

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
			return std::string("RyR");
		};
};


} // namespace neuro_collection
} // namespace ug

#endif // __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__RYR_H__

