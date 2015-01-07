/*
------------------------------------------------------------------------------------
General channel/pump transport interface class
------------------------------------------------------------------------------------
*/

#ifndef __IP3R_H__
#define __IP3R_H__

#include <iostream>
#include <cmath>
#include <string>


namespace ug{
namespace neuro_collection{



class IP3R : public IMembraneTransporter
{
	protected:
		const number R;			// universal gas constant
		const number T;			// temperature
		const number F;			// Faraday constant

		const number D1;		// IP3 binding (w/o Ca2+ inhibition)
		const number D2;		// Ca2+ inhibiting binding
		const number D3;		// IP3 binding (w/ Ca2+ inhibition)
		const number D5;		// Ca2+ activating binding
		const number RHO_IP3R;	// average channel density for IP3R channel
		const number MU_IP3R;	// IP3R channel conductance

		const number REF_CA_ER;	// reference endoplasmatic Ca2+ concentration (for conductances)

	public:
		// constructor
		IP3R() : IMembraneTransporter(),
		R(8.314), T(310.0), F(96485.0),
		D1(1.3e-7), D2(1.05e-6), D3(9.4e-7), D5(8.23e-8), RHO_IP3R(17.3), MU_IP3R(1.6e-12),
		REF_CA_ER(2.5e-4) {};

		// destructor
		virtual ~IP3R() {};

		// flux output functions
		virtual void flux(const std::vector<number>& u, std::vector<number>& flux)
		{
			number caCyt = u[0];	// cytosolic Ca2+ concentration
			number caER = u[1];		// ER Ca2+ concentration
			number ip3 = u[2];		// IP3 concentration

			// membrane current corresponding to diffusion pressure
			number current = R*T/(4*F*F) * MU_IP3R/REF_CA_ER * (caER - caCyt);

			// flux calculation
			number pOpen = pow(caCyt*ip3*D2 / ((caCyt*ip3 + ip3*D2 + D1*D2 + caCyt*D3) * (caCyt+D5)),3);

			flux[0] = pOpen * current;

			// dimensional correction: concentrations are mol/dm^3, but length unit is um
			flux[0] *= 1e15;
		}


		virtual void flux_derivative(const std::vector<number>& u, std::vector<std::vector<number> >& flux_derivs)
		{
			// get values of the unknowns in associated node
			number caCyt = u[0];
			number caER = u[1];
			number ip3 = u[2];

			// membrane potential
			number current = R*T/(4*F*F) * MU_IP3R/REF_CA_ER * (caER - caCyt);
			number dI_dCyt = - R*T/(4*F*F) * MU_IP3R/REF_CA_ER;
			number dI_dER = -dI_dCyt;

			// IP3R flux derivatives
			number schlonz1 = caCyt*ip3 + ip3*D2 + D1*D2 + caCyt*D3;
			number schlonz2 = schlonz1 * (caCyt+D5);
			number schlonz3 = (caCyt*ip3*D2) / schlonz2;
			number pOpen = pow(schlonz3,3);

			number dOpen_dCyt = 3.0*schlonz3*schlonz3*ip3*D2 * (1.0 - caCyt/schlonz2 * ( (ip3+D3)*(caCyt+D5) + schlonz1 )) / schlonz2;
			number dOpen_dIP3 = 3.0*schlonz3*schlonz3*caCyt*D2 * (1.0 - ip3/schlonz2 * (caCyt+D2)*(caCyt+D5)) / schlonz2;

			flux_derivs[0][0] = dOpen_dCyt * current + pOpen * dI_dCyt;
			flux_derivs[0][1] = pOpen * dI_dER;
			flux_derivs[0][2] = dOpen_dIP3 * current;

			// dimensional correction: concentrations are mol/dm^3, but length unit is um
			for (size_t i = 0; i < u.size(); i++)
				flux_derivs[0][i] *= 1e15;
		}


		// return number of unknowns this transport mechanism depends on
		virtual size_t n_dependencies()
		{
			return 3;
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
			return std::string("IP3R");
		};
};


} // namespace neuro_collection
} // namespace ug

#endif // __IP3R_H__

