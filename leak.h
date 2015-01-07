/*
 *	Discretization for the leakage flux through a membrane
 *
 *  Created on: 20.12.2011
 *      Author: mbreit
 */

#ifndef __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__LEAK_H__
#define __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__LEAK_H__

#include <string>


namespace ug{
namespace neuro_collection{


/// Discretization for the leakage flux through a membrane
/**
 * This class implements the leakage flux through a membrane.
 *
 * \tparam	TDomain		domain
 */

class Leak : public IMembraneTransporter
{
	public:
		// constructor
		Leak() : IMembraneTransporter() {};

		// destructor
		virtual ~Leak() {};

		// flux output functions
		virtual void flux(const std::vector<number>& u, std::vector<number>& flux)
		{
			number caCyt = u[0];	// cytosolic Ca2+ concentration
			number caER = u[1];		// ER Ca2+ concentration

			// membrane current corresponding to diffusion pressure
			// cheating a little here: utilizing density for leakage flux "density" constant,
			// since the leakage flux is not really caused by a density, but depends of course
			// on the densities of all the pumps/channels in the membrane and is therefore
			// position/subset dependent
			flux[0] = caER-caCyt;

			// dimensional correction: concentrations are mol/dm^3, but length unit is um
			flux[0] *= 1e15;
		}


		virtual void flux_derivative(const std::vector<number>& u, std::vector<std::vector<number> >& flux_derivs)
		{
			flux_derivs[0][0] = -1.0;
			flux_derivs[0][1] = 1.0;

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
			return std::string("Leak");
		};
};


} // namespace neuro_collection
} // namespace ug

#endif // __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__LEAK_H__

