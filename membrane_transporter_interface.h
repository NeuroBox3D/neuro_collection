/*
 *  General channel/pump transport interface class
 *
 *  Created on: 07.01.2015
 *     Authors: mbreit, mstepniewski
 */

#ifndef __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__MEMBRANE_TRANSPORTER_H__
#define __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__MEMBRANE_TRANSPORTER_H__


#include <string>


namespace ug{
namespace neuro_collection{


class IMembraneTransporter
{
	public:
		// constructor
		IMembraneTransporter() {};

		// destructor
		virtual ~IMembraneTransporter() {};

		// flux output functions
		virtual void flux(const std::vector<number>& u, std::vector<number>& flux) = 0;
		virtual void flux_derivative(const std::vector<number>& u, std::vector<std::vector<number> >& flux_derivs) = 0;

		// return number of unknowns this transport mechanism depends on
		virtual size_t n_dependencies() = 0;

		// return number of fluxes calculated by this machanism
		virtual size_t n_fluxes() = 0;

		// from where to where do the fluxes occur
		virtual std::pair<size_t,size_t> flux_from_to(size_t flux_i) = 0;

		virtual std::string name() = 0;
};


} // namespace neuro_collection
} // namespace ug

#endif // __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__MEMBRANE_TRANSPORTER_H__

