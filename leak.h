/*
 *	Discretization for the leakage flux through a membrane
 *
 *  Created on: 20.12.2011
 *      Author: mbreit
 */

#ifndef __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__LEAK_H__
#define __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__LEAK_H__

#include "membrane_transporter_interface.h"
#include "lib_disc/spatial_disc/elem_disc/inner_boundary/inner_boundary.h"


namespace ug{
namespace neuro_collection{


///@addtogroup plugin_neuro_collection
///@{


/// Discretization for the leakage flux through a membrane
/**
 * This class implements the leakage flux through a membrane. It assumes a linear dependency
 * of the flowing substance from the difference between source and target concentrations,
 * i.e. j = c * ([source]-[target]). The constant c can be set using parent class method
 * IMembraneTransporter::set_density_function(). Although this constitutes a slight misuse of
 * terminology (as the constant does not in fact represent a physical density of some kind of
 * "leakage channels" in the membrane), doing so allows a straightforward implementation.
 *
 * Units used in the implementation of this channel:
 * Concentrations:                    mM (= mol/m^3)
 * Leakage constant:                  m/s
 * Output unit of the "flux" method:  mM
 * Resulting flux density:            mol/(m^2 s)
 */

class Leak : public IMembraneTransporter
{
	public:
		enum{_S_=0, _T_};

	public:
	/// @copydoc IMembraneTransporter::IMembraneTransporter()
	Leak(std::vector<std::string> fcts);

	/// @copydoc IMembraneTransporter::IMembraneTransporter()
	virtual ~Leak();

	/// @copydoc IMembraneTransporter::calc_flux()
	virtual void calc_flux(const std::vector<number>& u, std::vector<number>& flux) const;

	/// @copydoc IMembraneTransporter::calc_flux_deriv()
	virtual void calc_flux_deriv(const std::vector<number>& u, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const;

	/// @copydoc IMembraneTransporter::n_dependencies()
	virtual const size_t n_dependencies() const;

	/// @copydoc IMembraneTransporter::n_fluxes()
	virtual size_t n_fluxes() const;

	/// @copydoc IMembraneTransporter::flux_from_to()
	virtual const std::pair<size_t,size_t> flux_from_to(size_t flux_i) const;

	/// @copydoc IMembraneTransporter::name()
	virtual const std::string name() const;

	/// @copydoc IMembraneTransporter::check_constant_allowed()
	virtual void check_constant_allowed(const size_t i, const number val) const;

	/// @copydoc IMembraneTransporter::print_units()
	virtual void print_units() const;
};

///@}


} // namespace neuro_collection
} // namespace ug

#endif // __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__LEAK_H__

