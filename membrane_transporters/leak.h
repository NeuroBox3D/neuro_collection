/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2011-12-20
 *
 * This file is part of NeuroBox, which is based on UG4.
 *
 * NeuroBox and UG4 are free software: You can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3
 * (as published by the Free Software Foundation) with the following additional
 * attribution requirements (according to LGPL/GPL v3 §7):
 *
 * (1) The following notice must be displayed in the appropriate legal notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 *
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 *
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating PDE based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * "Stepniewski, M., Breit, M., Hoffer, M. and Queisser, G.
 *   NeuroBox: computational mathematics in multiscale neuroscience.
 *   Computing and visualization in science (2019).
 * "Breit, M. et al. Anatomically detailed and large-scale simulations studying
 *   synapse loss and synchrony using NeuroBox. Front. Neuroanat. 10 (2016), 8"
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__LEAK_H
#define UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__LEAK_H

#include "membrane_transporter_interface.h"
#include "lib_disc/spatial_disc/elem_disc/inner_boundary/inner_boundary.h"


namespace ug{
namespace neuro_collection{


///@addtogroup plugin_neuro_collection
///@{


/// Discretization for the leakage flux through a membrane
/**
 * This class implements a leakage flux through a membrane. It assumes a current density
 * directed from source to target of Goldman-Hodgkin-Katz type, i.e.,
 *
 *     j = - p * zF/(RT) * V * (c_s - c_t * exp(zF/(RT)*V)) / (1 - exp(zF/(RT)*V)),
 *
 * with a permeability p, the source and target ionic concentrations c_s and c_t
 * and the membrane voltage V which is defined as V = (Phi_t - Phi_s) (which is not
 * necessarily the same way it is usually defined on the plasma membrane).
 *
 * Both the source and target concentrations as well as the source and target potentials
 * are functions that can be supplied.
 * If the potentials are not given, they are assumed to be zero, so the expression
 * for the current density is reduced to
 *
 *     j = p * (c_s - c_t).
 *
 * The constant p can either be set using the method
 *     set_permeability(),
 * then the output of the calc_flux() method is a molar flux density (units: mol/(m^2 s))
 * and the density of the mechanism needs to be set to 1 in MembraneTransportFV1
 * OR by leaving the permeability at its default value of 1 and instead setting the
 * density to its value using the method
 *     MembraneTransportFV1::set_density_function(),
 * whatever the user prefers.
 * Both ways constitute a slight misuse of terminology. :)
 *
 * Units used in the implementation of this channel:
 * Concentrations:                    mM (= mol/m^3)
 * Potentials:                        V
 * Leakage constant:                  m/s
 * Output unit of the "flux" method:  mol/(m^2 s) or mM
 * Resulting flux density:            mol/(m^2 s)
 */

class Leak : public IMembraneTransporter
{
	public:
		enum{_S_ = 0, _T_, _PHIS_, _PHIT_};

	public:
		/// @copydoc IMembraneTransporter::IMembraneTransporter(const std::vector<std::string)
		Leak(const std::vector<std::string>& fcts);

		/// @copydoc IMembraneTransporter::IMembraneTransporter()
		Leak(const char* fcts);

		/// @copydoc IMembraneTransporter::IMembraneTransporter()
		virtual ~Leak();

		void set_permeability(number p);

		void set_temperature(number t);

		void set_valency(int v);

		/// @copydoc IMembraneTransporter::calc_flux()
		virtual void calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const;

		/// @copydoc IMembraneTransporter::calc_flux_deriv()
		virtual void calc_flux_deriv(const std::vector<number>& u, GridObject* e, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const;

		/// @copydoc IMembraneTransporter::n_dependencies()
		virtual size_t n_dependencies() const;

		/// @copydoc IMembraneTransporter::n_fluxes()
		virtual size_t n_fluxes() const;

		/// @copydoc IMembraneTransporter::flux_from_to()
		virtual const std::pair<size_t,size_t> flux_from_to(size_t flux_i) const;

		/// @copydoc IMembraneTransporter::name()
		virtual const std::string name() const;

		/// @copydoc IMembraneTransporter::check_supplied_functions()
		virtual void check_supplied_functions() const;

		/// @copydoc IMembraneTransporter::print_units()
		virtual void print_units() const;

	protected:
		number m_perm;
		number m_temp;
		int m_z;
		bool m_bNoVoltage;
};

///@}


} // namespace neuro_collection
} // namespace ug

#endif // UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__LEAK_H

