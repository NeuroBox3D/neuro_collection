/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Martin Stepniewski
 * Creation date: 2015-06-30
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

#ifndef UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__MCU_H
#define UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__MCU_H

#include "membrane_transporter_interface.h"
#include "lib_disc/spatial_disc/elem_disc/inner_boundary/inner_boundary.h"


namespace ug{
namespace neuro_collection{


///@addtogroup plugin_neuro_collection
///@{


/// Discretization for the MCU
/**
 * This class implements the InnerBoundaryElemDisc to provide flux densities
 * and their derivatives for the mitochondrial uniporter MCU
 * in the mitochondrial membrane.
 *
 * Units used in the implementation of this channel:
 * [Ca_cyt]  M (= mol/m^3)
 * [Ca_mit]  M (= mol/m^3)
 *
 * Ca flux   mol/s
 */

/**
 * 	PARAMETER INFOS (s. "Characterization of Mg2+ inhibition of MCU", Pradhan et al. 2011)
 *
 *	Dissociation constants for Ca2+, Mg2+ and phosphate binding to the uniporter
 *		K_C 	= 3.965e-6;     	// in Mol
 *		K_M 	= 0.655e-3; 		// in Mol
 *		K_Pi 	= 0.2e-3 * 1e3;     // in Mol
 *
 *	Phophate effect on the uniporter function @ dissociation constants
 * 		K_CC 	= K_C * (1 + pi_cyt/(K_Pi+pi_cyt));
 *		K_MM 	= K_M / (1 + pi_cyt/(K_Pi+pi_cyt));
 *
 *	Binding affinity for Ca2+ and Mg2+ binding to the uniporter
 *		g(amma) = 3.9;				// unitless
 *
 *	Rate constant for limiting translocation of Ca2+ across the mitochondiral membrane
 *	    k 		= 27.94e-3;			// Reference (15) Scarpa and Graziotti in nmol/mg/s
 *		k 		= 16.51e-3;			// Reference (14) Vinogadrov and Scarpa in nmol/mg/s
 *		k 		= 1.27e-3;			// Reference  (9) Crompton et al. in nmol/mg/s
 *
 *	Deviation from linear GHK type formalism
 *		nH 		= 2.65;				// unitless
 *
 *  Cytosolic phosphate concentration
 * 		m_pi_cyt					// in Mol
 *
 *  Cytosolic Mg2+ concentration
 * 		m_mg_cyt					// in Mol
 *
 *  Mitochondrial Mg2+ concentration
 * 		m_mg_mit					// in Mol
 *
 *	Mitochondrial membrane potential
 *		m_psi 	 					//  in mV;
 */

class MCU : public IMembraneTransporter
{
	public:
        enum{_CCYT_=0, _CMIT_};


    protected:

		const number F;		  // Faraday's constant in kJ/mol/mV
		const number RT;	  // universal gas constant times temperature (310K) in kJ/mol

        const number K_C;     // in Mol
        const number K_M; 	  // in Mol
        const number K_Pi;    // in Mol

        const number g;	  	  // unitless
        const number nH;	  // unitless

        number m_k ;  		  // Reference  (9) Crompton et al. in nmol/mg/s
        					  // Reference (14) Vinogadrov and Scarpa in nmol/mg/s
        					  // Reference (15) Scarpa and Graziotti in nmol/mg/s

        number m_pi_cyt;	  // Phosphate concentration in Mol
        number m_mg_cyt;	  // Mg2+ concentration in Mol
        number m_mg_mit;	  // Mg2+ concentration in Mol

        number K_CC;		  // in Mol
        number K_MM;          // in Mol
        number m_psi;	      // in mV

        number m_mit_volume;  // in um^3
		number m_mit_surface; // in um^2



    public:
		/// @copydoc IMembraneTransporter::IMembraneTransporter(const std::vector<std::string)
		MCU(const std::vector<std::string>& fcts);

		/// @copydoc IMembraneTransporter::IMembraneTransporter()
		MCU(const char* fcts);

		/// @copydoc IMembraneTransporter::IMembraneTransporter()
		virtual ~MCU();

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


		// helper methods

		///	Sets mitochondrial volume
		void set_mit_volume(number mit_volume);

		/// Sets mitochondrial surface
		void set_mit_surface(number mit_surface);

		/// Sets cytosolic phosphate concentration
		void set_pi_cyt(number pi_cyt);

		/// Sets mitochondrial membrane potential
		void set_psi(number psi);

		/// Sets cytosolic Mg2+ concentration
		void set_mg_cyt(number mg_cyt);

		/// Sets mitochondrial Mg2+ concentration
		void set_mg_mit(number mg_mit);

		/// Sets rate constant k
		void set_rate_constant(number k);

		/// Debug method
		number get_flux(number ca_cyt, number ca_mit, number pi, number mg_cyt, number mg_mit, number psi);
};

///@}


} // namespace neuro_collection
} // namespace ug

#endif // UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__MCU_H

