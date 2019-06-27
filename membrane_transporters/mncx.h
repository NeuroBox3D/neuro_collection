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

/**
 * Complex NCX model from Papers:
 *
 * "A Biophysically Based Mathematical Model for the Kinetics of Mitochondrial
 * Na-Ca Antiporter" (using Model 1 and kinectic values of reference 1), Pradhan et al. 2010
 *
 * and
 *
 * "Analysis of cardiac mitochondrial Na Ca exchanger kinetics with a biophysical model of mitochondrial Ca handling
 * suggests a 3-1 stoichiometry", Pradhan et al. 2008
 *
 */

#ifndef UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__MNCX_H
#define UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__MNCX_H

#include "membrane_transporter_interface.h"
#include "lib_disc/spatial_disc/elem_disc/inner_boundary/inner_boundary.h"


namespace ug {
namespace neuro_collection {


///@addtogroup plugin_neuro_collection
///@{


/// Discretization for the MNCX
/**
 * This class implements the InnerBoundaryElemDisc to provide flux densities
 * and their derivatives for the MNCX
 *
 * Units used in the implementation of this channel:
 * [Ca_cyt]  M (= mol/m^3)
 * [Ca_mit]  M (= mol/m^3)
 * [Na_cyt]  M (= mol/m^3)
 * [Na_mit]  M (= mol/m^3)
 *
 * Ca flux   mol/s
 */

/**
 * 	PARAMETER INFOS (s. "A Biophysically Based Mathematical Model for the Kinetics of Mitochondrial
 * 						 Na-Ca Antiporter" (using Model 1 and kinectic values of reference 3), Pradhan et al. 2010
 *
 * 						 and
 *
 *				   		 "Analysis of cardiac mitochondrial Na Ca exchanger kinetics with a biophysical model of mitochondrial Ca handling
 * 						 suggests a 3-1 stoichiometry", Pradhan et al. 2008)
 *
 *	Dissociation constants for Ca2+ and Mg2+ binding to the antiporter
 *	Pradhan 2010:
 *		K_C 	= 2.28e-9	     	// in Mol
 *		K_N 	= 9.14e-3; 			// in Mol
 *	Pradhan 2008
 *		K_C 	= 2.1e-6	     	// in Mol
 *		K_N 	= 8.2e-3; 			// in Mol
 *
 *	Ratio of potential difference
 *	between Na2+ or Ca2+ bound to the site of antiporter
 *	facing the external (internal) side of the IMM and Na2+ or
 *	Ca2+ in the bulk phase to the total membrane potential
 *		alpha 	= 0.0;				// unitless
 *
 *	Displacement of external (internal)
 *	Na2+ or Ca2+ from the coordinate of maximum potential barrier
 *		beta 	= 0.5				// unitless
 *
 *	Rate constant for limiting translocation of Ca2+ across the mitochondiral membrane
 *		k 		= 4.9;			 	// Reference Paucek & Jaburek in umol/mg/min, Model 1 Pradhan 2010
 *				= 0.081666			// Reference Paucek & Jaburek in umol/mg/s,   Model 1 Pradhan 2010
 *		k 		= 1.41e-3;		 	// Reference Paucek & Jaburek in umol/mg/s,   Model   Pradhan 2008
 *              = 0.0846 			// Reference Paucek & Jaburek in umol/mg/min  Model   Pradhan 2008
 *
 *	Mitochondrial membrane potential
 *		m_psi 	 					// in mV;
 */

class MNCX : public IMembraneTransporter
{
	public:
        enum{_CCYT_=0, _CMIT_, _NCYT_, _NMIT_};


    protected:

		const number F;		  	// Faraday's constant in kJ/mol/mV
		const number RT;	  	// universal gas constant times temperature (310K) in kJ/mol

		const number K_C;		// in Mol
		const number K_N;		// in Mol

		const number k;			// Reference Paucek & Jaburek in umol/mg/min

		number m_psi;			// in mV

        number m_mit_volume;  	// in um^3
		number m_mit_surface; 	// in um^2

    public:
		/// @copydoc IMembraneTransporter::IMembraneTransporter(const std::vector<std::string)
        MNCX(const std::vector<std::string>& fcts);

		/// @copydoc IMembraneTransporter::IMembraneTransporter()
        MNCX(const char* fcts);

		/// @copydoc IMembraneTransporter::IMembraneTransporter()
		virtual ~MNCX();

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

		/// Sets mitochondrial membrane potential
		void set_psi(number psi);

		/// Debug method
		number get_flux(number ca_cyt, number ca_mit, number na_cyt, number na_mit, number psi);
};

///@}


} // namespace neuro_collection
} // namespace ug

#endif // UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__MNCX_H

