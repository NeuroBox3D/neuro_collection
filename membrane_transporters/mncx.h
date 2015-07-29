/*
 *	 Discretization for the mitochondrial NCX in the mitochondrial membrane
 *
 *  Created on: 30.06.2015
 *      Author: mstepnie
 */

/*
 * Complex NCX model from Paper: "A Biophysically Based Mathematical Model for the Kinetics of Mitochondrial
 * Na-Ca Antiporter" (using Model 1 and kinectic values of reference 1), Pradhan et al. 2010
 *
 */

#ifndef __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__MNCX_H__
#define __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__MNCX_H__

#include "membrane_transporter_interface.h"
#include "lib_disc/spatial_disc/elem_disc/inner_boundary/inner_boundary.h"


namespace ug{
namespace neuro_collection{


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
 * 						 Na-Ca Antiporter" (using Model 1 and kinectic values of reference 3), Pradhan et al. 2010)
 *
 *	Dissociation constants for Ca2+, Mg2+ and phosphate binding to the uniporter
 *		K_C 	= 2.28e-3	     	// in Mol
 *		K_N 	= 9.14e-3; 			// in Mol
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
 *		k 		= 4.9;			 	// Reference Paucek & Jaburek in umol/mg/min
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
		virtual const size_t n_dependencies() const;

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

#endif // __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__MNCX_H__

