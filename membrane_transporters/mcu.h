/*
 *	 Discretization for the mitochondrial uniporter MCU in the mitochondrial membrane
 *
 *  Created on: 30.06.2015
 *      Author: mstepnie
 */

#ifndef __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__MCU_H__
#define __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__MCU_H__

#include "membrane_transporter_interface.h"
#include "lib_disc/spatial_disc/elem_disc/inner_boundary/inner_boundary.h"


namespace ug{
namespace neuro_collection{


///@addtogroup plugin_neuro_collection
///@{


/// Discretization for the MCU
/**
 * This class implements the InnerBoundaryElemDisc to provide flux densities
 * and their derivatives for the MCU
 *
 * Units used in the implementation of this channel:
 * [Ca_cyt]  mM (= mol/m^3)
 * [Ca_out]  mM (= mol/m^3)
 *
 * Ca flux   mol/s
 */

class MCU : public IMembraneTransporter
{
	public:
        enum{_CCYT_=0, _CEXT_};


    protected:

// 	Todo: 	initialize const variables in constructor explicitly!

        const number F, RT, ca_valence;

        number K_C;     // in mMol
        number K_M; 	  // in mMol
        number gamma;
        number k ;	  // Reference  (9) Crompton et al. in mmol/mg/s
        					  // Reference (14) Vinogadrov and Scarpa in mmol/mg/s
        number nH;
        number K_Pi;    // in mMol
        number K_CC;
        number K_MM;
        number m_psi;

        number mg_int, mg_ext;

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
};

///@}


} // namespace neuro_collection
} // namespace ug

#endif // __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__MCU_H__

