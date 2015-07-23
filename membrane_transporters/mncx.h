/*
 *	 Discretization for the mitochondrial NCX in the mitochondrial membrane
 *
 *  Created on: 30.06.2015
 *      Author: mstepnie
 */

/*
 * Complex NCX model from Paper: "A Biophysically Based Mathematical Model for the Kinetics of Mitochondrial
 * Na-Ca Antiporter" (using Model 1 and kinectic values of reference 3)
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
 * [Ca_cyt]  mM (= mol/m^3)
 * [Ca_out]  mM (= mol/m^3)
 *
 * Ca flux   mol/s
 */

class MNCX : public IMembraneTransporter
{
	public:
        enum{_CCYT_=0, _CEXT_, _NCYT_, _NEXT_};


    protected:


        const number F, RT, ca_valence, na_valence;

		double m_psi;
		double alpha;
		double beta;

		double k_a;
		double K_N_0;
		double K_C_0;
		double K_C_int;
		double K_C_ext;
		double K_N_int;
		double K_N_ext;
		double K_H1;
		double K_H2;
		double H_int;
		double H_ext;

		double D_H_ext;
		double D_H_int;

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
};

///@}


} // namespace neuro_collection
} // namespace ug

#endif // __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__MNCX_H__

