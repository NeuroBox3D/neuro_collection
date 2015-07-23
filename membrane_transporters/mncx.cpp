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

#include "mncx.h"

namespace ug{
namespace neuro_collection{


MNCX::MNCX(const std::vector<std::string>& fcts) : IMembraneTransporter(fcts),
		   F(0.096484), RT(2.5775), ca_valence(2.0), na_valence(1.0)
{
// TODO some usefull setters
	m_psi 	= 0.0; //psi;
	alpha 	= 0.0;
	beta 	= 0.5;

	//k_a 	= 4.9;   		// in nmol/mg/min
	//k_a 	= 6.8;   		// in nmol/mg/min
	k_a = 0.0; //k;

	//K_N_0 	= 2.2e-3;  		// in nMol
	//K_N_0 	= 9.14e-3;  		// in nMol
	K_N_0 = 0.0; //K_N;

	//K_C_0 	= 2.28e-9;    	// in nMol
	K_C_0 = 0.0; //K_C;

	K_H1 	= 6.45e1;  		// in nMol
	K_H2 	= 1.39e2;  		// in nMol

//  pH-Werte
	H_ext = 7.3e-10; // pH von 6
	H_int = 7.3e-10; // pH von 6

//  pH Abhängigkeit verändert Konstante
//  berechnet die pH-Abhängigkeitsfaktoren (da kein dynamischer pH erwünscht im Konstruktor)
	/*
	double tmpFactor_int = (H_int/K_H1 + 1 + K_H2/H_int);
	double tmpFactor_ext = (H_ext/K_H1 + 1 + K_H2/H_ext);

	D_H_int = tmpFactor_int * (H_int/K_H2);
	D_H_ext = tmpFactor_ext * (H_ext/K_H2);

	K_C_int = tmpFactor_int * K_C_0;
	K_C_ext = tmpFactor_ext * K_C_0;
	*/

	K_C_int = K_C_0;
	K_C_ext = K_C_0;

	K_N_int = K_N_0;
	K_N_ext = K_N_0;

	D_H_int = 1.0;
	D_H_ext = 1.0;
}

MNCX::MNCX(const char* fcts) : IMembraneTransporter(fcts),
		   F(0.096484), RT(2.5775), ca_valence(2.0), na_valence(1.0)
{
//	IMPLEMENT INITIALIZATION

	m_psi 	= 0.0; //psi;
	alpha 	= 0.0;
	beta 	= 0.5;

	//k_a 	= 4.9;   		// in nmol/mg/min
	//k_a 	= 6.8;   		// in nmol/mg/min
	k_a = 0.0; //k;

	//K_N_0 	= 2.2e-3;  		// in nMol
	//K_N_0 	= 9.14e-3;  		// in nMol
	K_N_0 = 0.0; //K_N;

	//K_C_0 	= 2.28e-9;    	// in nMol
	K_C_0 = 0.0; //K_C;

	K_H1 	= 6.45e1;  		// in nMol
	K_H2 	= 1.39e2;  		// in nMol

//  pH-Werte
	H_ext = 7.3e-10; // pH von 6
	H_int = 7.3e-10; // pH von 6

//  pH Abhängigkeit verändert Konstante
//  berechnet die pH-Abhängigkeitsfaktoren (da kein dynamischer pH erwünscht im Konstruktor)
	/*
	double tmpFactor_int = (H_int/K_H1 + 1 + K_H2/H_int);
	double tmpFactor_ext = (H_ext/K_H1 + 1 + K_H2/H_ext);

	D_H_int = tmpFactor_int * (H_int/K_H2);
	D_H_ext = tmpFactor_ext * (H_ext/K_H2);

	K_C_int = tmpFactor_int * K_C_0;
	K_C_ext = tmpFactor_ext * K_C_0;
	*/

	K_C_int = K_C_0;
	K_C_ext = K_C_0;

	K_N_int = K_N_0;
	K_N_ext = K_N_0;

	D_H_int = 1.0;
	D_H_ext = 1.0;
}


MNCX::~MNCX()
{
	// nothing to do
}


void MNCX::calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const
{
// 	get values of the unknowns in associated node
	number caCyt = u[_CCYT_];	// cytosolic Ca2+ concentration
	number caExt = u[_CEXT_];	// extracellular Ca2+ concentration

	number naCyt = u[_NCYT_];	// cytosolic Na+ concentration
	number naExt = u[_NEXT_];	// extracellular Na2+ concentration

	double phi = F * m_psi / RT;

	double KNe = K_N_ext * exp(-alpha*na_valence*phi);
	double KNx = K_N_int * exp( alpha*na_valence*phi);
	double KCe = K_C_ext * exp(-alpha*ca_valence*phi);
	double KCx = K_C_int * exp( alpha*ca_valence*phi);

	//double k_in  = k_a * exp(( 3*beta*na_valence - beta*ca_valence)*phi);
	//double k_out = k_a * exp((-3*beta*na_valence + beta*ca_valence)*phi);
	double k_in  = k_a * exp(0.5*phi);
	double k_out = k_a * exp(-0.5*phi);

	double in_term  = 1 + pow(naExt, 3)/pow(KNe, 3) + caCyt/KCx + caCyt*pow(naExt, 3)/(KCx*(pow(KNe, 3)));
	double out_term = 1 + pow(naCyt, 3)/pow(KNx, 3) + caExt/KCe + caExt*pow(naCyt, 3)/(KCe*(pow(KNx, 3)));
	//double out_term = pow(naCyt, 3)/pow(KNx, 3) + caExt/KCe + caExt*pow(naCyt, 3)/(KCe*(pow(KNx, 3)));
	double D = D_H_int*in_term + D_H_ext*out_term - 1;
	//double D = in_term + out_term;

	double flux_in =  caCyt*pow(naExt, 3)/(KCx*pow(KNe, 3)) * k_in * D_H_int;
	double flux_out = caExt*pow(naCyt, 3)/(KCe*pow(KNx, 3)) * k_out * D_H_ext;

	//double flux_in =  caCyt*pow(naExt, 3)/(KCx*pow(KNe, 3));
	//double flux_out = caExt*pow(naCyt, 3)/(KCe*pow(KNx, 3));

	//double flux = 5.28/D*(flux_in - flux_out);
	double fluxes = 1/D*(flux_in - flux_out);

//  for flux in nmol per second and not minutes!
	//flux *= (1/60.0);

	// TODO right flux unit
	flux[0] = fluxes;
	flux[1] = fluxes*2;
}


void MNCX::calc_flux_deriv(const std::vector<number>& u, GridObject* e, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const
{
	// 	get values of the unknowns in associated node
		number caCyt = u[_CCYT_];	// cytosolic Ca2+ concentration
		number caExt = u[_CEXT_];	// extracellular Ca2+ concentration

		number naCyt = u[_NCYT_];	// cytosolic Na+ concentration
		number naExt = u[_NEXT_];	// extracellular Na2+ concentration

	// Setting vars for derivs
		double phi = F * m_psi / RT;

		double KNe = K_N_ext * exp(-alpha*na_valence*phi);
		double KNx = K_N_int * exp( alpha*na_valence*phi);
		double KCe = K_C_ext * exp(-alpha*ca_valence*phi);
		double KCx = K_C_int * exp( alpha*ca_valence*phi);

		double k_in  = k_a * exp(( 3*beta*na_valence - beta*ca_valence)*phi);
		double k_out = k_a * exp((-3*beta*na_valence + beta*ca_valence)*phi);

		double in_term  = 1 + pow(naExt, 3)/pow(KNe, 3) + caCyt/KCx + caCyt*pow(naExt, 3)/(KCx*(pow(KNe, 3)));
		double out_term = 1 + pow(naCyt, 3)/pow(KNx, 3) + caExt/KCe + caExt*pow(naCyt, 3)/(KCe*(pow(KNx, 3)));

		double D = D_H_int*in_term + D_H_ext*out_term - 1;

		double flux_in =  caCyt*pow(naExt, 3)/(KCx*pow(KNe, 3)) * k_in * D_H_int;
		double flux_out = caExt*pow(naCyt, 3)/(KCe*pow(KNx, 3)) * k_out * D_H_ext;


	// calculation really derivs of caInt
		double dD_caCyt = D_H_int/KCx + D_H_int*pow(naExt, 3)/(KCx*pow(KNe, 3));

		double dFlux_caCyt = k_in*D_H_int*pow(naExt, 3)/(KCx*pow(KNe, 3))/D - ((flux_in-flux_out)*dD_caCyt)/(D*D);

	//  for flux in nmol per second and not minutes!
		dFlux_caCyt *= (1/60.0);

		flux_derivs[0][0].first = local_fct_index(_CCYT_);
		flux_derivs[0][0].second = dFlux_caCyt; //TODO right flux unit



	// calculation of derivs of caExt
		double dD_caExt = D_H_ext/KCe + D_H_ext*pow(naCyt, 3)/(KCe*pow(KNx, 3));

		double dFlux_caExt = -k_out*D_H_ext*pow(naCyt, 3)/(KCe*pow(KNx, 3))/D - ((flux_in-flux_out)*dD_caExt)/(D*D);

	//  for flux in nmol per second and not minutes!
		dFlux_caExt *= (1/60.0);

		flux_derivs[0][1].first = local_fct_index(_CEXT_);
		flux_derivs[0][1].second = dFlux_caExt; //TODO right flux unit



	// calculation of derivs of naInt
		double dD_naCyt = D_H_int/pow(KNe, 3) + D_H_int*3*pow(naExt, 2)*caCyt/(KCx*pow(KNe, 3));

		double dFlux_naCyt = k_in*D_H_int*caCyt*3*pow(naExt, 2)/(KCx*pow(KNe, 3))/D - ((flux_in-flux_out)*dD_naCyt)/(D*D);

	//  for flux in nmol per second and not minutes!
		dFlux_naCyt *= (1/60.0);

		flux_derivs[0][2].first = local_fct_index(_NCYT_);
		flux_derivs[0][2].second = dFlux_naCyt; //TODO right flux unit



	// calculation of derivs of naExt

		double dD_naExt = D_H_ext/pow(KNx, 3) + D_H_ext*3*pow(naCyt, 2)*caExt/(KCe*pow(KNx, 3));

		double dFlux_naExt = -k_out*D_H_ext*caExt*3*pow(naCyt, 2)/(KCe*pow(KNx, 3))/D - ((flux_in-flux_out)*dD_naExt)/(D*D);

	//  for flux in nmol per second and not minutes!
		dFlux_naExt *= (1/60.0);

		flux_derivs[0][3].first = local_fct_index(_NEXT_);
		flux_derivs[0][3].second = dFlux_naExt; //TODO right flux unit


}


const size_t MNCX::n_dependencies() const
{
	return 2;
}


size_t MNCX::n_fluxes() const
{
	return 2;
};


const std::pair<size_t,size_t> MNCX::flux_from_to(size_t flux_i) const
{
	size_t from, to;
	if (allows_flux(_CCYT_)) from = local_fct_index(_CCYT_); else from = InnerBoundaryConstants::_IGNORE_;
	if (allows_flux(_CEXT_)) to = local_fct_index(_CEXT_); else to = InnerBoundaryConstants::_IGNORE_;

	if (allows_flux(_NCYT_)) from = local_fct_index(_NCYT_); else from = InnerBoundaryConstants::_IGNORE_;
	if (allows_flux(_NEXT_)) to = local_fct_index(_NEXT_); else to = InnerBoundaryConstants::_IGNORE_;

	return std::pair<size_t, size_t>(from, to);
}


const std::string MNCX::name() const
{
	return std::string("MNCX");
};


void MNCX::check_supplied_functions() const
{
	// Check that not both, inner and outer calcium concentrations are not supplied;
	// in that case, calculation of a flux would be of no consequence.
	if (!allows_flux(_CCYT_) && !allows_flux(_CEXT_))
	{
		UG_THROW("Supplying neither cytosolic nor extracellular calcium concentrations is not allowed.\n"
			 	 "This would mean that the flux calculation would be of no consequence\n"
				 "and this pump mechanism would not do anything.");
	}

	if (!allows_flux(_NCYT_) && !allows_flux(_NEXT_))
	{
		UG_THROW("Supplying neither cytosolic nor extracellular natrium concentrations is not allowed.\n"
			 	 "This would mean that the flux calculation would be of no consequence\n"
				 "and this pump mechanism would not do anything.");
	}
}


void MNCX::print_units() const
{
	std::string nm = name();
	size_t n = nm.size();
	UG_LOG(std::endl);
	UG_LOG("+------------------------------------------------------------------------------+"<< std::endl);
	UG_LOG("|  Units used in the implementation of " << nm << std::string(n>=40?0:40-n, ' ') << "|" << std::endl);
	UG_LOG("|------------------------------------------------------------------------------|"<< std::endl);
	UG_LOG("|    Input                                                                     |"<< std::endl);
	UG_LOG("|      [Ca_cyt]  mM (= mol/m^3)                                                |"<< std::endl);
	UG_LOG("|      [Ca_out]  mM (= mol/m^3)                                                |"<< std::endl);
	UG_LOG("|      [Na_cyt]  mM (= mol/m^3)                                                |"<< std::endl);
	UG_LOG("|      [Na_out]  mM (= mol/m^3)                                                |"<< std::endl);
	UG_LOG("|                                                                              |"<< std::endl);
	UG_LOG("|    Output                                                                    |"<< std::endl);
	UG_LOG("|      Ca flux   mol/s                                                         |"<< std::endl);
	UG_LOG("|      Na flux   mol/s                                                         |"<< std::endl);
	UG_LOG("+------------------------------------------------------------------------------+"<< std::endl);
	UG_LOG(std::endl);
}


} // namespace neuro_collection
} // namespace ug


