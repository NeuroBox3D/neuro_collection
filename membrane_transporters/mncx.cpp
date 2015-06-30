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


MNCX::MNCX(const std::vector<std::string>& fcts) : IMembraneTransporter(fcts)
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

MNCX::MNCX(const char* fcts) : IMembraneTransporter(fcts)
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
//	number caCyt = u[_CCYT_];	// cytosolic Ca2+ concentration
//	number caExt = u[_CEXT_];	// extracellular Ca2+ concentration

//	IMPLEMENT ME
}


void MNCX::calc_flux_deriv(const std::vector<number>& u, GridObject* e, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const
{
// 	get values of the unknowns in associated node
//	number caCyt = u[_CCYT_];	// cytosolic Ca2+ concentration
//	number caExt = u[_CEXT_];	// extracellular Ca2+ concentration

//	IMPLEMENT ME
}


const size_t MNCX::n_dependencies() const
{
//	Todo
	return 1;
}


size_t MNCX::n_fluxes() const
{
//	Todo
	return 1;
};


const std::pair<size_t,size_t> MNCX::flux_from_to(size_t flux_i) const
{
	size_t from, to;
	if (allows_flux(_CCYT_)) from = local_fct_index(_CCYT_); else from = InnerBoundaryConstants::_IGNORE_;
	if (allows_flux(_CEXT_)) to = local_fct_index(_CEXT_); else to = InnerBoundaryConstants::_IGNORE_;

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
	UG_LOG("|                                                                              |"<< std::endl);
	UG_LOG("|    Output                                                                    |"<< std::endl);
	UG_LOG("|      Ca flux   mol/s                                                         |"<< std::endl);
	UG_LOG("+------------------------------------------------------------------------------+"<< std::endl);
	UG_LOG(std::endl);
}


} // namespace neuro_collection
} // namespace ug


