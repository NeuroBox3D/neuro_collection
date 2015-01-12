/*
 *  FV1CalciumERElemDisc.h
 *
 *  Created on: 20.12.2011
 *      Author: markusbreit
 */

#ifndef __H__UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__TWO_SIDED_MEMBRANE_TRANSPORT_FV1__
#define __H__UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__TWO_SIDED_MEMBRANE_TRANSPORT_FV1__


#include "lib_disc/spatial_disc/elem_disc/inner_boundary/inner_boundary.h"
#include "bindings/lua/lua_user_data.h"
#include "common/util/smart_pointer.h"
#include "membrane_transporter_interface.h"


namespace ug
{
namespace neuro_collection
{

// forward declaration of IMembraneTransporter
class IMembraneTransporter;

///@addtogroup plugin_neuro_collection
///@{

/// Finite Volume element discretization for the inner BndCond on a two-sided membrane
/**
 * This class implements the InnerBoundary interface to provide element local
 * assemblings for the unknown-dependent Neumann flux over a membrane, where the flowing
 * unknowns are present on both sides of the membrane.
 */

template<typename TDomain>
class TwoSidedMembraneTransportFV1
: public FV1InnerBoundaryElemDisc<TDomain>
{
	protected:
		const number R;			// universal gas constant
		const number T;			// temperature
		const number F;			// Faraday constant

		typedef typename FV1InnerBoundaryElemDisc<TDomain>::FluxCond FluxCond;
		typedef typename FV1InnerBoundaryElemDisc<TDomain>::FluxDerivCond FluxDerivCond;

	public:
	///	world dimension
		static const int dim = FV1InnerBoundaryElemDisc<TDomain>::dim;

	public:
	/// constructor (can be deleted after successful implementation of unified membrane transport) TODO
		TwoSidedMembraneTransportFV1(const char* functions, const char* subsets)
					: FV1InnerBoundaryElemDisc<TDomain>(functions, subsets),
					  R(8.314), T(310.0), F(96485.0){};

	/// constructor
		TwoSidedMembraneTransportFV1(const char* subsets, SmartPtr<IMembraneTransporter> mt)
					: FV1InnerBoundaryElemDisc<TDomain>(),
					  R(8.314), T(310.0), F(96485.0), m_spMembraneTransporter(mt)
		{
			// check validity of transporter setup and then lock
			mt->check_and_lock();

			static_cast<IElemDisc<TDomain>*>(this)->set_subsets(subsets);
			static_cast<IElemDisc<TDomain>*>(this)->set_functions(mt->symb_fcts());
		};

	/// destructor
		virtual ~TwoSidedMembraneTransportFV1() {};

	public:
	/// adding density information for pumps/channels in membrane
		void set_density_function(SmartPtr<CplUserData<number,dim> > densityFct)
		{
			this->m_spDensityFct = densityFct;
		}

	/// adding density information for pumps/channels in membrane
		void set_density_function(const char* name)
		{
			// name must be a valid lua function name conforming to LuaUserNumber specs
			if (LuaUserData<number, dim>::check_callback_returns(name))
			{
				set_density_function(LuaUserDataFactory<number, dim>::create(name));
				return;
			}

			// no match found
			if (!CheckLuaCallbackName(name))
				UG_THROW("Lua-Callback with name '" << name << "' does not exist.");

			// name exists, but wrong signature
			UG_THROW("Cannot find matching callback signature. Use:\n"
					"Number - Callback\n" << (LuaUserData<number, dim>::signature()) << "\n");
		}

		void set_density_function(const number dens)
		{
			set_density_function(make_sp(new ConstUserNumber<dim>(dens)));
		}

	/// set transport mechanism
		void set_membrane_transporter(SmartPtr<IMembraneTransporter> mt)
		{
			m_spMembraneTransporter = mt;
		}

	/// flux assembling routines inherited from FV1InnerBoundaryElemDisc
		virtual bool fluxDensityFct(const std::vector<LocalVector::value_type>& u, const MathVector<dim>& coords, int si, FluxCond& fc)
		{
			size_t n_flux = m_spMembraneTransporter->n_fluxes();

			// calculate single-channel flux
			fc.flux.resize(n_flux);
			fc.from.resize(n_flux);
			fc.to.resize(n_flux);

			m_spMembraneTransporter->flux(u, fc.flux);

			// get density in membrane
			if (!this->m_spDensityFct.valid())
			{
				UG_THROW("No density information available for " << m_spMembraneTransporter->name()
						 << " membrane transport mechanism. Please set using set_density_function().");
			}
			number density;
			(*this->m_spDensityFct)(density, coords, this->time(), si);

			for (size_t i = 0; i < n_flux; i++)
			{
				fc.flux[i] *= density;
				fc.from[i] = m_spMembraneTransporter->flux_from_to(i).first;
				fc.to[i] = m_spMembraneTransporter->flux_from_to(i).second;
			}

			return true;
		}

		bool fluxDensityDerivFct(const std::vector<LocalVector::value_type>& u, const MathVector<dim>& coords, int si, FluxDerivCond& fdc)
		{
			size_t n_dep = m_spMembraneTransporter->n_dependencies();
			size_t n_flux = m_spMembraneTransporter->n_fluxes();

			// calculate single-channel flux
			fdc.fluxDeriv.resize(n_flux);
			fdc.from.resize(n_flux);
			fdc.to.resize(n_flux);
			for (size_t i = 0; i < n_flux; i++)
				fdc.fluxDeriv[i].resize(n_dep);

			m_spMembraneTransporter->flux_deriv(u, fdc.fluxDeriv);

			number density;
			if (this->m_spDensityFct.valid())
				(*this->m_spDensityFct)(density, coords, this->time(), si);
			else
			{
				UG_THROW("No density information available for " << m_spMembraneTransporter->name()
						<< " membrane transport mechanism. Please set using set_density_function().");
			}

			for (size_t i = 0; i < n_flux; i++)
			{
				for (size_t j = 0; j < n_dep; j++)
					fdc.fluxDeriv[i][j].second *= density;
				fdc.from[i] = m_spMembraneTransporter->flux_from_to(i).first;
				fdc.to[i] = m_spMembraneTransporter->flux_from_to(i).second;
			}

			return true;
		}

	protected:
		SmartPtr<CplUserData<number,dim> > m_spDensityFct;
		SmartPtr<IMembraneTransporter> m_spMembraneTransporter;
};



/// Finite Volume Element Discretization for the inner BndCond on an ER membrane
/**
 * This class implements the InnerBoundary interface to provide element local
 * assemblings for the unknown-dependent Calcium-Neumann-flux over the ER.
 */

template<typename TDomain>
class TwoSidedERCalciumTransportFV1
: public TwoSidedMembraneTransportFV1<TDomain>
{
	protected:
		using TwoSidedMembraneTransportFV1<TDomain>::R;
		using TwoSidedMembraneTransportFV1<TDomain>::T;
		using TwoSidedMembraneTransportFV1<TDomain>::F;

		const number REF_CA_ER;	// reference endoplasmatic Ca2+ concentration (for conductances)
		const number CF;		// correction factor: the ER being more or less a smooth cable
								// in our simulations, the real surface is not adequately represented;
								// therefore, the channel/pump densities in the ER membrane are
								// multiplied by this factor

	public:
	///	world dimension
		static const int dim = TwoSidedMembraneTransportFV1<TDomain>::dim;

	public:
	/// constructor
		TwoSidedERCalciumTransportFV1(const char* functions, const char* subsets)
					: TwoSidedMembraneTransportFV1<TDomain>(functions, subsets),
					  REF_CA_ER(2.5e-4), CF(1.0) {};
};




/// Discretization for the IP3R calcium channel in the ER membrane
/**
 * This class implements the InnerBoundaryElemDisc to provide flux densities
 * and their derivatives for the De Young & Keizer (1992) model of IP3R channels.
 *
 * \tparam	TDomain		domain
 */
template<typename TDomain>
class TwoSidedIP3RFV1
: public TwoSidedERCalciumTransportFV1<TDomain>
{
	protected:
		using TwoSidedERCalciumTransportFV1<TDomain>::R;
		using TwoSidedERCalciumTransportFV1<TDomain>::T;
		using TwoSidedERCalciumTransportFV1<TDomain>::F;
		using TwoSidedERCalciumTransportFV1<TDomain>::REF_CA_ER;
		using TwoSidedERCalciumTransportFV1<TDomain>::CF;

		const number D1;		// IP3 binding (w/o Ca2+ inhibition)
		const number D2;		// Ca2+ inhibiting binding
		const number D3;		// IP3 binding (w/ Ca2+ inhibition)
		const number D5;		// Ca2+ activating binding
		const number RHO_IP3R;	// average channel density for IP3R channel
		const number MU_IP3R;	// IP3R channel conductance

	private:
		typedef typename FV1InnerBoundaryElemDisc<TDomain>::FluxCond FluxCond;
		typedef typename FV1InnerBoundaryElemDisc<TDomain>::FluxDerivCond FluxDerivCond;

	public:
		///	World dimension
		static const int dim = TwoSidedERCalciumTransportFV1<TDomain>::dim;

	public:

		TwoSidedIP3RFV1(const char* functions, const char* subsets)
			: TwoSidedERCalciumTransportFV1<TDomain>(functions, subsets),
			  D1(1.3e-7), D2(1.05e-6), D3(9.4e-7), D5(8.23e-8), RHO_IP3R(17.3), MU_IP3R(1.6e-12) {};

	private:
		// virtual functions inherited from FV1InnerBoundaryElemDisc
		bool fluxDensityFct(const std::vector<LocalVector::value_type>& u, const MathVector<dim>& coords, int si, FluxCond& fc)
		{
			if (u.size() < 3)
			{
				UG_LOG( "ERROR in 'FV1IP3RElemDisc:fluxFct':"
						" LocalVector u does not have enough functions"
						" (has " << u.size() << ", but needs at least 3).");
				return false;
			}

			number caCyt = u[0];	// cytosolic Ca2+ concentration
			number caER = u[1];		// ER Ca2+ concentration
			number ip3 = u[2];		// IP3 concentration

			// membrane current corresponding to diffusion pressure
			number current = R*T/(4*F*F) * MU_IP3R/REF_CA_ER * (caER - caCyt);

			// flux calculation
			number pOpen = pow(caCyt*ip3*D2 / ((caCyt*ip3 + ip3*D2 + D1*D2 + caCyt*D3) * (caCyt+D5)),3);

			number density;
			if (this->m_spDensityFct.valid())
				(*this->m_spDensityFct)(density, coords, this->time(), si);
			else
				density = CF * RHO_IP3R;

			number flux = density * pOpen * current;

			// dimensional correction: concentrations are mol/dm^3, but length unit is um
			flux *= 1e15;

			fc.flux.resize(1);	fc.flux[0] = flux;
			fc.from.resize(1);	fc.from[0] = 1;		// ER
			fc.to.resize(1);	fc.to[0] = 0;		// cytosol

			return true;
		}

		bool fluxDensityDerivFct(const std::vector<LocalVector::value_type>& u, const MathVector<dim>& coords, int si, FluxDerivCond& fdc)
		{
			// get values of the unknowns in associated node
			number caCyt = u[0];
			number caER = u[1];
			number ip3 = u[2];

			// membrane potential
			number current = R*T/(4*F*F) * MU_IP3R/REF_CA_ER * (caER - caCyt);
			number dI_dCyt = - R*T/(4*F*F) * MU_IP3R/REF_CA_ER;
			number dI_dER = -dI_dCyt;

			// IP3R flux derivatives
			number schlonz1 = caCyt*ip3 + ip3*D2 + D1*D2 + caCyt*D3;
			number schlonz2 = schlonz1 * (caCyt+D5);
			number schlonz3 = (caCyt*ip3*D2) / schlonz2;
			number pOpen = pow(schlonz3,3);

			number density;
			if (this->m_spDensityFct.valid())
				(*this->m_spDensityFct)(density, coords, this->time(), si);
			else
				density = CF * RHO_IP3R;

			number dOpen_dCyt = 3.0*schlonz3*schlonz3*ip3*D2 * (1.0 - caCyt/schlonz2 * ( (ip3+D3)*(caCyt+D5) + schlonz1 )) / schlonz2;
			number dOpen_dIP3 = 3.0*schlonz3*schlonz3*caCyt*D2 * (1.0 - ip3/schlonz2 * (caCyt+D2)*(caCyt+D5)) / schlonz2;

			number d_dCyt = density * (dOpen_dCyt * current + pOpen * dI_dCyt);
			number d_dER  = density * pOpen * dI_dER;
			number d_dIP3 = density * dOpen_dIP3 * current;

			// dimensional correction: concentrations are mol/dm^3, but length unit is um
			d_dCyt *= 1e15;
			d_dER *= 1e15;
			d_dIP3 *= 1e15;

			// add to Jacobian
			fdc.fluxDeriv.resize(1);
			fdc.fluxDeriv[0].resize(3);

			for (int i = 0; i < 3; i++)
				fdc.fluxDeriv[0][i].first = i;

			fdc.fluxDeriv[0][0].second = d_dCyt;
			fdc.fluxDeriv[0][1].second = d_dER;
			fdc.fluxDeriv[0][2].second = d_dIP3;

			fdc.from.resize(1);	fdc.from[0] = 1;	// ER
			fdc.to.resize(1);	fdc.to[0] = 0;		// cytosol

			return true;
		}
};


/// Discretization for the RyR calcium channel in the ER membrane
/**
 * This class implements the InnerBoundaryElemDisc to provide flux densities
 * and their derivatives for the Keizer & Levine (1996) RyR model.
 *
 * \tparam	TDomain		domain
 */
template<typename TDomain>
class TwoSidedRyRFV1
: public TwoSidedERCalciumTransportFV1<TDomain>
{
	protected:
		using TwoSidedERCalciumTransportFV1<TDomain>::R;
		using TwoSidedERCalciumTransportFV1<TDomain>::T;
		using TwoSidedERCalciumTransportFV1<TDomain>::F;
		using TwoSidedERCalciumTransportFV1<TDomain>::REF_CA_ER;
		using TwoSidedERCalciumTransportFV1<TDomain>::CF;

		const number KA;		// calcium binding (C1 <--> O1)
		const number KB;		// calcium binding (O1 <--> O2)
		const number KC;		// O1 <--> C2
		const number RHO_RYR;	// average channel density for RyR channel
		const number MU_RYR;	// RyR channel conductance

	private:
		typedef typename FV1InnerBoundaryElemDisc<TDomain>::FluxCond FluxCond;
		typedef typename FV1InnerBoundaryElemDisc<TDomain>::FluxDerivCond FluxDerivCond;

	public:
		///	World dimension
		static const int dim = TwoSidedERCalciumTransportFV1<TDomain>::dim;

	public:

		TwoSidedRyRFV1(const char* functions, const char* subsets)
			: TwoSidedERCalciumTransportFV1<TDomain>(functions, subsets),
			  KA(5.21e25), KB(3.89e18), KC(17.5), RHO_RYR(0.86), MU_RYR(5.0e-11) {};

	private:
		// virtual functions inherited from FV1InnerBoundaryElemDisc
		bool fluxDensityFct(const std::vector<LocalVector::value_type>& u, const MathVector<dim>& coords, int si, FluxCond& fc)
		{
			if (u.size() < 2)
			{
				UG_LOG( "ERROR in 'FV1RyRElemDisc:fluxFct':"
						" LocalVector u does not have enough functions"
						" (has " << u.size() << ", but needs at least 2).");
				return false;
			}

			number caCyt = u[0];	// cytosolic Ca2+ concentration
			number caER = u[1];		// ER Ca2+ concentration

			// RyR flux calculation
			number current = R*T/(4*F*F) * MU_RYR/REF_CA_ER * (caER - caCyt);
			number pOpen =  (1.0 + KB*pow(caCyt,3)) / (1.0 + KC + 1.0/(KA*pow(caCyt,4)) + KB*pow(caCyt,3));

			number density;
			if (this->m_spDensityFct.valid())
				(*this->m_spDensityFct)(density, coords, this->time(), si);
			else
				density = CF * RHO_RYR;

			number flux = density * pOpen * current;

			// dimensional correction: concentrations are mol/dm^3, but length unit is um
			flux *= 1e15;

			fc.flux.resize(1, 0.0);	fc.flux[0] = flux;
			fc.from.resize(1);	fc.from[0] = 1;		// ER
			fc.to.resize(1);	fc.to[0] = 0;		// cytosol

			return true;
		}

		bool fluxDensityDerivFct(const std::vector<LocalVector::value_type>& u, const MathVector<dim>& coords, int si, FluxDerivCond& fdc)
		{
			// get values of the unknowns in associated node
			number caCyt = u[0];
			number caER = u[1];

			// membrane potential
			number current = R*T/(4*F*F) * MU_RYR/REF_CA_ER * (caER - caCyt);
			number dI_dCyt = - R*T/(4*F*F) * MU_RYR/REF_CA_ER;
			number dI_dER = -dI_dCyt;

			// RyR flux derivs
			number schlonz1 = 1.0 + KB*pow(caCyt,3);
			number schlonz2 = 1.0 + KC + 1.0/(KA*pow(caCyt,4)) + KB*pow(caCyt,3);
			number pOpen = schlonz1 / schlonz2;

			number density;
			if (this->m_spDensityFct.valid())
				(*this->m_spDensityFct)(density, coords, this->time(), si);
			else
				density = CF * RHO_RYR;

			number dOpen_dCyt = (3.0*KB*caCyt*caCyt + schlonz1/schlonz2*(4.0/(KA*pow(caCyt,5)) - 3.0*KB*caCyt*caCyt)) / schlonz2;

			number d_dCyt = density * (dOpen_dCyt * current + pOpen * dI_dCyt);
			number d_dER  = density * pOpen * dI_dER;

			// dimensional correction: concentrations are mol/dm^3, but length unit is um
			d_dCyt *= 1e15;
			d_dER *= 1e15;

			// add to Jacobian
			fdc.fluxDeriv.resize(1);
			fdc.fluxDeriv[0].resize(2);

			for (int i = 0; i < 2; i++)
				fdc.fluxDeriv[0][i].first = i;

			fdc.fluxDeriv[0][0].second = d_dCyt;
			fdc.fluxDeriv[0][1].second = d_dER;

			fdc.from.resize(1);	fdc.from[0] = 1;	// ER
			fdc.to.resize(1);	fdc.to[0] = 0;		// cytosol

			return true;
		}
};


/// Discretization for the SERCA calcium pump in the ER membrane
/**
 * This class implements the InnerBoundaryElemDisc to provide flux densities
 * and their derivatives for the Sneyd et al. (2003) SERCA model.
 *
 * \tparam	TDomain		domain
 */
template<typename TDomain>
class TwoSidedSERCAFV1
: public TwoSidedERCalciumTransportFV1<TDomain>
{
	public:
		using TwoSidedERCalciumTransportFV1<TDomain>::CF;

		const number VS;			// maxmimal transport power per pump (divided by [ca_ER])
		const number KS;			// concentration at which halfmaximal pumping occurs
		const number RHO_SERCA;	// average pump density (fit to zero-flux at resting)

	private:
		typedef typename FV1InnerBoundaryElemDisc<TDomain>::FluxCond FluxCond;
		typedef typename FV1InnerBoundaryElemDisc<TDomain>::FluxDerivCond FluxDerivCond;

	public:
		///	World dimension
		static const int dim = TwoSidedERCalciumTransportFV1<TDomain>::dim;

	public:

		TwoSidedSERCAFV1(const char* functions, const char* subsets)
			: TwoSidedERCalciumTransportFV1<TDomain>(functions, subsets),
			  VS(6.5e-27), KS(1.8e-7), RHO_SERCA(1973.0) {};

	private:
		// virtual functions inherited from FV1InnerBoundaryElemDisc
		bool fluxDensityFct(const std::vector<LocalVector::value_type>& u, const MathVector<dim>& coords, int si, FluxCond& fc)
		{
			if (u.size() < 2)
			{
				UG_LOG( "ERROR in 'FV1SERCAElemDisc:fluxFct':"
						" LocalVector u does not have enough functions"
						" (has " << u.size() << ", but needs at least 2).");
				return false;
			}

			number caCyt = u[0];	// cytosolic Ca2+ concentration
			number caER = u[1];		// ER Ca2+ concentration

			// SERCA parts
			number density;
			if (this->m_spDensityFct.valid())
				(*this->m_spDensityFct)(density, coords, this->time(), si);
			else
				density = CF * RHO_SERCA;
			number flux = -density * VS*caCyt / ((KS + caCyt) * caER);

			// dimensional correction: concentrations are mol/dm^3, but length unit is um
			flux *= 1e15;

			fc.flux.resize(1, 0.0);	fc.flux[0] = flux;
			fc.from.resize(1);	fc.from[0] = 1;		// ER
			fc.to.resize(1);	fc.to[0] = 0;		// cytosol

			return true;
		}

		bool fluxDensityDerivFct(const std::vector<LocalVector::value_type>& u, const MathVector<dim>& coords, int si, FluxDerivCond& fdc)
		{
			// get values of the unknowns in associated node
			number caCyt = u[0];
			number caER = u[1];

			// SERCA flux derivs
			number density;
			if (this->m_spDensityFct.valid())
				(*this->m_spDensityFct)(density, coords, this->time(), si);
			else
				density = CF * RHO_SERCA;

			number d_dCyt = - density * (VS*KS) / (pow(KS+caCyt,2)*caER);
			number d_dER  = density * VS*caCyt / ((KS+caCyt)*caER*caER);

			// dimensional correction: concentrations are mol/dm^3, but length unit is um
			d_dCyt *= 1e15;
			d_dER *= 1e15;

			// add to Jacobian
			fdc.fluxDeriv.resize(1);
			fdc.fluxDeriv[0].resize(2);

			for (int i = 0; i < 2; i++)
				fdc.fluxDeriv[0][i].first = i;

			fdc.fluxDeriv[0][0].second = d_dCyt;
			fdc.fluxDeriv[0][1].second = d_dER;

			fdc.from.resize(1);	fdc.from[0] = 1;	// ER
			fdc.to.resize(1);	fdc.to[0] = 0;		// cytosol

			return true;
		}
};


/// Discretization for the leak flux through the ER membrane
/**
 * This class implements the leak flux through the ER membrane -
 * beware: constitutes but a wild guess atm.
 *
 * \tparam	TDomain		domain
 */
template<typename TDomain>
class TwoSidedERCalciumLeakFV1
: public TwoSidedERCalciumTransportFV1<TDomain>
{
	private:
		using TwoSidedERCalciumTransportFV1<TDomain>::CF;

		const number LEAK_ER;	// leak flux constant (leakFluxDensity = LEAK * (ca_er-ca_cyt))

	private:
		typedef typename FV1InnerBoundaryElemDisc<TDomain>::FluxCond FluxCond;
		typedef typename FV1InnerBoundaryElemDisc<TDomain>::FluxDerivCond FluxDerivCond;

	public:
		///	World dimension
		static const int dim = TwoSidedERCalciumTransportFV1<TDomain>::dim;

	public:

		TwoSidedERCalciumLeakFV1(const char* functions, const char* subsets)
			: TwoSidedERCalciumTransportFV1<TDomain>(functions, subsets), LEAK_ER(3.4e-17) {};

	private:
		// virtual functions inherited from FV1InnerBoundaryElemDisc
		bool fluxDensityFct(const std::vector<LocalVector::value_type>& u, const MathVector<dim>& coords, int si, FluxCond& fc)
		{
			if (u.size() < 2)
			{
				UG_LOG( "ERROR in 'FV1ERLeakElemDisc:fluxFct':"
						" LocalVector u does not have enough functions"
						" (has " << u.size() << ", but needs at least 2).");
				return false;
			}

			number caCyt = u[0];	// cytosolic Ca2+ concentration
			number caER = u[1];		// ER Ca2+ concentration

			// membrane current corresponding to diffusion pressure
			// cheating a little here: utilizing density for leakage flux "density" constant,
			// since the leakage flux is not really caused by a density, but depends of course
			// on the densities of all the pumps/channels in the membrane and is therefore
			// position/subset dependent
			number density;
			if (this->m_spDensityFct.valid())
				(*this->m_spDensityFct)(density, coords, this->time(), si);
			else
				density = CF * LEAK_ER;
			number flux = density * (caER-caCyt);

			// dimensional correction: concentrations are mol/dm^3, but length unit is um
			flux *= 1e15;

			fc.flux.resize(1);	fc.flux[0] = flux;
			fc.from.resize(1);	fc.from[0] = 1;		// ER
			fc.to.resize(1);	fc.to[0] = 0;		// cytosol

			return true;
		}

		bool fluxDensityDerivFct(const std::vector<LocalVector::value_type>& u, const MathVector<dim>& coords, int si, FluxDerivCond& fdc)
		{
			number density;
			if (this->m_spDensityFct.valid())
				(*this->m_spDensityFct)(density, coords, this->time(), si);
			else
				density = CF * LEAK_ER;

			number d_dER  = density;
			number d_dCyt = -density;

			// dimensional correction: concentrations are mol/dm^3, but length unit is um
			d_dCyt *= 1e15;
			d_dER *= 1e15;

			// add to Jacobian
			fdc.fluxDeriv.resize(1);
			fdc.fluxDeriv[0].resize(2);

			for (int i = 0; i < 2; i++)
				fdc.fluxDeriv[0][i].first = i;

			fdc.fluxDeriv[0][0].second = d_dCyt;
			fdc.fluxDeriv[0][1].second = d_dER;

			fdc.from.resize(1);	fdc.from[0] = 1;	// ER
			fdc.to.resize(1);	fdc.to[0] = 0;		// cytosol

			return true;
		}
};

///@}

} // end namespace neuro_collection
} // end namespace ug


#endif /*__H__UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__TWO_SIDED_MEMBRANE_TRANSPORT_FV1__*/
