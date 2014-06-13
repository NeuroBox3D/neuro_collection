/*
 *  one_sided_membrane_transport_fv1.h
 *
 *  Created on: 20.12.2011
 *      Author: markusbreit
 */

#ifndef __H__UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__ONE_SIDED_MEMBRANE_TRANSPORT_FV1__
#define __H__UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__ONE_SIDED_MEMBRANE_TRANSPORT_FV1__


#include "dependent_neumann_boundary_fv1.h"
#include "bindings/lua/lua_user_data.h"
#include "common/util/smart_pointer.h"



namespace ug
{
namespace neuro_collection
{

/// finite volume element discretizations for boundary conditions on an outer membrane
/**
 * This class implements the dependent_neumann_boundary interface to provide element local
 * assemblings for the unknown-dependent Neumann flux over an outer membrane; i.e. a
 * membrane where the flowing unknowns are only present on one side.
 */
template<typename TDomain>
class OneSidedMembraneTransportFV1
: public DependentNeumannBoundaryFV1<TDomain>
{
	protected:
		const number R;		// universal gas constant
		const number T;		// temperature
		const number F;		// Faraday constant


	///	world dimension
		static const int dim = DependentNeumannBoundaryFV1<TDomain>::dim;

	public:
	/// constructor
		OneSidedMembraneTransportFV1(const char* functions, const char* subsets)
					: DependentNeumannBoundaryFV1<TDomain>(functions, subsets),
					  R(8.314), T(310.0), F(96485.0) {};

	public:
	/// adding density information for pumps/channels in membrane
		void set_density_function(SmartPtr<UserData<number,dim> > densityFct)
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

	protected:
		SmartPtr<UserData<number,dim> > m_spDensityFct;
};



/// Discretization for the PMCA pump in the PM
/**
 * This class implements the InnerBoundaryElemDisc to provide flux densities
 * and their derivatives for the Graupner (2003) model of PMCA pumps (saturated CaM).
 *
 * \tparam	TDomain		domain
 */
template<typename TDomain>
class OneSidedPMCAFV1
: public OneSidedMembraneTransportFV1<TDomain>
{
	private:
		const number KD_P;					// mol*dm^-3 (Elwess et al.)
		//const number KD_P = 3.4e-07;		// mol*dm^-3 (Graupner)
		const number IMAX_P;				// mol*s^-1
		const number RHO_P;		//3.0e03;	// um^-2

		typedef typename DependentNeumannBoundaryFV1<TDomain>::NFluxCond NFluxCond;
		typedef typename DependentNeumannBoundaryFV1<TDomain>::NFluxDerivCond NFluxDerivCond;

	public:
		///	World dimension
		static const int dim = OneSidedMembraneTransportFV1<TDomain>::dim;

	public:

		OneSidedPMCAFV1(const char* functions, const char* subsets)
			: OneSidedMembraneTransportFV1<TDomain>(functions, subsets),
			  KD_P(6.0e-8), IMAX_P(1.7e-23), RHO_P(500.0) {};

	private:
		// virtual functions inherited from FV1MyNeumannBoundaryElemDisc
		bool fluxDensityFct(const std::vector<LocalVector::value_type>& u, const MathVector<dim>& coords, int si, NFluxCond& fc)
		{
			if (u.size() < 1)
			{
				UG_LOG( "ERROR in 'FV1PMCAElemDisc:fluxFct':"
						" LocalVector u does not have enough functions"
						" (has " << u.size() << ", but needs at least 1).");
				return false;
			}

			number caCyt = u[0];	// cytosolic Ca2+ concentration

			number gatingFactor = caCyt*caCyt / (KD_P*KD_P + caCyt*caCyt);

			number density;
			if (this->m_spDensityFct.valid())
				(*this->m_spDensityFct)(density, coords, this->time(), si);
			else
				density = RHO_P;

			number flux = - density * gatingFactor * IMAX_P;

			// dimensional correction: concentrations are mol/dm^3, but length unit is um
			flux *= 1e15;

			fc.flux.resize(1, 0.0);	fc.flux[0] = flux;
			fc.to.resize(1);		fc.to[0] = 0;		// cytosol

			return true;
		}

		bool fluxDensityDerivFct(const std::vector<LocalVector::value_type>& u, const MathVector<dim>& coords, int si, NFluxDerivCond& fdc)
		{
			// get values of the unknowns in associated node
			number caCyt = u[0];

			number dGating_dCyt = 2*KD_P*KD_P*caCyt / std::pow(KD_P*KD_P + caCyt*caCyt, 2);

			number density;
			if (this->m_spDensityFct.valid())
				(*this->m_spDensityFct)(density, coords, this->time(), si);
			else
				density = RHO_P;

			number d_dCyt = - density * dGating_dCyt * IMAX_P;

			// dimensional correction: concentrations are mol/dm^3, but length unit is um
			d_dCyt *= 1e15;

			// add to Jacobian
			fdc.fluxDeriv.resize(1);
			fdc.fluxDeriv[0].resize(1, 0.0);

			fdc.fluxDeriv[0][0] = d_dCyt;

			fdc.to.resize(1);	fdc.to[0] = 0;		// cytosol

			return true;
		}
};


/// Discretization for the NCX pump in the PM
/**
 * This class implements the InnerBoundaryElemDisc to provide flux densities
 * and their derivatives for the Graupner (2003) model of NCX pumps.
 *
 * \tparam	TDomain		domain
 */
template<typename TDomain>
class OneSidedNCXFV1
: public OneSidedMembraneTransportFV1<TDomain>
{
	private:
		const number KD_N;			// mol*dm^-3
		const number IMAX_N;		// mol*s^-1
		const number RHO_N;//1.0e02	// um^-2

		typedef typename DependentNeumannBoundaryFV1<TDomain>::NFluxCond NFluxCond;
		typedef typename DependentNeumannBoundaryFV1<TDomain>::NFluxDerivCond NFluxDerivCond;

	public:
		///	World dimension
		static const int dim = OneSidedMembraneTransportFV1<TDomain>::dim;

	public:

		OneSidedNCXFV1(const char* functions, const char* subsets)
			: OneSidedMembraneTransportFV1<TDomain>(functions, subsets),
			  KD_N(1.8e-6), IMAX_N(2.5e-21), RHO_N(15.0) {};

	private:
		// virtual functions inherited from FV1MyNeumannBoundaryElemDisc
		bool fluxDensityFct(const std::vector<LocalVector::value_type>& u, const MathVector<dim>& coords, int si, NFluxCond& fc)
		{
			if (u.size() < 1)
			{
				UG_LOG( "ERROR in 'FV1NCXElemDisc:fluxFct':"
						" LocalVector u does not have enough functions"
						" (has " << u.size() << ", but needs at least 1).");
				return false;
			}

			number caCyt = u[0];	// cytosolic Ca2+ concentration

			number gatingFactor = caCyt / (KD_N + caCyt);

			number density;
			if (this->m_spDensityFct.valid())
				(*this->m_spDensityFct)(density, coords, this->time(), si);
			else
				density = RHO_N;

			number flux = - density * gatingFactor * IMAX_N;

			// dimensional correction: concentrations are mol/dm^3, but length unit is um
			flux *= 1e15;

			fc.flux.resize(1, 0.0);	fc.flux[0] = flux;
			fc.to.resize(1);		fc.to[0] = 0;		// cytosol

			return true;
		}

		bool fluxDensityDerivFct(const std::vector<LocalVector::value_type>& u, const MathVector<dim>& coords, int si, NFluxDerivCond& fdc)
		{
			// get values of the unknowns in associated node
			number caCyt = u[0];

			number dGating_dCyt = KD_N / std::pow(KD_N + caCyt, 2);

			number density;
			if (this->m_spDensityFct.valid())
				(*this->m_spDensityFct)(density, coords, this->time(), si);
			else
				density = RHO_N;

			number d_dCyt = - density * dGating_dCyt * IMAX_N;

			// dimensional correction: concentrations are mol/dm^3, but length unit is um
			d_dCyt *= 1e15;

			// add to Jacobian
			fdc.fluxDeriv.resize(1);
			fdc.fluxDeriv[0].resize(1, 0.0);

			fdc.fluxDeriv[0][0] = d_dCyt;

			fdc.to.resize(1);	fdc.to[0] = 0;		// cytosol

			return true;
		}
};


/// Discretization for the leak flux through the plasma membrane
/**
 * This class implements the leak flux through the plasma membrane -
 * beware: constitutes but a wild guess atm.
 *
 * \tparam	TDomain		domain
 */
template<typename TDomain>
class OneSidedPMCalciumLeakFV1
: public OneSidedMembraneTransportFV1<TDomain>
{
	private:
		const number LEAK_PM;	//2.1e-20;	// leakage flux

		typedef typename DependentNeumannBoundaryFV1<TDomain>::NFluxCond NFluxCond;
		typedef typename DependentNeumannBoundaryFV1<TDomain>::NFluxDerivCond NFluxDerivCond;

	public:
		///	World dimension
		static const int dim = OneSidedMembraneTransportFV1<TDomain>::dim;

	public:

		OneSidedPMCalciumLeakFV1(const char* functions, const char* subsets)
			: OneSidedMembraneTransportFV1<TDomain>(functions, subsets), LEAK_PM(6.85e-22) {};

	private:
		// virtual functions inherited from FV1MyNeumannBoundaryElemDisc
		bool fluxDensityFct(const std::vector<LocalVector::value_type>& u, const MathVector<dim>& coords, int si, NFluxCond& fc)
		{
			// cheating a little here: utilizing density for leakage flux "density" constant,
			// since the leakage flux is not really caused by a density, but depends of course
			// on the densities of all the pumps/channels in the membrane and is therefore
			// position/subset dependent
			number flux;
			if (this->m_spDensityFct.valid())
				(*this->m_spDensityFct)(flux, coords, this->time(), si);
			else
				flux = LEAK_PM;

			// dimensional correction: concentrations are mol/dm^3, but length unit is um
			flux *= 1e15;

			fc.flux.resize(1, 0.0);	fc.flux[0] = flux;
			fc.to.resize(1);	fc.to[0] = 0;		// cytosol

			return true;
		}

		bool fluxDensityDerivFct(const std::vector<LocalVector::value_type>& u, const MathVector<dim>& coords, int si, NFluxDerivCond& fdc)
		{
			// add to Jacobian
			fdc.fluxDeriv.resize(1);
			fdc.fluxDeriv[0].resize(1);

			fdc.fluxDeriv[0][0] = 0.0;

			fdc.to.resize(1);	fdc.to[0] = 0;		// cytosol

			return true;
		}
};

} // end namespace neuro_collection
} // end namespace ug


#endif /*__H__UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__ONE_SIDED_MEMBRANE_TRANSPORT_FV1__*/
