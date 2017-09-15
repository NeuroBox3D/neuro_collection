/*
 *  membrane_transport_fv1.cpp
 *
 *  Created on: 20.12.2011
 *      Author: mbreit
 */

#include "membrane_transport_fv1.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"

namespace ug {
namespace neuro_collection {


template<typename TDomain>
MembraneTransportFV1<TDomain>::MembraneTransportFV1(const char* subsets, SmartPtr<IMembraneTransporter> mt)
: FV1InnerBoundaryElemDisc<TDomain>(),
  R(8.314), T(310.0), F(96485.0), m_spMembraneTransporter(mt), m_bNonRegularGrid(false)
{
	// check validity of transporter setup and then lock
	mt->check_and_lock();

	// elem discs (subsets and) functions need only be set after the previous check
	static_cast<IElemDisc<TDomain>*>(this)->set_subsets(subsets);
	static_cast<IElemDisc<TDomain>*>(this)->set_functions(mt->symb_fcts());
}

template<typename TDomain>
MembraneTransportFV1<TDomain>::MembraneTransportFV1(const std::vector<std::string>& subsets, SmartPtr<IMembraneTransporter> mt)
: FV1InnerBoundaryElemDisc<TDomain>(),
  R(8.314), T(310.0), F(96485.0), m_spMembraneTransporter(mt), m_bNonRegularGrid(false)
{
	// check validity of transporter setup and then lock
	mt->check_and_lock();

	// elem discs (subsets and) functions need only be set after the previous check
	static_cast<IElemDisc<TDomain>*>(this)->set_subsets(subsets);
	static_cast<IElemDisc<TDomain>*>(this)->set_functions(mt->symb_fcts());
}

template<typename TDomain>
MembraneTransportFV1<TDomain>::~MembraneTransportFV1()
{
	// nothing to do
}


template<typename TDomain>
void MembraneTransportFV1<TDomain>::set_density_function(SmartPtr<CplUserData<number,dim> > densityFct)
{
	this->m_spDensityFct = densityFct;
}

template<typename TDomain>
void MembraneTransportFV1<TDomain>::set_density_function(const number dens)
{
	set_density_function(make_sp(new ConstUserNumber<dim>(dens)));
}

template<typename TDomain>
void MembraneTransportFV1<TDomain>::set_density_function(const char* name)
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

template<typename TDomain>
void MembraneTransportFV1<TDomain>::set_membrane_transporter(SmartPtr<IMembraneTransporter> mt)
{
	m_spMembraneTransporter = mt;
}


template<typename TDomain>
bool MembraneTransportFV1<TDomain>::fluxDensityFct
(
	const std::vector<LocalVector::value_type>& u,
	GridObject* e,
	const MathVector<dim>& coords,
	int si,
	FluxCond& fc
)
{
	size_t n_flux = m_spMembraneTransporter->n_fluxes();

	// calculate single-channel flux
	fc.flux.resize(n_flux);
	fc.from.resize(n_flux);
	fc.to.resize(n_flux);

	m_spMembraneTransporter->flux(u, e, fc.flux);

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


template<typename TDomain>
bool MembraneTransportFV1<TDomain>::fluxDensityDerivFct
(
	const std::vector<LocalVector::value_type>& u,
	GridObject* e,
	const MathVector<dim>& coords,
	int si,
	FluxDerivCond& fdc
)
{
	size_t n_dep = m_spMembraneTransporter->n_dependencies();
	size_t n_flux = m_spMembraneTransporter->n_fluxes();

	// calculate single-channel flux
	fdc.fluxDeriv.resize(n_flux);
	fdc.from.resize(n_flux);
	fdc.to.resize(n_flux);
	for (size_t i = 0; i < n_flux; i++)
		fdc.fluxDeriv[i].resize(n_dep);

	m_spMembraneTransporter->flux_deriv(u, e, fdc.fluxDeriv);

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

template<typename TDomain>
void MembraneTransportFV1<TDomain>::prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	// remember
	m_bNonRegularGrid = bNonRegularGrid;

	// set assembling functions from base class first
	this->FV1InnerBoundaryElemDisc<TDomain>::prepare_setting(vLfeID, bNonRegularGrid);

	// update assemble functions
	register_all_fv1_funcs();
}


template<typename TDomain>
void MembraneTransportFV1<TDomain>::prep_timestep
(
    number future_time,
    number time,
    VectorProxyBase* upb
)
{
	m_spMembraneTransporter->prep_timestep(future_time, time, upb);
}


template<typename TDomain>
void MembraneTransportFV1<TDomain>::register_all_fv1_funcs()
{
	// register prep_timestep function for all known algebra types
	Register<bridge::CompileAlgebraList>(this);
}





template<typename TDomain>
MembraneTransport1d<TDomain>::MembraneTransport1d(const char* subsets, SmartPtr<IMembraneTransporter> mt)
: IElemDisc<TDomain>("", ""), m_radius(0.0), m_spDensityFct(SPNULL), m_spMembraneTransporter(mt), m_currSI(-1)
{
	// check validity of transporter setup and then lock
	mt->check_and_lock();

	static_cast<IElemDisc<TDomain>*>(this)->set_subsets(subsets);
	static_cast<IElemDisc<TDomain>*>(this)->set_functions(mt->symb_fcts());
}

template<typename TDomain>
MembraneTransport1d<TDomain>::MembraneTransport1d(const std::vector<std::string>& subsets, SmartPtr<IMembraneTransporter> mt)
: IElemDisc<TDomain>("", ""), m_radius(0.0), m_spDensityFct(SPNULL), m_spMembraneTransporter(mt), m_currSI(-1)
{
	// check validity of transporter setup and then lock
	mt->check_and_lock();

	static_cast<IElemDisc<TDomain>*>(this)->set_subsets(subsets);
	static_cast<IElemDisc<TDomain>*>(this)->set_functions(mt->symb_fcts());
}


template<typename TDomain>
void MembraneTransport1d<TDomain>::set_density_function(SmartPtr<CplUserData<number,dim> > densityFct)
{
	this->m_spDensityFct = densityFct;
}

template<typename TDomain>
void MembraneTransport1d<TDomain>::set_density_function(const number dens)
{
	set_density_function(make_sp(new ConstUserNumber<dim>(dens)));
}

template<typename TDomain>
void MembraneTransport1d<TDomain>::set_density_function(const char* name)
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

template<typename TDomain>
void MembraneTransport1d<TDomain>::set_radius(number r)
{
	m_radius = r;
}


template<typename TDomain>
bool MembraneTransport1d<TDomain>::fluxDensityFct
(
	const std::vector<LocalVector::value_type>& u,
	GridObject* e,
	const MathVector<dim>& coords,
	int si,
	FluxCond& fc
)
{
	size_t n_flux = m_spMembraneTransporter->n_fluxes();

	// calculate single-channel flux
	fc.flux.resize(n_flux);
	fc.from.resize(n_flux);
	fc.to.resize(n_flux);

	m_spMembraneTransporter->flux(u, e, fc.flux);

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
		fc.flux[i] *= density * 2.0*PI*m_radius;
		fc.from[i] = m_spMembraneTransporter->flux_from_to(i).first;
		fc.to[i] = m_spMembraneTransporter->flux_from_to(i).second;
	}

	return true;
}


template<typename TDomain>
bool MembraneTransport1d<TDomain>::fluxDensityDerivFct
(
	const std::vector<LocalVector::value_type>& u,
	GridObject* e,
	const MathVector<dim>& coords,
	int si,
	FluxDerivCond& fdc
)
{
	size_t n_dep = m_spMembraneTransporter->n_dependencies();
	size_t n_flux = m_spMembraneTransporter->n_fluxes();

	// calculate single-channel flux
	fdc.fluxDeriv.resize(n_flux);
	fdc.from.resize(n_flux);
	fdc.to.resize(n_flux);
	for (size_t i = 0; i < n_flux; i++)
		fdc.fluxDeriv[i].resize(n_dep);

	m_spMembraneTransporter->flux_deriv(u, e, fdc.fluxDeriv);

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
			fdc.fluxDeriv[i][j].second *= density * 2.0*PI*m_radius;
		fdc.from[i] = m_spMembraneTransporter->flux_from_to(i).first;
		fdc.to[i] = m_spMembraneTransporter->flux_from_to(i).second;
	}

	return true;
}

template<typename TDomain>
void MembraneTransport1d<TDomain>::
prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	// check that Lagrange 1st order
	for (size_t i = 0; i < vLfeID.size(); ++i)
		if (vLfeID[i].type() != LFEID::LAGRANGE || vLfeID[i].order() != 1)
			UG_THROW("MembraneTransport1d: 1st order Lagrange expected.");

	// update assemble functions
	register_assembling_funcs();
}

template<typename TDomain>
bool MembraneTransport1d<TDomain>::
use_hanging() const
{
	return false;
}

template<typename TDomain>
void MembraneTransport1d<TDomain>::prep_timestep
(
    number future_time,
    number time,
    VectorProxyBase* upb
)
{
	m_spMembraneTransporter->prep_timestep(future_time, time, upb);
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void MembraneTransport1d<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{
	m_currSI = si;
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void MembraneTransport1d<TDomain>::
fsh_elem_loop()
{}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void MembraneTransport1d<TDomain>::
prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
{
	// update geometry for this element
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();
	try {geo.update(elem, vCornerCoords, &(this->subset_handler()));}
	UG_CATCH_THROW("FV1InnerBoundaryElemDisc::prep_elem: "
						"Cannot update Finite Volume Geometry.");
}

// assemble stiffness part of Jacobian
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void MembraneTransport1d<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// get finite volume geometry
	const static TFVGeom& fvgeom = GeomProvider<TFVGeom>::get();

	for (size_t i = 0; i < fvgeom.num_scv(); ++i)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = fvgeom.scv(i);

		// get associated node
		const int co = scv.node_id();

		// get solution at the corner of the scv
		size_t nFct = u.num_fct();
		std::vector<LocalVector::value_type> uAtCorner(nFct);
		for (size_t fct = 0; fct < nFct; fct++)
			uAtCorner[fct] = u(fct,co);

		// get corner coordinates
		const MathVector<dim>& cc = scv.global_corner(0);

		FluxDerivCond fdc;
		if (!fluxDensityDerivFct(uAtCorner, elem, cc, m_currSI, fdc))
			UG_THROW("FV1InnerBoundaryElemDisc::add_jac_A_elem:"
							" Call to fluxDensityDerivFct resulted did not succeed.");

		// scale with volume of SCV
		number scale = scv.volume();
		for (size_t j=0; j<fdc.fluxDeriv.size(); j++)
			for (size_t k=0; k<fdc.fluxDeriv[j].size(); k++)
				fdc.fluxDeriv[j][k].second *= scale;

		// add to Jacobian
		for (size_t j=0; j<fdc.fluxDeriv.size(); j++)
		{
			for (size_t k=0; k<fdc.fluxDeriv[j].size(); k++)
			{
				if (fdc.from[j] != InnerBoundaryConstants::_IGNORE_)
					J(fdc.from[j], co, fdc.fluxDeriv[j][k].first, co) += fdc.fluxDeriv[j][k].second;
				if (fdc.to[j] != InnerBoundaryConstants::_IGNORE_)
					J(fdc.to[j], co, fdc.fluxDeriv[j][k].first, co)	-= fdc.fluxDeriv[j][k].second;
			}
		}
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void MembraneTransport1d<TDomain>::
add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void MembraneTransport1d<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// get finite volume geometry
	static TFVGeom& fvgeom = GeomProvider<TFVGeom>::get();

	// loop Boundary Faces
	for (size_t i = 0; i < fvgeom.num_scv(); ++i)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = fvgeom.scv(i);

		// get associated node
		const int co = scv.node_id();

		// get solution at the corner of the scv
		size_t nFct = u.num_fct();
		std::vector<LocalVector::value_type> uAtCorner(nFct);
		for (size_t fct = 0; fct < nFct; fct++)
			uAtCorner[fct] = u(fct,co);

		// get corner coordinates
		const MathVector<dim>& cc = scv.global_corner(0);

		// get flux densities in that node
		FluxCond fc;
		if (!fluxDensityFct(uAtCorner, elem, cc, m_currSI, fc))
		{
			UG_THROW("FV1InnerBoundaryElemDisc::add_def_A_elem:"
						" Call to fluxDensityFct did not succeed.");
		}

		// scale with volume of SCV
		number scale = scv.volume();

		for (size_t j=0; j<fc.flux.size(); j++)
			fc.flux[j] *= scale;

		// add to defect
		for (size_t j=0; j<fc.flux.size(); j++)
		{
			if (fc.from[j] != InnerBoundaryConstants::_IGNORE_) d(fc.from[j], co) += fc.flux[j];
			if (fc.to[j] != InnerBoundaryConstants::_IGNORE_) d(fc.to[j], co) -= fc.flux[j];
		}
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void MembraneTransport1d<TDomain>::
add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void MembraneTransport1d<TDomain>::
add_rhs_elem(LocalVector& rhs, GridObject* elem, const MathVector<dim> vCornerCoords[])
{}




template<typename TDomain>
void MembraneTransport1d<TDomain>::register_assembling_funcs()
{
	// register prep_timestep function for all known algebra types
	RegisterPrepTimestepFct<bridge::CompileAlgebraList>(this);

	// register assembling functionality
	typedef RegularEdge TElem;
	typedef FV1Geometry<RegularEdge, dim> TFVGeom;
	typedef MembraneTransport1d<TDomain> T;
	const ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

	this->set_prep_elem_loop_fct(id, &T::template prep_elem_loop<TElem, TFVGeom>);
	this->set_prep_elem_fct(id, &T::template prep_elem<TElem, TFVGeom>);
	this->set_add_def_A_elem_fct(id, &T::template add_def_A_elem<TElem, TFVGeom>);
	this->set_add_def_M_elem_fct(id, &T::template add_def_M_elem<TElem, TFVGeom>);
	this->set_add_rhs_elem_fct(id, &T::template add_rhs_elem<TElem, TFVGeom>);
	this->set_add_jac_A_elem_fct(id, &T::template add_jac_A_elem<TElem, TFVGeom>);
	this->set_add_jac_M_elem_fct(id, &T::template add_jac_M_elem<TElem, TFVGeom>);
	this->set_fsh_elem_loop_fct(id, &T::template fsh_elem_loop<TElem, TFVGeom>);
}


// explicit template specializations
#ifdef UG_DIM_1
	template class MembraneTransportFV1<Domain1d>;
	template class MembraneTransport1d<Domain1d>;
#endif
#ifdef UG_DIM_2
	template class MembraneTransportFV1<Domain2d>;
	template class MembraneTransport1d<Domain2d>;
#endif
#ifdef UG_DIM_3
	template class MembraneTransportFV1<Domain3d>;
	template class MembraneTransport1d<Domain3d>;
#endif


} // end namespace neuro_collection
} // end namespace ug

