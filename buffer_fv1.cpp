/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2012-11-07
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

#include "buffer_fv1.h"

namespace ug {
namespace neuro_collection {


template<typename TDomain>
BufferFV1<TDomain>::BufferFV1(const char* subsets)
: IElemDisc<TDomain>(NULL, subsets), m_bNonRegularGrid(false)
{
	register_all_fv1_funcs();
	m_reactions.reserve(1);
}

template<typename TDomain>
BufferFV1<TDomain>::~BufferFV1()
{
	// nothing to do
}

template<typename TDomain>
void BufferFV1<TDomain>::set_num_reactions(size_t n)
{
	UG_COND_THROW(m_reactions.size(),
		"Number of buffering reactions must be set before adding the first reaction.");

	// we need to ensure the vector is big enough for all elements,
	// as we cannot allow re-allocation of its elements
	// (pointers to DataImports are held elsewhere through the register_import() method)
	m_reactions.reserve(n);
}


template<typename TDomain>
void BufferFV1<TDomain>::prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	// check number
	if (vLfeID.size() != this->num_fct())
		UG_THROW("FV1BufferElemDisc: needs exactly " << this->num_fct() << " functions.");

	// check that Lagrange 1st order
	for (size_t i = 0; i < vLfeID.size(); ++i)
		if (vLfeID[i].order() != 1 || vLfeID[i].type() != LFEID::LAGRANGE)
			UG_THROW("FV1BufferElemDisc FV scheme only implemented for 1st order Lagrange,\n"
				"but function " << i << " is of type " << vLfeID[i].type()
				<< " and order " << vLfeID[i].order() << ".");

	// remember
	m_bNonRegularGrid = bNonRegularGrid;

	// update assemble functions
	register_all_fv1_funcs();
}

template<typename TDomain>
bool BufferFV1<TDomain>::use_hanging() const
{
	return true;
}


template<typename TDomain>
void
BufferFV1<TDomain>::
add_reaction(const char* fct1, const char* fct2,
			 SmartPtr<CplUserData<number, dim> > tbc,
			 SmartPtr<CplUserData<number, dim> > k1,
			 SmartPtr<CplUserData<number, dim> > k2)
{
	// check that we still have enough space in the reaction vector
	UG_COND_THROW(m_reactions.size() >= m_reactions.capacity(),
		"Reaction cannot be added, admissible number of reaction terms ("
		<< m_reactions.capacity() << ") is reached.\n"
		"Use the set_num_reactions() method before adding the first reaction to accommodate more.");

	// determine function indices
	// check if agents already in functions schema; if not: add
	std::vector<std::string> fcts = this->symb_fcts();
	size_t fctIndex1 = 0; bool found1 = false;
	size_t fctIndex2 = 0; bool found2 = false;
	for (size_t i = 0; i < fcts.size(); i++)
	{
		if (!fcts[i].compare(fct1)) // compare evaluates to 0 iff equal
		{
			fctIndex1 = i;
			found1 = true;
			if (found2) break;
		}

		if (!fcts[i].compare(fct2))
		{
			fctIndex2 = i;
			found2 = true;
			if (found1) break;
		}
	}

	if (!found1)
	{
		fctIndex1 = fcts.size();
		fcts.push_back(std::string(fct1));
	}
	if (!found2)
	{
		fctIndex2 = fcts.size();
		fcts.push_back(std::string(fct2));
	}

	// save functions
	this->set_functions(fcts);

	// set entry in reactions vector
	m_reactions.push_back(ReactionInfo<dim>(fctIndex1, fctIndex2, tbc, k1, k2));

	// register new DataImports
	this->register_import(m_reactions[m_reactions.size()-1].tot_buffer);
	this->register_import(m_reactions[m_reactions.size()-1].k_bind);
	this->register_import(m_reactions[m_reactions.size()-1].k_unbind);
}


#ifdef UG_FOR_LUA
template<typename TDomain>
void
BufferFV1<TDomain>::
add_reaction(const char* fct1, const char* fct2, number tbc, const char* k1, const char* k2)
{
	add_reaction(fct1, fct2,
				 make_sp(new ConstUserNumber<dim>(tbc)),
				 LuaUserDataFactory<number,dim>::create(k1),
				 LuaUserDataFactory<number,dim>::create(k2));
}
#endif


template<typename TDomain>
void
BufferFV1<TDomain>::
add_reaction(const char* fct1, const char* fct2, number tbc, number k1, number k2)
{
	add_reaction(fct1, fct2,
				 make_sp(new ConstUserNumber<dim>(tbc)),
				 make_sp(new ConstUserNumber<dim>(k1)),
				 make_sp(new ConstUserNumber<dim>(k2)));
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void BufferFV1<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{
	//	set local positions
	if (!TFVGeom::usesHangingNodes)
	{
		static const int refDim = TElem::dim;
		static const TFVGeom& geo = GeomProvider<TFVGeom>::get();
		const MathVector<refDim>* vSCVip = geo.scv_local_ips();
		const size_t numSCVip = geo.num_scv_ips();

		for (size_t i=0; i<m_reactions.size(); i++)
		{
			m_reactions[i].tot_buffer.template set_local_ips<refDim>(vSCVip, numSCVip, false);
			m_reactions[i].k_bind.template set_local_ips<refDim>(vSCVip, numSCVip, false);
			m_reactions[i].k_unbind.template set_local_ips<refDim>(vSCVip, numSCVip, false);
		}
	}
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void BufferFV1<TDomain>:: fsh_elem_loop()
{}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void BufferFV1<TDomain>::
prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
{
	// update geometry for this element
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();
	try {geo.update(elem, vCornerCoords, &(this->subset_handler()));}
	UG_CATCH_THROW("FV1BufferElemDisc::prep_elem: Cannot update Finite Volume Geometry.");

	// set local positions
	if (TFVGeom::usesHangingNodes)
	{
		const int refDim = TElem::dim;
		const MathVector<refDim>* vSCVip = geo.scv_local_ips();
		const size_t numSCVip = geo.num_scv_ips();

		for (size_t i=0; i<m_reactions.size(); i++)
		{
			m_reactions[i].tot_buffer.template set_local_ips<refDim>(vSCVip, numSCVip);
			m_reactions[i].k_bind.template set_local_ips<refDim>(vSCVip, numSCVip);
			m_reactions[i].k_unbind.template set_local_ips<refDim>(vSCVip, numSCVip);
		}
	}

	//	set global positions
	const MathVector<dim>* vSCVip = geo.scv_global_ips();
	const size_t numSCVip = geo.num_scv_ips();

	for (size_t i=0; i<m_reactions.size(); i++)
	{
		m_reactions[i].tot_buffer.set_global_ips(vSCVip, numSCVip);
		m_reactions[i].k_bind.set_global_ips(vSCVip, numSCVip);
		m_reactions[i].k_unbind.set_global_ips(vSCVip, numSCVip);
	}
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void BufferFV1<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// get finite volume geometry
	static const TFVGeom& fvgeom = GeomProvider<TFVGeom>::get();

	// loop subcontrol volumes
	for (size_t ip = 0; ip < fvgeom.num_scv(); ++ip)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = fvgeom.scv(ip);

		// get associated node
		const int co = scv.node_id();

		// loop reactions
		for (size_t j = 0; j < m_reactions.size(); j++)
		{
			// compute local defect
			const struct ReactionInfo<dim>& r = m_reactions[j];
			UG_ASSERT(r.k_bind.data_given() && r.k_unbind.data_given() &&  r.tot_buffer.data_given(),
					  "Data import for buffering reaction has no data.");
			number def =   r.k_bind[ip] * u(r.buffer, co) * u(r.buffered, co)
						 - r.k_unbind[ip] * (r.tot_buffer[ip] - u(r.buffer, co));

			// scale with scv volume and add to defect
			d(r.buffer, co) +=	def * scv.volume();
			d(r.buffered, co) += def * scv.volume();
		}
	}
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void BufferFV1<TDomain>::
add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{}


// assemble stiffness part of Jacobian
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void BufferFV1<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// get finite volume geometry
	static const TFVGeom& fvgeom = GeomProvider<TFVGeom>::get();

	// loop scvs
	for (size_t ip = 0; ip < fvgeom.num_scv(); ++ip)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = fvgeom.scv(ip);

		// get associated node
		const int co = scv.node_id();

		// loop reactions
		for (size_t j = 0; j < m_reactions.size(); j++)
		{
			// compute local derivatives
			const struct ReactionInfo<dim>& r = m_reactions[j];
			UG_ASSERT(r.k_bind.data_given() && r.k_unbind.data_given(), "Data import for buffering reaction has no data.");
			number d_dBuff  = r.k_bind[ip] * u(r.buffered, co) + r.k_unbind[ip];
			number d_dBuffd = r.k_bind[ip] * u(r.buffer, co);

			// scale with scv volume and add to Jacobian
			J(r.buffer, co, r.buffer, co)     += d_dBuff * scv.volume();
			J(r.buffer, co, r.buffered, co)   += d_dBuffd * scv.volume();
			J(r.buffered, co, r.buffer, co)   += d_dBuff * scv.volume();
			J(r.buffered, co, r.buffered, co) += d_dBuffd * scv.volume();
		}
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void BufferFV1<TDomain>::
add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void BufferFV1<TDomain>::
add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
{}


// ///////////////////////////////
//   error estimation (begin)   //

//	prepares the loop over all elements of one type for the computation of the error estimator
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void BufferFV1<TDomain>::
prep_err_est_elem_loop(const ReferenceObjectID roid, const int si)
{
//	get the error estimator data object and check that it is of the right type
//	we check this at this point in order to be able to dispense with this check later on
//	(i.e. in prep_err_est_elem and compute_err_est_A_elem())
	if (this->m_spErrEstData.get() == NULL)
	{
		UG_THROW("No ErrEstData object has been given to this ElemDisc!");
	}

	err_est_type* err_est_data = dynamic_cast<err_est_type*>(this->m_spErrEstData.get());

	if (!err_est_data)
	{
		UG_THROW("Dynamic cast to MultipleSideAndElemErrEstData failed."
				<< std::endl << "Make sure you handed the correct type of ErrEstData to this discretization.");
	}

	if (!err_est_data->equal_elem_order())
	{
		UG_THROW("The underlying error estimator data objects of this discretization's "
				 "error estimator do not all have the same integration orders. This case "
				 "is not supported by the implementation. If you need it, implement!");
	}

	if (m_reactions.size() > 0 && err_est_data->num() < 1)
	{
		UG_THROW("No underlying error estimator data objects present. No IPs can be determined.");
	}

//	set local positions
	if (!TFVGeom::usesHangingNodes)
	{
		static const int refDim = TElem::dim;

		// I think, this is not strictly necessary. This will probably never happen.
		// And even if it did, I suppose, it is no harm.
		if (dim != refDim)
		{
			UG_THROW("Dimension of the element this disc is assembled for is not the same as world dimension. "
					"This should not happen.");
		}

		// get local IPs
		size_t numElemIPs;
		const MathVector<refDim>* elemIPs;
		try
		{
			// take IPs from first underlying object (we have forced the others to be the same)
			numElemIPs = err_est_data->get(0)->num_elem_ips(roid);
			elemIPs = err_est_data->get(0)->template elem_local_ips<refDim>(roid);

			if (!elemIPs) return;	// (is NULL if TElem is not of the same dim as TDomain
		}
		UG_CATCH_THROW("Integration points for error estimator cannot be set.");

		// set local IPs in imports
		for (size_t i=0; i<m_reactions.size(); i++)
		{
			m_reactions[i].tot_buffer.template set_local_ips<refDim>(elemIPs, numElemIPs, false);
			m_reactions[i].k_bind.template set_local_ips<refDim>(elemIPs, numElemIPs, false);
			m_reactions[i].k_unbind.template set_local_ips<refDim>(elemIPs, numElemIPs, false);
		}

		// store values of shape functions in local IPs
		LagrangeP1<typename reference_element_traits<TElem>::reference_element_type> trialSpace
					= Provider<LagrangeP1<typename reference_element_traits<TElem>::reference_element_type> >::get();

		m_shapeValues.resize(numElemIPs, trialSpace.num_sh());
		for (size_t ip = 0; ip < numElemIPs; ip++)
			trialSpace.shapes(m_shapeValues.shapesAtElemIP(ip), elemIPs[ip]);
	}
};

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void BufferFV1<TDomain>::
prep_err_est_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
//	error estimator
	err_est_type* err_est_data = dynamic_cast<err_est_type*>(this->m_spErrEstData.get());

//	roid
	ReferenceObjectID roid = elem->reference_object_id();

//	set local positions
	if (TFVGeom::usesHangingNodes)
	{
		static const int refDim = TElem::dim;

		// I think, this is not strictly necessary. This will probably never happen.
		// And even if it did, I suppose, it is no harm.
		if (dim != refDim)
		{
			UG_THROW("Dimension of the element this disc is assembled for is not the same as world dimension. "
					"This should not happen.");
		}

		size_t numElemIPs;
		const MathVector<refDim>* elemIPs;
		try
		{
			// take IPs from first underlying object (we have forced the others to be the same)
			numElemIPs = err_est_data->get(0)->num_elem_ips(roid);
			elemIPs = err_est_data->get(0)->template elem_local_ips<refDim>(roid);

			if (!elemIPs) return;	// is NULL if TElem is not of the same dim as TDomain
		}
		UG_CATCH_THROW("Integration points for error estimator cannot be set.");

		// set local IPs in imports
		for (size_t i=0; i<m_reactions.size(); i++)
		{
			m_reactions[i].tot_buffer.template set_local_ips<refDim>(elemIPs, numElemIPs);
			m_reactions[i].k_bind.template set_local_ips<refDim>(elemIPs, numElemIPs);
			m_reactions[i].k_unbind.template set_local_ips<refDim>(elemIPs, numElemIPs);
		}

		// store values of shape functions in local IPs
		LagrangeP1<typename reference_element_traits<TElem>::reference_element_type> trialSpace
					= Provider<LagrangeP1<typename reference_element_traits<TElem>::reference_element_type> >::get();

		m_shapeValues.resize(numElemIPs, trialSpace.num_sh());
		for (size_t ip = 0; ip < numElemIPs; ip++)
			trialSpace.shapes(m_shapeValues.shapesAtElemIP(ip), elemIPs[ip]);
	}

//	set global positions
	size_t numElemIPs;
	MathVector<dim>* elemIPs;

	try
	{
		numElemIPs = err_est_data->get(0)->num_elem_ips(roid);
		elemIPs = err_est_data->get(0)->elem_global_ips(elem, vCornerCoords);
	}
	UG_CATCH_THROW("Global integration points for error estimator cannot be set.");

	// set local IPs in imports
	for (size_t i=0; i<m_reactions.size(); i++)
	{
		m_reactions[i].tot_buffer.set_global_ips(&elemIPs[0], numElemIPs);
		m_reactions[i].k_bind.set_global_ips(&elemIPs[0], numElemIPs);
		m_reactions[i].k_unbind.set_global_ips(&elemIPs[0], numElemIPs);
	}
}


//	computes the error estimator contribution for one element
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void BufferFV1<TDomain>::
compute_err_est_A_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale)
{
	// we have only elem parts here, no integral over a side

	// err est data object
	err_est_type* err_est_data = dynamic_cast<err_est_type*>(this->m_spErrEstData.get());
	if (err_est_data->num() < this->symb_fcts().size())
	{
		UG_THROW("MultipleSideAndElemErrEstData object does not contain enough error estimators:"
				<< std::endl << "Needs at least " << this->symb_fcts().size() <<
				", but has only " << err_est_data->num() << ".");
	}

	if (err_est_data->get(0)->surface_view().get() == NULL) {UG_THROW("Error estimator has NULL surface view.");}
	MultiGrid* pErrEstGrid = (MultiGrid*) (err_est_data->get(0)->surface_view()->subset_handler()->multi_grid());

	typename MultiGrid::traits<typename SideAndElemErrEstData<TDomain>::elem_type>::secure_container elem_list;
	pErrEstGrid->associated_elements_sorted(elem_list, (TElem*) elem);
	if (elem_list.size() != 1)
	{
		UG_THROW ("There should only be exactly one elem to be processed.");
	}

	try
	{
		// loop ips for both
		for (size_t ip = 0; ip < err_est_data->get(0)->num_elem_ips(elem->reference_object_id()); ip++)
		{
			// loop reactions
			for (size_t j = 0; j < m_reactions.size(); j++)
			{
				const struct ReactionInfo<dim>& r = m_reactions[j];
				if (!r.k_bind.data_given() || !r.k_unbind.data_given() ||  !r.tot_buffer.data_given())
					UG_THROW("Data import for buffering reaction does not have sufficient data.");

				number val_c = 0.0;
				number val_b = 0.0;
				for (size_t sh = 0; sh < m_shapeValues.num_sh(); sh++)
				{
					val_c += u(r.buffered,sh) * m_shapeValues.shapeAtElemIP(sh,ip);
					val_b += u(r.buffer,sh) * m_shapeValues.shapeAtElemIP(sh,ip);
				}

				number val = r.k_bind[ip]*val_b*val_c - r.k_unbind[ip]*(r.tot_buffer[ip] - val_b);

				// add to correct error values
				(*err_est_data->get(this->m_fctGrp[r.buffer])) (elem_list[0],ip) += scale * val;
				(*err_est_data->get(this->m_fctGrp[r.buffered])) (elem_list[0],ip) += scale * val;
			}
		}
	}
	UG_CATCH_THROW("Values for the error estimator could not be assembled at every IP." << std::endl
			<< "Maybe wrong type of ErrEstData object? This implementation needs: SideAndElemErrEstData.");
}

//	postprocesses the loop over all elements of one type in the computation of the error estimator
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void BufferFV1<TDomain>::
fsh_err_est_elem_loop()
{
	// finish the element loop in the same way as the actual discretization
	this->template fsh_elem_loop<TElem, TFVGeom> ();
}

//   error estimation (end)     //
// ///////////////////////////////


////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

// register
template<typename TDomain>
void BufferFV1<TDomain>::
register_all_fv1_funcs()
{
//	get all grid element types in this
	typedef typename domain_traits<dim>::DimElemList ElemList;

//	switch assemble functions
	boost::mpl::for_each<ElemList>(RegisterFV1(this));
}



// explicit template specializations
#ifdef UG_DIM_1
	template class BufferFV1<Domain1d>;
#endif
#ifdef UG_DIM_2
	template class BufferFV1<Domain2d>;
#endif
#ifdef UG_DIM_3
	template class BufferFV1<Domain3d>;
#endif


} // namespace neuro_collection
} // namespace ug
