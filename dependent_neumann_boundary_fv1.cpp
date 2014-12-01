/*
 * dependent_neumann_boundary_fv1.cpp
 *
 *  Created on: 18.01.2013
 *      Author: markusbreit
 */

#include "dependent_neumann_boundary_fv1.h"

namespace ug
{
namespace neuro_collection
{

template<typename TDomain>
void DependentNeumannBoundaryFV1<TDomain>::prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	// check number
	if (vLfeID.size() != this->num_fct())
		UG_THROW("DependentNeumannBoundaryFV1: needs exactly " << this->num_fct() << " functions.");

	// check that Lagrange 1st order
	for (std::size_t i = 0; i < vLfeID.size(); ++i)
		if (vLfeID[i] != LFEID(LFEID::LAGRANGE, dim, 1))
			UG_THROW("DependentNeumannBoundaryFV1: Only first order implemented.");

	// remember
	m_bNonRegularGrid = bNonRegularGrid;

	// update assemble functions
	register_all_fv1_funcs();
}


template<typename TDomain>
bool DependentNeumannBoundaryFV1<TDomain>::use_hanging() const
{
	return true;
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void DependentNeumannBoundaryFV1<TDomain>::prep_elem_loop(const ReferenceObjectID roid, const int si)
{}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void DependentNeumannBoundaryFV1<TDomain>::fsh_elem_loop()
{}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void DependentNeumannBoundaryFV1<TDomain>::prep_elem
(
	const LocalVector& u,
	GridObject* elem,
	const ReferenceObjectID roid,
	const MathVector<dim> vCornerCoords[]
)
{
	// update geometry for this element
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();
	try {
		geo.update(elem, vCornerCoords, &(this->subset_handler()));
	}
	UG_CATCH_THROW("Cannot update Finite Volume Geometry.\n");

	// TODO: The following does not appear to conform to the current paradigm
	// of where geometric data is stored and passed.
	// It works, however. Still, that is...
	// A more conforming way of getting the vertex information to its desired target
	// (e.g. the Borg-Graham implementation) is to be found.
	UG_ASSERT(dynamic_cast<TElem*>(elem) != NULL, "Wrong element type.");
	TElem* pElem = static_cast<TElem*>(elem);

	m_vVertices.clear();

	// resize container
	const std::size_t numVertex = pElem->num_vertices();

	// add vertex pointers to container
	for (std::size_t i = 0; i < numVertex; ++i)
		m_vVertices.push_back(pElem->vertex(i));
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void DependentNeumannBoundaryFV1<TDomain>::add_def_A_elem
(
	LocalVector& d,
	const LocalVector& u,
	GridObject* elem,
	const MathVector<dim> vCornerCoords[]
)
{
	// get finite volume geometry
	const static TFVGeom& fvgeom = GeomProvider<TFVGeom>::get();

	// loop Boundary Faces
	for (std::size_t i = 0; i < fvgeom.num_bf(); ++i)
	{
		// get current BF
		const typename TFVGeom::BF& bf = fvgeom.bf(i);

		// get associated node and subset index
		const int co = bf.node_id();
		int si = fvgeom.subset_index();

		// get solution at the corner of the bf
		std::size_t nFct = u.num_fct();
		std::vector<LocalVector::value_type> uAtCorner = std::vector<LocalVector::value_type>(nFct);
		for (std::size_t fct = 0; fct < nFct; fct++)
			uAtCorner[fct] = u(fct, co);

		// get corner coordinates
		const MathVector<dim>& cc = vCornerCoords[co];

		// store current vertex
		// TODO: This is kind of an ugly quick fix; the vertex is only needed in the BorgGraham class.
		// Maybe make this method virtual here and implement it in BorgGraham for the special case.
		// (The same goes for the derivative and the prep_elem methods.)
		m_currVertex = m_vVertices[co];

		// get flux densities in that node
		NFluxCond fc;
		if (!fluxDensityFct(uAtCorner, elem, cc, si, fc))
			UG_THROW("Call to fluxDensityFct did not succeed.");

		// scale with volume of BF
		for (std::size_t j = 0; j < fc.flux.size(); j++)
			fc.flux[j] *= bf.volume();

		// add to defect
		for (std::size_t j = 0; j < fc.flux.size(); j++)
			d(fc.to[j], co) -= fc.flux[j];
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void DependentNeumannBoundaryFV1<TDomain>::add_def_M_elem
(
	LocalVector& d,
	const LocalVector& u,
	GridObject* elem,
	const MathVector<dim> vCornerCoords[]
)
{}

// assemble stiffness part of Jacobian
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void DependentNeumannBoundaryFV1<TDomain>::add_jac_A_elem
(
	LocalMatrix& J,
	const LocalVector& u,
	GridObject* elem,
	const MathVector<dim> vCornerCoords[]
)
{
	// get finite volume geometry
	const static TFVGeom& fvgeom = GeomProvider<TFVGeom>::get();

	for (std::size_t i = 0; i < fvgeom.num_bf(); ++i)
	{
		// get current BF
		const typename TFVGeom::BF& bf = fvgeom.bf(i);

		// get associated node and subset index
		const int co = bf.node_id();
		int si = fvgeom.subset_index();

		// get solution at the corner of the bf
		std::size_t nFct = u.num_fct();
		std::vector<LocalVector::value_type> uAtCorner = std::vector<
				LocalVector::value_type>(nFct);
		for (std::size_t fct = 0; fct < nFct; fct++)
			uAtCorner[fct] = u(fct, co);

		// get corner coordinates
		const MathVector<dim>& cc = vCornerCoords[co];

		// store current vertex
		m_currVertex = m_vVertices[co];

		NFluxDerivCond fdc;
		if (!fluxDensityDerivFct(uAtCorner, elem, cc, si, fdc))
			UG_THROW("Call to fluxDensityDerivFct did not succeed.");

		// scale with volume of BF
		for (std::size_t j = 0; j < fdc.fluxDeriv.size(); j++)
			for (std::size_t k = 0; k < fdc.fluxDeriv[j].size(); k++)
				fdc.fluxDeriv[j][k] *= bf.volume();

		// add to Jacobian
		for (std::size_t j = 0; j < fdc.fluxDeriv.size(); j++)
			for (std::size_t k = 0; k < fdc.fluxDeriv[j].size(); k++)
				J(fdc.to[j], co, k, co) -= fdc.fluxDeriv[j][k];
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void DependentNeumannBoundaryFV1<TDomain>::add_jac_M_elem
(
	LocalMatrix& J,
	const LocalVector& u,
	GridObject* elem,
	const MathVector<dim> vCornerCoords[]
)
{}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void DependentNeumannBoundaryFV1<TDomain>::add_rhs_elem
(
	LocalVector& d,
	GridObject* elem,
	const MathVector<dim> vCornerCoords[]
)
{}

// ///////////////////////////////
//   error estimation (begin)   //

//	prepares the loop over all elements of one type for the computation of the error estimator
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void DependentNeumannBoundaryFV1<TDomain>::prep_err_est_elem_loop(const ReferenceObjectID roid, const int si)
{
	//	get the error estimator data object and check that it is of the right type
	//	we check this at this point in order to be able to dispense with this check later on
	//	(i.e. in prep_err_est_elem and compute_err_est_A_elem())
	if (this->m_spErrEstData.get() == NULL)
		{UG_THROW("No ErrEstData object has been given to this ElemDisc!");}

	err_est_type* err_est_data = dynamic_cast<err_est_type*>(this->m_spErrEstData.get());

	if (!err_est_data)
	{
		UG_THROW("Dynamic cast to MultipleSideAndElemErrEstData failed." << std::endl
				  << "Make sure you handed the correct type of ErrEstData to this discretization.");
	}

	if (!err_est_data->equal_side_order())
	{
		UG_THROW("The underlying error estimator data objects of this discretization's "
				 "error estimator do not all have the same integration orders. This case "
				 "is not supported by the implementation. If you need it, implement!");
	}

	if (err_est_data->num() < 1)
		{UG_THROW("No underlying error estimator data objects present. No IPs can be determined.");}

	//	set local positions
	if (!TFVGeom::usesHangingNodes)
	{
		static const int refDim = TElem::dim;

		// get local IPs
		std::size_t numSideIPs;
		const MathVector<refDim>* sideIPs;
		try
		{
			numSideIPs = err_est_data->get(0)->num_side_ips(roid);
			sideIPs = err_est_data->get(0)->template side_local_ips<refDim>(roid);
		}
		UG_CATCH_THROW("Integration points for error estimator cannot be set.");

		// store values of shape functions in local IPs
		LagrangeP1<typename reference_element_traits<TElem>::reference_element_type> trialSpace =
			Provider<LagrangeP1<typename reference_element_traits<TElem>::reference_element_type> >::get();

		m_shapeValues.resize(numSideIPs, trialSpace.num_sh());
		for (std::size_t ip = 0; ip < numSideIPs; ip++)
			trialSpace.shapes(m_shapeValues.shapesAtSideIP(ip), sideIPs[ip]);
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void DependentNeumannBoundaryFV1<TDomain>::prep_err_est_elem
(
	const LocalVector& u,
	GridObject* elem,
	const MathVector<dim> vCornerCoords[]
)
{
	// get error estimator
	err_est_type* err_est_data = dynamic_cast<err_est_type*>(this->m_spErrEstData.get());

	// set local positions
	if (TFVGeom::usesHangingNodes)
	{
		static const int refDim = TElem::dim;
		ReferenceObjectID roid = elem->reference_object_id();

		std::size_t numSideIPs;
		const MathVector<refDim>* sideIPs;
		try
		{
			numSideIPs = err_est_data->get(0)->num_side_ips(roid);
			sideIPs = err_est_data->get(0)->template side_local_ips<refDim>(roid);
		}
		UG_CATCH_THROW("Integration points for error estimator cannot be set.");

		// store values of shape functions in local IPs
		LagrangeP1<typename reference_element_traits<TElem>::reference_element_type> trialSpace =
			Provider<LagrangeP1<typename reference_element_traits<TElem>::reference_element_type> >::get();

		m_shapeValues.resize(numSideIPs, trialSpace.num_sh());
		for (std::size_t ip = 0; ip < numSideIPs; ip++)
			trialSpace.shapes(m_shapeValues.shapesAtSideIP(ip), sideIPs[ip]);
	}
}

//	computes the error estimator contribution for one element
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void DependentNeumannBoundaryFV1<TDomain>::compute_err_est_A_elem
(
	const LocalVector& u,
	GridObject* elem,
	const MathVector<dim> vCornerCoords[],
	const number& scale
)
{
	// get error estimator
	err_est_type* err_est_data = dynamic_cast<err_est_type*>(this->m_spErrEstData.get());

	// cast this elem to side_type of error estimator
	typename SideAndElemErrEstData<TDomain>::side_type* side =
		dynamic_cast<typename SideAndElemErrEstData<TDomain>::side_type*>(elem);
	if (!side)
	{
		UG_THROW("Error in DependentNeumannBoundaryFV1<TDomain>::compute_err_est_A_elem():\n"
				 "Element that error assembling routine is called for has the wrong type.");
	}

	// global IPs
	ReferenceObjectID roid = elem->reference_object_id();
	std::size_t numSideIPs = err_est_data->get(0)->num_side_ips(roid);
	MathVector<dim>* globIPs = err_est_data->get(0)->side_global_ips(elem, vCornerCoords);

	// loop IPs
	try
	{
		for (std::size_t sip = 0; sip < numSideIPs; sip++)
		{
			// get values of u at ip (interpolate)
			std::size_t nFct = u.num_fct();
			std::vector<LocalVector::value_type> uAtIP(nFct);

			for (std::size_t fct = 0; fct < nFct; fct++)
			{
				uAtIP[fct] = 0.0;
				for (std::size_t sh = 0; sh < m_shapeValues.num_sh(); sh++)
					uAtIP[fct] += u(fct, sh) * m_shapeValues.shapeAtSideIP(sh, sip);
			}

			// ip coordinates
			const MathVector<dim>& ipCoords = globIPs[sip];

			// elem subset
			int si = this->subset_handler().get_subset_index(elem);

			NFluxCond fc;
			if (!fluxDensityFct(uAtIP, elem, ipCoords, si, fc))
			{
				UG_THROW("DependentNeumannBoundaryFV1::compute_err_est_A_elem:"
						" Call to fluxDensityFct did not succeed.");
			}

			// add to estimator values
			for (std::size_t j = 0; j < fc.flux.size(); j++)
				(*err_est_data->get(fc.to[j]))(side, sip) += scale * fc.flux[j];
		}
	}
	UG_CATCH_THROW("Values for the error estimator could not be assembled at every IP."
		<< std::endl << "Maybe wrong type of ErrEstData object? "
		"This implementation needs: SideAndElemErrEstData.");
}

//	postprocesses the loop over all elements of one type in the computation of the error estimator
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void DependentNeumannBoundaryFV1<TDomain>::fsh_err_est_elem_loop()
{
	// nothing to do
}

//   error estimation (end)     //
// ///////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

// register for 1D
template<typename TDomain>
void DependentNeumannBoundaryFV1<TDomain>::register_all_fv1_funcs()
{
//	get all grid element types in this dimension and below
	typedef typename domain_traits<dim>::ManifoldElemList ElemList;

//	switch assemble functions
	boost::mpl::for_each<ElemList>(RegisterFV1(this));
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void DependentNeumannBoundaryFV1<TDomain>::register_fv1_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	this->clear_add_fct(id);
	this->set_prep_elem_loop_fct(id, &T::template prep_elem_loop<TElem, TFVGeom>);
	this->set_prep_elem_fct(id, &T::template prep_elem<TElem, TFVGeom>);
	this->set_fsh_elem_loop_fct(id, &T::template fsh_elem_loop<TElem, TFVGeom>);
	this->set_add_jac_A_elem_fct(id, &T::template add_jac_A_elem<TElem, TFVGeom>);
	this->set_add_jac_M_elem_fct(id, &T::template add_jac_M_elem<TElem, TFVGeom>);
	this->set_add_def_A_elem_fct(id, &T::template add_def_A_elem<TElem, TFVGeom>);
	this->set_add_def_M_elem_fct(id, &T::template add_def_M_elem<TElem, TFVGeom>);
	this->set_add_rhs_elem_fct(id, &T::template add_rhs_elem<TElem, TFVGeom>);

	// error estimator parts
	this->set_prep_err_est_elem_loop(id, &T::template prep_err_est_elem_loop<TElem, TFVGeom>);
	this->set_prep_err_est_elem(id, &T::template prep_err_est_elem<TElem, TFVGeom>);
	this->set_compute_err_est_A_elem(id, &T::template compute_err_est_A_elem<TElem, TFVGeom>);
	this->set_fsh_err_est_elem_loop(id, &T::template fsh_err_est_elem_loop<TElem, TFVGeom>);
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class DependentNeumannBoundaryFV1<Domain1d>;
#endif
#ifdef UG_DIM_2
template class DependentNeumannBoundaryFV1<Domain2d>;
#endif
#ifdef UG_DIM_3
template class DependentNeumannBoundaryFV1<Domain3d>;
#endif


} // end namespace neuro_collection
} // end namspace ug
