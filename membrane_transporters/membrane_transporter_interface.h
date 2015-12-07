/*
 *  General channel/pump transport interface class
 *
 *  Created on: 07.01.2015
 *     Authors: mbreit, mstepniewski
 */

#ifndef __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__MEMBRANE_TRANSPORTER_INTERFACE_H__
#define __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__MEMBRANE_TRANSPORTER_INTERFACE_H__


#include "common/common.h"	// for UG_LOG, UG_THROW and others
#include "lib_grid/grid/grid_base_objects.h"
#include "lib_disc/common/local_algebra.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"	// VectorProxyBase

#include <utility>      	// for std::pair
#include <string>
#include <sstream>
#include <vector>
#include <map>



namespace ug{
namespace neuro_collection{


///@addtogroup plugin_neuro_collection
///@{


/**
 * @brief Interface class for membrane transport mechanisms
 *
 * This class is supposed to be the unifying base class for any membrane transport mechanism,
 * i.e. for any type of pumps or channels that transport a quantity (typically ions) through
 * a membrane.
 *
 * It can be used both for unilateral and bilateral exchange mechanisms across a membrane,
 * i.e.
 *   - for situations where both sides of the membrane are being simulated and the flux needs
 *     to be added on one and subtracted on the other side
 *   - for situations where only one side of the membrane is being simulated and on the other
 *     side one has a more or less constant quantity; in such a setting, the flux is only
 *     added or subtracted on the side the simulation is done for and ignored on the other
 *     (infinite reservoir or sink).
 *
 * As a matter of principle, each derivation of this interface should always be able to handle
 * both situations. In the latter case, constant values for the non-simulated side can be set
 * using the set_constant() method.
 * You may, however, forbid certain configurations by implementing the check_constant_allowed()
 * method.
 *
 * @authors mbreit, mstepniewski
 * @date    07.01.2015
 */
class IMembraneTransporter
{
	public:
		/**
		 * @brief Constructor
		 *
		 * Constructs a membrane transporter object that depends or has an effect on each
		 * of the functions given in vFct.
		 * The exact order of the functions is determined by the actual transporter class
		 * that derives from this interface.
		 * If a function the membrane transporter mechanism depends on is not available as
		 * grid function it can be omitted by passing "" (the empty string) as argument in
		 * its place. Note, however, that such missing functions need to be replaced by
		 * constant values via the set_constant() method.
		 *
		 * @param vFct   vector of function names
		 */
		IMembraneTransporter(const std::vector<std::string>& vFct);

		/**
		 * @brief Constructor
		 *
		 * Constructs a membrane transporter object that depends or has an effect on each
		 * of the functions given in vFct.
		 * The exact order of the functions is determined by the actual transporter class
		 * that derives from this interface.
		 * If a function the membrane transporter mechanism depends on is not available as
		 * grid function it can be omitted by passing "" (the empty string) as argument in
		 * its place. Note, however, that such missing functions need to be replaced by
		 * constant values via the set_constant() method.
		 *
		 * @param fct   function names as string, separated by ','
		 */
		IMembraneTransporter(const char* vFct);

		/// Destructor
		virtual ~IMembraneTransporter();

		/**
		 * @brief Prepares a transport mechanism for time step assemblings
		 *
		 * This method will be called by the MembraneTransportFV1 method of the same name.
		 * It can be used if any modifications have to be made before any time step
		 * before the fluxes for this time step can be calculated.
		 * This is especially useful for the implementation of the Borg-Graham type VDCCs using
		 * grid attachments.
		 *
		 * @param time            new point in time
		 * @param u               local vector of unknowns
		 * @param elem            the element modifications can be made for
		 * @param vCornerCoords   the element's corner coordinates
		 */
		virtual void prep_timestep
		(
			const number time, VectorProxyBase* upb
		);

		/**
		 * @brief Calculates the fluxes through this mechanism (same for all mechanisms)
		 *
		 * This method is called by MembraneTransportFV1::fluxDensityFct() and will
		 * receive all of the values this element discretization knows (which are exactly
		 * the given functions from the constructor, as MembraneTransportFV1 is
		 * constructed requiring an object of IMembraneTransporter and gets its functions
		 * from it).
		 * They will be complemented by constant values supplied by the user in such a way
		 * that for any function involved there is a value. These values are then passed to
		 * the calc_flux() method of a derived object which will eventually calculate the
		 * flux.
		 *
		 * This method also takes care of scaling the inputs and output fluxes.
		 *
		 * @param u      vector containing values from known grid functions (at a specific
		 * 				 location)
		 * @param e      element this flux is assembled on
		 * @param flux   output vector containing the calculated fluxes
		 */
		void flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const;

		/**
		 * @brief Calculates the derivatives of the fluxes through this mechanism (same for all mechanisms)
		 *
		 * This method is called by MembraneTransportFV1::fluxDensityDerivFct() and
		 * will receive all of the values this element discretization knows (which are exactly
		 * the given functions from the constructor, as MembraneTransportFV1 is
		 * constructed requiring an object of IMembraneTransporter and gets its functions
		 * from it).
		 * They will be complemented by constant values supplied by the user in such a way
		 * that for any function involved there is a value. These values are then passed to
		 * the calc_flux_deriv() method of a derived object which will eventually calculate the
		 * flux derivatives.
		 *
		 * This method also takes care of scaling the inputs and output flux derivatives.
		 *
		 * @param u             vector containing values from known grid functions (at a
		 *                      specific location)
		 * @param e             element these flux derivatives are assembled on
		 * @param flux_derivs   output matrix containing the calculated flux derivatives as
		 *                      pairs of (index, value), where index specifies with regard
		 *                      to which unknown the derivative is taken and value is the
		 *                      derivative
		 */
		void flux_deriv(const std::vector<number>& u, GridObject* e, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const;

		/**
		 * @brief Calculates the flux through a membrane transport system (system-specific)
		 *
		 * @param u      vector with values for all involved unknowns (created by flux())
		 * @param e      element this flux is assembled on
		 * @param flux   output vector containing the calculated fluxes (not yet scaled)
		 */
		virtual void calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const = 0;

#if 0
		/**
		 * @brief Calculates a specific flux through a membrane transport system
		 *
		 * @param u      vector with values for all involved unknowns (created by flux())
		 * @param index  index of the flux to be calculated
		 * @param flux   output vector containing the calculated fluxes (not yet scaled)
		 */
		virtual number calc_flux(const std::vector<number>& u, size_t index) const {return 0.0;};
#endif

		/**
		 * @brief Calculates the flux derivatives through a membrane transport system (system-specific)
		 *
		 * @param u             vector with values for all involved unknowns (created by flux_deriv())
		 * @param e             element these flux derivatives are assembled on
		 * @param flux_derivs   output matrix containing the calculated flux derivatives as
		 *                      pairs of (index, value), where index specifies with regard
		 *                      to which unknown the derivative is taken and value is the
		 *                      derivative (not yet scaled)
		 */
		virtual void calc_flux_deriv(const std::vector<number>& u, GridObject* e, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const = 0;

		/**
		 * @brief Gives information about how many variables the flux depends on
		 *
		 * That is, the number returned here is the exact number of unknowns that one needs to
		 * calculate derivatives for. Take care not to include the unknowns set constant into
		 * this number!
		 *
		 * @return   number of unknowns this transport mechanism depends on
		 */
		virtual size_t n_dependencies() const = 0;

		/**
		 * @brief Number of fluxes realized by this mechanism
		 * @return   the number of fluxes going through this membrane transport system
		 */
		virtual size_t n_fluxes() const = 0;

		/**
		 * @brief Information on the direction of the i-th flux
		 *
		 * The ordering is the same as for the constructor for parameter i.
		 * The ordering of the return indices is that of m_mfInd.
		 * If one of the indices is to be ignored (in the case of a unilateral flux) this can be
		 * achieved by setting the index to the value InnerBoundaryConstants::_IGNORE_.
		 *
		 * @param i   index the information is to be retrieved for
		 * @return    a pair of indices (from, to), where 'from' is the index of the grid function
		 *            the flux comes from and 'to' is the index of the function it goes to
		 */
		virtual const std::pair<size_t,size_t> flux_from_to(size_t i) const = 0;

		/**
		 * @brief Name of the membrane transport system
		 *
		 * This method can be used by the base class in error outputs in order to specify a concrete
		 * derived class where an error has occurred.
		 *
		 * @return name of the derived mechanism
		 */
		virtual const std::string name() const = 0;

		/**
		 * @brief Return supplied function names
		 *
		 * Supplied functions are those which are not passed to the constructor as "".
		 * This method is called by the constructor MembraneTransportFV1::MembraneTransportFV1().
		 * @return supplied function names
		 */
		const std::vector<std::string>& symb_fcts() const;

		/**
		 * @brief Information on the local function index of the given unknown.
		 *
		 * This method is to be called to retrieve the local function index of an unknown of
		 * a transport mechanism. If none of the unknowns is passed as "" in the constructor
		 * this is the same index as in the constructor.
		 * Otherwise the method will invoke an error if the specified index i does not belong
		 * to a supplied function. If it does belong to one then the index according to the
		 * local map m_mfInd is returned.
		 *
		 * This method is useful when implementing the flux_from_to() method for any transport
		 * mechanism.
		 *
		 * --------
		 * Example
		 * --------
		 * Suppose the membrane transport mechanism has been constructed with {"", "a", "", "b"}.
		 * Then local_fct_index(3) would return 1, since "b" is the second supplied function,
		 * local_fct_index(2) would throw an error as the third unknown is not a supplied function.
		 *
		 * @param i   index (from constructor) which a supplied function index is requested for
		 * @return    local function index for the requested unknown index
		 */
		size_t local_fct_index(const size_t i) const;

		/**
		 * @brief Check whether the combination of supplied functions is allowed
		 *
		 * UG_THROWs, if not allowed (for example when both the inner and outer concentration of the
		 * transported species are not supplied - which would make no sense).
		 */
		virtual void check_supplied_functions() const;

		/**
		 * @brief Set a constant value for one of the unknowns
		 *
		 * This constant will be used in all the calculations instead of any other possibly
		 * existing value.
		 * Note that it is not mandatory for the thereby replaced unknown to have been passed
		 * as "" to the constructor. However, replacing a supplied function by a constant means
		 * that the dependency on this function is effectively eliminated!
		 *
		 * The ordering of indices corresponds to that of the constructor.
		 *
		 * @param i     index of the unknown to be set constant
		 * @param val   constant value to be set
		 */
		void set_constant(const size_t i, const number val);

		/**
		 * @brief Check whether the unknown of an index is set constant
		 *
		 * If it is, the constant value set is returned in the second argument.
		 * Otherwise, the second argument remains unchanged.
		 *
		 * The ordering of indices corresponds to that of the constructor.
		 *
		 * @param i     index to be checked
		 * @param val   constant value (if set)
		 * @return      true iff function of desired index is set constant
		 */
		bool has_constant_value(const size_t i, number& val) const;

		/**
		 * @brief Check whether the unknown of an index is set constant
		 *
		 * The ordering of indices corresponds to that of the constructor.
		 *
		 * @param i     index to be checked
		 * @return      true iff function of specified index is set constant
		 */
		bool has_constant_value(const size_t i) const;

		/**
		 * @brief Check whether the unknown of an index is a supplied function
		 *
		 * Supplied functions are those which are not passed to the constructor as "".
		 * Of course, only supplied functions allow flux being directed to or from them.
		 * Thus, this method allows to perform a check on whether a flux is to be a pure
		 * influx or efflux (without source or sink, resp., e.g. a boundary condition)
		 * or whether it is an exchange flux between to sides.
		 *
		 * The ordering of indices corresponds to that of the constructor.
		 *
		 * @param i   index to be checked
		 * @return    true iff the specified index belongs to a supplied function
		 */
		bool allows_flux(const size_t i) const;

		/**
		 * @brief Prints the units this implementation uses for inputs and flux outputs
		 */
		virtual void print_units() const;

		/**
		 * @brief Scaling of the inputs
		 *
		 * This method can be used to specify scaling factors that are used to adapt the units
		 * of involved functions (which may differ from the units required in the implementation).
		 *
		 * If you want to specify one scaling factor you need to specify them all.
		 * The ordering of indices corresponds to that of the constructor.
		 *
		 * Default values are 1.0 for each factor.
		 *
		 * @param scale   vector of scaling factors
		 */
		void set_scale_inputs(const std::vector<number>& scale);

		/**
		 * @brief Scaling of a single input
		 *
		 * This method can be used to specify a scaling factor for the unit adaptation of a specific
		 * involved function (which may differ from the unit required in the implementation).
		 *
		 * The ordering of indices corresponds to that of the constructor.
		 *
		 * @param i		  index of the input to be scaled
		 * @param scale   scaling factor
		 */
		void set_scale_input(const size_t i, const number scale);

		/**
		 * @brief Scaling of the outputs
		 *
		 * This method can be used to specify scaling factors that are used to adapt the units
		 * of the calculated fluxes (which may differ from the units required in the rest of a
		 * user's discretization).
		 *
		 * If you want to specify one scaling factor you need to specify them all.
		 * The ordering of indices corresponds to that of the constructor.
		 *
		 * @param scale   vector of scaling factors
		 */
		void set_scale_fluxes(const std::vector<number>& scale);

		/**
		 * @brief Scaling of the outputs
		 *
		 * This method can be used to specify a scaling factor for the unit adaptation of a specific
		 * flux calculated by this mechanism (which may differ from the units required in the rest
		 * of a user's discretization).
		 *
		 * The ordering of indices corresponds to that of the constructor.
		 *
		 * @param i		 index of the flux to be scaled
		 * @param scale  scaling factor
		 */
		void set_scale_flux(const size_t i, const number scale);

		/**
		 * @brief Check that all values are either given as unknowns or constants
		 *
		 * If the check is performed with a positive result then the settings will be locked.
		 * This entails that unknowns can no longer be set constant.
		 */
		void check_and_lock();

		/**
		 * @brief Check lock status
		 *
		 * Check if the check_and_lock() method has already locked the settings.
		 *
		 * @return true iff lock has been set
		 */
		bool is_locked() const;

	private:
		/**
		 * @brief Add values set constant to supplied values
		 *
		 * This is a private helper function used in flux() and flux_deriv().
		 *
		 * @param u      vector of given function values
		 * @param u_wc   output vector complemented by constant values
		 *               (or replaced by them if both supplied and constant values exist)
		 */
		void create_local_vector_with_constants(const std::vector<number>& u, std::vector<number>& u_wc) const;

	private:
		/// local vector of supplied function names
		std::vector<std::string> m_vFct;

		/// indices of unknowns in local vector (of supplied functions), if existent
		std::map<size_t, size_t> m_mfInd;

		/// constant values map
		std::map<size_t, number> m_mConstVal;

		/// number of functions in total (supplied or constant)
		const size_t n_fct;

		/// scaling factor for inputs
		std::vector<number> m_vScaleInputs;

		/// scaling factor for calculated fluxes
		std::vector<number> m_vScaleFluxes;

		/// lock status
		bool m_bLocked;
};

///@}

} // namespace neuro_collection
} // namespace ug

#endif // __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__MEMBRANE_TRANSPORTER_INTERFACE_H__

