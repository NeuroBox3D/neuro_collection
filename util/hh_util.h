/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Myra Huymayer
 * Creation date: 2017-09-21
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

#ifndef UG__PLUGINS__NEURO_COLLECTION__UTIL__HH_UTIL_H
#define UG__PLUGINS__NEURO_COLLECTION__UTIL__HH_UTIL_H


#include "common/types.h"  // for number

namespace ug {
namespace neuro_collection {


/// @addtogroup neuro_collection
/// @{

/**
 * @brief n opening rate
 * @param vm membrane potential (in units of mV)
 * @return opening rate (in units of ms^-1)
 */
number alpha_n(number vm);

/**
 * @brief n closing rate
 * @param vm membrane potential (in units of mV)
 * @return closing rate (in units of ms^-1)
 */
number beta_n(number vm);

/**
 * @brief n limit value
 * @param vm membrane potential (in units of mV)
 * @return limit value (in units of 1)
 */
number n_infty(number vm);

/**
 * @brief n time constant
 * @param vm membrane potential (in units of mV)
 * @return time constant (in units of ms)
 */
number tau_n(number vm);


/**
 * @brief m opening rate
 * @param vm membrane potential (in units of mV)
 * @return opening rate (in units of ms^-1)
 */
number alpha_m(number vm);

/**
 * @brief m closing rate
 * @param vm membrane potential (in units of mV)
 * @return closing rate (in units of ms^-1)
 */
number beta_m(number vm);

/**
 * @brief m limit value
 * @param vm membrane potential (in units of mV)
 * @return limit value (in units of 1)
 */
number m_infty(number vm);

/**
 * @brief m time constant
 * @param vm membrane potential (in units of mV)
 * @return time constant (in units of ms)
 */
number tau_m(number vm);


/**
 * @brief h opening rate
 * @param vm membrane potential (in units of mV)
 * @return opening rate (in units of ms^-1)
 */
number alpha_h(number vm);

/**
 * @brief h closing rate
 * @param vm membrane potential (in units of mV)
 * @return closing rate (in units of ms^-1)
 */
number beta_h(number vm);

/**
 * @brief h limit value
 * @param vm membrane potential (in units of mV)
 * @return limit value (in units of 1)
 */
number h_infty(number vm);

/**
 * @brief h time constant
 * @param vm membrane potential (in units of mV)
 * @return time constant (in units of ms)
 */
number tau_h(number vm);

/// @}

} // namspace neuro_collection
} // namespace ug


#endif // UG__PLUGINS__NEURO_COLLECTION__UTIL__HH_UTIL_H
