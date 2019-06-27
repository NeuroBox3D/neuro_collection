/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Martin Stepniewski
 * Creation date: 2015-06-10
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


/****************************************************************************/
/*                                                                          */
/* Important Channel Parameters:                                            */
/*                                                                          */
/*                              Reverse      absolute						*/
/*								Potential    Conductance       Permeability */
/* Channel    Gating Equation   E_rev (mV)   g (microSiemens)  p (cm^3/s)	*/
/*  Ca_N       m^2*h             135          					1e-8		*/
/*  Ca_L       m^2               135          					3e-9		*/
/*  Ca_T       m^2*h			 135          					3e-9		*/
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/* Important Gating Particle Parameters:                                    */
/*                                                                          */
/*										   particle							*/
/*								effective  voltage							*/
/*								particle   sensor							*/
/*								valence	   asymmetry						*/
/* Channel   Particle   order   z          gamma      K    V_12   tau_0		*/
/*  Ca_N	  m          2        3.4                       -21     1.5		*/
/*			  h          1       -2                         -40    75		*/
/*  Ca_L      m		     2        4.6                       - 1     1.5		*/
/*  Ca_T	  m          2        3                         -36     1.5		*/
/*			  h          1       -5.2                       -68    10		*/
/*                                                                          */
/* _INFO_:                                                                  */
/*                                                                          */
/* K (1/ms)	 :	leading coefficient for alpha_prime_x, beta_prime_x			*/
/* V_12	(mV) :	voltage for which alpha_prime_x and beta_prime_x are equal	*/
/* tau_0 (ms):	rate-limiting step in state transition						*/
/* 1uS		 :  = 1 nA/mV                                                   */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/* Other Parameters:														*/
/*																			*/
/* universal gas constant	:	8.31447 J/(mol*Kelvin) = C*mV/(mmol*Kelvin)	*/
/* Faraday constant			:	9.6485*10^4 C/mol = 96.485 C/mmol			*/
/* absolute temperature		:	300 Kelvin = 27 deg C						*/
/*																			*/
/****************************************************************************/


#include "action_potential_train.h"


namespace ug{
namespace neuro_collection{


ActionPotentialTrain::ActionPotentialTrain()
{
	m_stimBegin = 0.0;
	m_stimEnd = 0.0;
	m_AP_duration = 0.01;
	m_basicVoltage = -65.0;
}


ActionPotentialTrain::ActionPotentialTrain(	number stimBegin, number stimEnd,
											number stimFreq, number basicVoltage)
{
	m_stimBegin = stimBegin;
	m_stimEnd = stimEnd;
	m_AP_duration = 1./stimFreq;
	m_basicVoltage = basicVoltage;

	if(m_AP_duration < 0.01)
			UG_THROW("ActionPotentialTrain: Value for AP_duration only allowed to be >= 0.01.");
}


ActionPotentialTrain::~ActionPotentialTrain()
{
	// do nothing
}


number ActionPotentialTrain::membrane_potential(number time)
{
	double vm;

	if(time < m_stimBegin || time > m_stimEnd)
		vm = m_basicVoltage;
	else
		vm = AP_voltage_trace(time);

	return vm;
}


number ActionPotentialTrain::global_to_local_AP_time(number time)
{
	int nCompletedAPs = time / m_AP_duration;
	return time - m_AP_duration * nCompletedAPs;
}


number ActionPotentialTrain::AP_voltage_trace(number time)
{	
//	Translate global time into time inside one AP interval
	double t = global_to_local_AP_time(time);
	
	/*
	 * De-, re and hyperpolarization always takes place within 0.01s,
	 * AP_duration > 0.010s doesn't stretch the characteristic AP trace,
	 * only the basicVoltage duration and gap to the next AP.
	 *
	 * Avoid AP_duration < 0.010s. AP duration of 0.01s means a max.
	 * stimulation frequency of 100Hz.
	 */

	if(m_AP_duration < 0.01)
		UG_THROW("ActionPotentialTrain: Value for AP_duration only allowed to be >= 0.01.");

	if( t < 0.001 || t > 0.008 )
        return m_basicVoltage;
		
    if( t < 0.0025 ) 
		return m_basicVoltage + 8e3 * ( t - 0.001 );
	
    if( t < 0.003 )
        return -53.0 + 160e3 * ( t - 0.0025 );
	
    if( t < 0.0045 )
        return 27.0 - 62e3 * ( t - 0.003 );
	
    if( t < 0.007 )
        return -66.0 + 0.4e3 * ( t - 0.0045 );

    /*
     * Alternative AP trace
     *
    if( t < 0.001 || t > 0.005 )
            return m_basicVoltage;

	if( t < 0.0015 )
		return m_basicVoltage + 30 * ( t - 0.001 );

	if( t < 0.002 )
		return -55.0 + 190 * ( t - 0.0015 );

	if( t < 0.003 )
		return 40.0 - 120 * ( t - 0.002 );

	if( t < 0.0035 )
		return -80.0 - 20 * ( t - 0.003 );

	if( t < 0.005 )
		return -90.0 + 40./3 * ( t - 0.0035 );
	*/

    return m_basicVoltage;
}


} // neuro_collection
} // namespace ug
























