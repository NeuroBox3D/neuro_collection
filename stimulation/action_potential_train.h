/**
* @file 	action_potential_train.h
*
* @brief 	Header file for the action potential train stimulous.
*
* @author 	Martin Stepniewski <br>
* 			Goethe Center for Scientific Computing GCSC <br>
* 			Kettenhofweg 139 <br>
* 			University Frankfurt <br>
* 			Frankfurt, Germany <br>
* 			email:	martin.stepniewski@gcsc.uni-frankfurt.de <br>
*
* @date		10/06/2015
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



#ifndef __ACTION_POTENTIAL_TRAIN_H__
#define __ACTION_POTENTIAL_TRAIN_H__

#include "common/common.h"
#include <cmath>


namespace ug{
namespace neuro_collection{

///@addtogroup plugin_neuro_collection
///@{

/// Class for action potential train stimulation
/** This class establishes an (isopotential) action potential stimulation
 *  for a given start and end time, AP duration and basic voltage. It
 *  then provides the current membrane potential.
**/
class ActionPotentialTrain
{
	public:
		/// standard constructor
		ActionPotentialTrain();

		/**
		 * @brief constructor with vectors
		 *
		 * @param stimBegin		stimulation begin time in s
		 * @param stimEnd		stimulation end time in s
		 * @param AP_duration	action potential duration in s
		 * @param basicVoltage	resting potential in mV
		 */
		ActionPotentialTrain(	number stimBegin, number stimEnd,
								number stimFreq, number basicVoltage);

		/// destructor
		~ActionPotentialTrain();

		///	Get current membrane potential inside specified frequency stimulation interval
		number membrane_potential(number time);

	private:
		///	Translate global time to local time inside one AP interval
		number global_to_local_AP_time(number time);

		///	characteristic action potential voltage trace
		number AP_voltage_trace(number time);

		number m_stimBegin;
		number m_stimEnd;
		number m_AP_duration;
		number m_basicVoltage;
};

///@}

} // namespace neuro_collection
} // namespace ug

#endif // __ACTION_POTENTIAL_TRAIN_H__

