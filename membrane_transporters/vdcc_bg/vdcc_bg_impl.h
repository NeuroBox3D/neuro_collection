/*
 * vdcc_bg_impl.h
 *
 *  Created on: 05.02.2013
 *      Author: mbreit
 */

#include "vdcc_bg.h"

namespace ug{
namespace neuro_collection{


template<typename TDomain>
template<int TType>
void VDCC_BG<TDomain>::set_channel_type()
{
	if (m_initiated)
	{
		UG_THROW("Borg-Graham channel type can not be changed after initialization.\n"
				 "Call set_channel_type() BEFORE init().");
	}

	switch (TType)
	{
		// - all gating params according to table 5, page 98 of the Borg-Graham article
		//   "Interpretations of data and mechanisms for hippocampal pyramidal cell models"
		// - conductance value lambda of 5.5 pS from Gollasch et al. - "Single calcium
		//   channel currents of arterial smooth muscle at physiological calcium concentrations"
		//	 and used to derive permeability values p_T, p_N and p_L inspired by
		//   Fisher et al. - "Properties and distribution of single
		//   voltage-gated calcium channels in adult hippocampal neurons" in an offset way:
		//
		//	 p_T = 1.9e-19 m^3/s; p_N = 2*p_T; p_L = 3*p_T
		//
		//	 where p = lambda * R * T / (valence^2 * F^2 * Ca_o)
		//
		case BG_Ntype:
			m_gpMGate = GatingParams(3.4, -21.0, 1.5);
			m_gpHGate = GatingParams(-2.0, -40.0, 75.0);
			m_perm = 3.8e-19;
			m_mp = 2;
			m_hp = 1;
			break;
		case BG_Ltype:
			m_gpMGate = GatingParams(4.6, -1.0, 1.5);
			m_gpHGate = GatingParams(0.0, 0.0, 0.0);
			m_perm = 5.7e-19;
			m_mp = 2;
			m_hp = 0;
			break;
		case BG_Ttype:
			m_gpMGate = GatingParams(3.0, -36.0, 1.5);
			m_gpHGate = GatingParams(-5.2, -68.0, 10.0);
			m_perm = 1.9e-19;
			m_mp = 2;
			m_hp = 1;
			break;
		default:
			UG_THROW("Type of Borg-Graham channel does not match any of the pre-implemented.");
			break;
	}
}

} // neuro_collection
} // namespace ug
