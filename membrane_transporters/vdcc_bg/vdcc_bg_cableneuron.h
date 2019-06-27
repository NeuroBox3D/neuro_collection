/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2016-12-20
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

#ifndef UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__VDCC_BG__VDCC_BG_CABLENEURON_H
#define UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__VDCC_BG__VDCC_BG_CABLENEURON_H

#include "vdcc_bg.h"
#include "../../hybrid_neuron_communicator.h"

#include "lib_disc/time_disc/theta_time_step.h"
#include "lib_disc/operator/linear_operator/assembled_linear_operator.h"
#include "lib_algebra/operator/preconditioner/ilu.h"
#include "lib_algebra/operator/linear_solver/linear_solver.h"
#include "lib_disc/io/vtkoutput.h"

namespace ug {
namespace neuro_collection {


/// @addtogroup plugin_neuro_collection
/// @{

/// Borg Graham type VGCCs with membrane potential supply by UG's cable_neuron project.
/** This class is a specialization of the Borg-Graham interface.
 *  It supplies the channel with the necessary membrane potential values by an object of
 *  cable_neuron::CableEquation object discretizing a 1d representation of a 3d geometry.
**/

template <typename TDomain>
class VDCC_BG_CN
: public VDCC_BG<TDomain>
{
    protected:
        typedef typename GeomObjBaseTypeByDim<TDomain::dim>::base_obj_type elem_t;
        typedef typename elem_t::side side_t;
        typedef CPUAlgebra algebra_t;
        typedef algebra_t::vector_type vec_t;

    public:
        /// constructor
        VDCC_BG_CN
        (
            const std::vector<std::string>& fcts,
            const std::vector<std::string>& subsets,
            SmartPtr<ApproximationSpace<TDomain> > approx1d,
            SmartPtr<ApproximationSpace<TDomain> > approx3d,
            const std::string& potFctName
        );

        /// destructor
        virtual ~VDCC_BG_CN();

        // the following methods MUST be called before any object of this class can work properly
        /// set the domain disc for the cable problem
        void set_domain_disc_1d(SmartPtr<DomainDiscretization<TDomain, CPUAlgebra> > domDisc);

        /// set the elem disc for the cable problem
        void set_cable_disc(SmartPtr<cable_neuron::CableEquation<TDomain> > ce);

        /// set the neuron IDs of the 3d represented neuron within the network
        void set_3d_neuron_ids(const std::vector<size_t>& ids);

        /// set equilibrium potential
        void set_initial_values(const std::vector<number>& v_initVals);

        // the remaining methods are optional
        /**
         * @brief Set scale factor from 3d to 1d.
         * This is required as 3d and 1d morphology might be differently scaled.
         * For example, 1d geometries used in the cable_neuron code should be scaled to length unit [m],
         * while 3d simulation length scale is often [um]. In that case, we need a scaling factor of 1e-6.
         */
        void set_coordinate_scale_factor_3d_to_1d(number scale);

        /// set whether solver output is to be verbose
        void set_solver_output_verbose(bool verbose);

        /// set vtk output is required and where and when to write to
        void set_vtk_output(const std::string& fileName, number plotStep);

        /// set time step
        void set_time_steps_for_simulation_and_potential_update(number dtSim, number dtPot);

        /// set a communicator object for hybrid neuron treatment
        //unused
        //void set_hybrid_neuron_communicator(SmartPtr<HybridNeuronCommunicator<TDomain> > spHNC);


        /// @copydoc VDCC_BG::init()
        virtual void init(number time);

        /// @copydoc IMembraneTransporter::prepare_timestep()
        virtual void prepare_timestep(number future_time, const number time, VectorProxyBase* upb);

        /// @copydoc VDCC_BG::update_potential()
        virtual void update_potential(side_t* elem);

        /// @copydoc IMembraneTransporter::name()
        virtual const std::string name() const;

        /// @copydoc IMembraneTransporter::print_units()
        virtual void print_units() const;

    private:
        SmartPtr<ApproximationSpace<TDomain> > m_spApprox1d;
        SmartPtr<ApproximationSpace<TDomain> > m_spApprox3d;
        SmartPtr<HybridNeuronCommunicator<TDomain> > m_spHNC;
        SmartPtr<DomainDiscretization<TDomain, algebra_t> > m_spDomDiscCable;
        SmartPtr<cable_neuron::CableEquation<TDomain> > m_spCableDisc;
        SmartPtr<AssemblingTuner<algebra_t> > m_spAssTuner;
        SmartPtr<ThetaTimeStep<algebra_t> > m_spTimeDisc;
        SmartPtr<AssembledLinearOperator<algebra_t> > m_spLinOp;
        SmartPtr<ILU<algebra_t> > m_spILU;
        SmartPtr<LinearSolver<vec_t> > m_spLinearSolver;
        SmartPtr<GridFunction<TDomain, algebra_t> > m_spU;
        SmartPtr<vec_t> m_spUOld;
        SmartPtr<GridFunction<TDomain, algebra_t> > m_spB;
        SmartPtr<VectorTimeSeries<vec_t> > m_spSolTimeSeries;
        SmartPtr<VTKOutput<TDomain::dim> > m_spVTKOut;

        size_t m_potFctInd;
        number m_scaleFactor3dTo1d;
        bool m_bSolverVerboseOutput;
        bool m_bVTKOutput;
        std::string m_vtkFileName;
        number m_pstep;
        std::vector<number> m_vInitVals;

        std::vector<uint> m_vNID;

        number m_curTime;
        number m_dt;
        size_t m_stepLv;
        size_t m_StepCheckBackCounter[16];

        number m_dt_potentialUpdate;
        number m_timeSinceLastPotentialUpdate;
};

} // namespace neuro_collection
} // namespace ug

#endif // UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__VDCC_BG__VDCC_BG_CABLENEURON_H
