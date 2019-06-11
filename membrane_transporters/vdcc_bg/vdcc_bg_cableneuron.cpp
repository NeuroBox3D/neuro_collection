/*
 * vdcc_bg_cableneuron.cpp
 *
 *  Created on: 20.12.2016
 *      Author: mbreit
 */

#include "vdcc_bg_cableneuron.h"

#include "lib_grid/algorithms/debug_util.h" // ElementDebugInfo
#include "lib_disc/operator/composite_conv_check.h"
#include "lib_algebra/operator/operator_util.h" // ApplyLinearSolver
#include "lib_disc/function_spaces/interpolate.h"

namespace ug {
namespace neuro_collection {


template <typename TDomain>
VDCC_BG_CN<TDomain>::VDCC_BG_CN
(
    const std::vector<std::string>& fcts,
    const std::vector<std::string>& subsets,
    SmartPtr<ApproximationSpace<TDomain> > approx1d,
    SmartPtr<ApproximationSpace<TDomain> > approx3d,
    const std::string& potFctName
)
: VDCC_BG<TDomain>(fcts, subsets, approx3d),
  m_spApprox1d(approx1d), m_spApprox3d(approx3d),
  m_spHNC(SPNULL), m_spDomDiscCable(SPNULL), m_spCableDisc(SPNULL),
  m_spAssTuner(SPNULL), m_spTimeDisc(SPNULL), m_spLinOp(SPNULL),
  m_spILU(SPNULL), m_spLinearSolver(SPNULL),
  m_spU(SPNULL), m_spUOld(SPNULL), m_spB(SPNULL),
  m_spSolTimeSeries(SPNULL), m_spVTKOut(SPNULL),
  m_potFctInd(0),
  m_scaleFactor3dTo1d(1.0),
  m_bSolverVerboseOutput(false), m_bVTKOutput(false),
  m_vtkFileName(std::string("")), m_pstep(1e-3),
  m_vNID(1,0),
  m_curTime(0.0), m_dt(1e-4), m_stepLv(0),
  m_dt_potentialUpdate(1e-4), m_timeSinceLastPotentialUpdate(0.0)
{
    m_potFctInd = approx1d->fct_id_by_name(potFctName.c_str());
}


template <typename TDomain>
VDCC_BG_CN<TDomain>::~VDCC_BG_CN()
{}


template <typename TDomain>
void VDCC_BG_CN<TDomain>::set_domain_disc_1d(SmartPtr<DomainDiscretization<TDomain, CPUAlgebra> > domDisc)
{
    m_spDomDiscCable = domDisc;
}


template <typename TDomain>
void VDCC_BG_CN<TDomain>::set_cable_disc(SmartPtr<cable_neuron::CableEquation<TDomain> > ce)
{
    m_spCableDisc = ce;
}


template <typename TDomain>
void VDCC_BG_CN<TDomain>::set_3d_neuron_ids(const std::vector<size_t>& ids)
{
	// uint is not registered, we therefore use size_t as param type
	// but we need to convert this here
	m_vNID.resize(ids.size());
	for (size_t i = 0; i < ids.size(); ++i)
		m_vNID[i] = (uint) ids[i];

    if (m_spHNC.valid())
    	m_spHNC->set_neuron_ids(m_vNID);
}

#if 0
// unused
template <typename TDomain>
void VDCC_BG_CN<TDomain>::set_hybrid_neuron_communicator(SmartPtr<HybridNeuronCommunicator<TDomain> > spHNC)
{
    m_spHNC = spHNC;
    m_spHNC->set_potential_subsets(this->m_vSubset);
}
#endif

template <typename TDomain>
void VDCC_BG_CN<TDomain>::set_coordinate_scale_factor_3d_to_1d(number scale)
{
    m_scaleFactor3dTo1d = scale;
}


template <typename TDomain>
void VDCC_BG_CN<TDomain>::set_solver_output_verbose(bool verbose)
{
    m_bSolverVerboseOutput = verbose;
}


template <typename TDomain>
void VDCC_BG_CN<TDomain>::set_vtk_output(const std::string& fileName, number plotStep)
{
    m_bVTKOutput = true;
    m_vtkFileName = fileName;
    m_pstep = plotStep;
}


template <typename TDomain>
void VDCC_BG_CN<TDomain>::set_initial_values(const std::vector<number>& v_initVals)
{
    // check that number of initial values matches number of functions in approx space
    UG_COND_THROW(m_spApprox1d->dd(GridLevel(), false)->num_fct() != v_initVals.size(),
        "Number of provided initial values (" << v_initVals.size() << ") does not match "
        "number of functions contained in approximation space (" << m_spApprox1d->dd(GridLevel(), false)->num_fct() << ").");

    m_vInitVals = v_initVals;
}


template <typename TDomain>
void VDCC_BG_CN<TDomain>::set_time_steps_for_simulation_and_potential_update(number dtSim, number dtPot)
{
    m_dt = dtSim;
    m_dt_potentialUpdate = dtPot;

    if (m_dt_potentialUpdate < m_dt)
    {
        UG_LOGN("Warning: Time step size chosen for potential update is smaller than that chosen for the cable simulation.\n"
            "This is not allowed. Increasing potential update interval.");
        m_dt_potentialUpdate = m_dt;
    }
}



template <typename TDomain>
void VDCC_BG_CN<TDomain>::init(number time)
{
    // check that all necessary data inputs have been set
    UG_COND_THROW(!m_spDomDiscCable.valid(), "No valid 1d domain disc set. "
        "Do so using the method 'set_domain_disc_1d'.");
    UG_COND_THROW(!m_vInitVals.size(), "No initial values set. Do so by using the method 'set_initial_values'.");
    UG_COND_THROW(!m_spCableDisc.valid(), "No valid 1d cable discretization set. "
        "Do so by using the method 'set_cable_disc'.");

    // save ass tuner for later
    m_spAssTuner = m_spDomDiscCable->ass_tuner();

    // time discretization
    m_spTimeDisc = make_sp(new ThetaTimeStep<algebra_t>(m_spDomDiscCable));
    m_spTimeDisc->set_theta(1.0);

    // create instationary operator
    m_spLinOp = make_sp(new AssembledLinearOperator<algebra_t>(m_spTimeDisc));

    // solver setup
    //SmartPtr<GridFunctionDebugWriter<TDomain, algebra_t> > dbgWriter
    //   = make_sp(new GridFunctionDebugWriter<TDomain, algebra_t>(m_spApprox1d));
    //dbgWriter->set_vtk_output(true);

    // todo: maybe have the user set a convergence check if needed for more flexibility
    SmartPtr<CompositeConvCheck<vec_t, TDomain> > linConvCheck
       (new CompositeConvCheck<vec_t, TDomain>(m_spApprox1d, 2, 1e-21, 1e-12));
    //linConvCheck->set_component_check("v", 1e-21, 1e-12);
    linConvCheck->set_verbose(m_bSolverVerboseOutput);

    m_spILU = make_sp(new ILU<algebra_t>());
    m_spLinearSolver = make_sp(new LinearSolver<vec_t>());
    m_spLinearSolver->set_preconditioner(m_spILU);
    m_spLinearSolver->set_convergence_check(linConvCheck);
    //cgSolver:set_debug(dbgWriter)


    // initialize solution and time stepping
    m_curTime = time;

    m_spU = make_sp(new GridFunction<TDomain, algebra_t>(m_spApprox1d));
    m_spB = make_sp(new GridFunction<TDomain, algebra_t>(m_spApprox1d));
    m_spU->set(0.0);

    size_t nFct = m_spApprox1d->dd(GridLevel(), false)->num_fct();
    for (size_t fct = 0; fct < nFct; ++fct)
    {
        const std::string& fctName = m_spU->name(fct);
        Interpolate(m_vInitVals[fct], m_spU, fctName.c_str());
    }

    // write start solution to vtk file if required
    if (m_bVTKOutput)
    {
        m_spVTKOut = make_sp(new VTKOutput<TDomain::dim>());
        try {m_spVTKOut->print(m_vtkFileName.c_str(), *m_spU, 0, time);}
        UG_CATCH_THROW("VTK output file prefixed '" << m_vtkFileName << "' could not be written to.");
    }

    // store grid function in vector of old solutions
    m_spUOld = m_spU->clone();
    m_spSolTimeSeries = make_sp(new VectorTimeSeries<vec_t>());
    m_spSolTimeSeries->push(m_spUOld, time);

    // prepare step size adaptivity
    m_stepLv = 0;
    m_StepCheckBackCounter[m_stepLv] = 0;



    // finally, init hybrid neuron communicator ...
    m_spHNC = make_sp(new HybridNeuronCommunicator<TDomain>(m_spApprox3d, m_spApprox1d));
    m_spHNC->set_coordinate_scale_factor_3d_to_1d(m_scaleFactor3dTo1d);
    m_spHNC->set_potential_subsets(this->m_vSubset);
    m_spHNC->set_solution_and_potential_index(m_spU, m_potFctInd);
	m_spHNC->set_neuron_ids(m_vNID);

    // ... and communicate current membrane potential values
    m_spHNC->coordinate_potential_values();

    // call init of base class
    VDCC_BG<TDomain>::init(time);
}



template <typename TDomain>
void VDCC_BG_CN<TDomain>::prep_timestep(number future_time, const number time, VectorProxyBase* upb)
{
    // initiate if this has not already been done
    if (!this->m_initiated)
        init(time);

    // do nothing if new time is older than we are
    // (might happen in case of external reset due to non-convergence)
    if (future_time <= m_curTime)
    {
        // we may need to update the potential mapping (after redistribution)
        m_spHNC->coordinate_potential_values();

        return;
    }

    // the 1d simulation needs to be updated to the given time
    // modify time step size if needed
    number nsteps = floor((future_time - m_curTime) / m_dt);
    if (future_time - (m_curTime + nsteps*m_dt) > 1e-3)
        m_dt = (future_time - m_curTime) / nsteps;

    while (future_time - m_curTime > 1e-4*m_dt)
    {
        // setup time disc for old solutions and time step
        m_spTimeDisc->prepare_step(m_spSolTimeSeries, m_dt);

        // reduce time step if cfl < m_dt
        // (this needs to be done AFTER prepare_step as channels are updated there)
        bool dtChanged = false;
        number cfl = m_spCableDisc->template estimate_cfl_cond<vec_t>(m_spSolTimeSeries->latest());
        if (m_bSolverVerboseOutput)
            UG_LOGN("estimated CFL condition: dt < " << cfl)
        while (m_dt > cfl)
        {
            m_dt = m_dt/2.0;

            UG_COND_THROW(m_stepLv+1 > 15, "Time step for 1d cable simulation too small.");

            ++m_stepLv;
            m_StepCheckBackCounter[m_stepLv] = 0;
            if (m_bSolverVerboseOutput)
                UG_LOGN("estimated CFL condition: dt < " << cfl << " - reducing time step to " << m_dt);
            dtChanged = true;
        }

        // increase time step if cfl > m_dt / 2.0 (and if time is aligned with new bigger step size)
        while (m_dt*2.0 < cfl && m_stepLv > 0 && m_StepCheckBackCounter[m_stepLv] % 2 == 0)
        {
            m_dt *= 2.0;
            --m_stepLv;
            m_StepCheckBackCounter[m_stepLv] += m_StepCheckBackCounter[m_stepLv+1]/2.0;
            m_StepCheckBackCounter[m_stepLv+1] = 0;
            if (m_bSolverVerboseOutput)
                UG_LOGN("estimated CFL condition: dt < " << cfl << " - increasing time step to " << m_dt);
            dtChanged = true;
        }

        if (m_bSolverVerboseOutput)
            UG_LOGN("++++++ POINT IN TIME " << floor((m_curTime+m_dt)/m_dt+0.5)*m_dt << " BEGIN ++++++");

        // prepare again with new time step size
        if (dtChanged)
            m_spTimeDisc->prepare_step(m_spSolTimeSeries, m_dt);

        // assemble linear problem
        bool matrixIsConst = m_curTime != 0.0 && !dtChanged; // TODO: what if start time is != 0 !?
        m_spAssTuner->set_matrix_is_const(matrixIsConst);
        AssembleLinearOperatorRhsAndSolution(*m_spLinOp, *m_spU, *m_spB);

        // apply linear solver
        m_spILU->set_disable_preprocessing(matrixIsConst);
        UG_COND_THROW(!ApplyLinearSolver<vec_t>(m_spLinOp, *m_spU, *m_spB, m_spLinearSolver),
            "Could not apply linear solver for the 1d cable problem.");

        // update to new time
        m_curTime = m_spSolTimeSeries->time(0) + m_dt;
        m_timeSinceLastPotentialUpdate += m_dt;

        // update potential values and integrate gating params
        if (m_timeSinceLastPotentialUpdate >= m_dt_potentialUpdate)
        {
            // communicate current membrane potential values
            m_spHNC->coordinate_potential_values();

            // call regular prep timestep from base class
            VDCC_BG<TDomain>::prep_timestep(future_time, time, upb);

            m_timeSinceLastPotentialUpdate -= m_dt_potentialUpdate;
        }

        // vtk output, if required
        if (m_bVTKOutput && fabs(m_curTime / m_pstep - floor(m_curTime / m_pstep + 0.5)) < 1e-5)
            m_spVTKOut->print(m_vtkFileName.c_str(), *m_spU, floor(m_curTime / m_pstep + 0.5), m_curTime);

        // update time series (reuse memory)
        m_spUOld = m_spSolTimeSeries->oldest();
#ifdef UG_PARALLEL
        VecAssign(*m_spUOld, *m_spU);
#else
        *m_spUOld = *m_spU;
#endif
        m_spSolTimeSeries->push_discard_oldest(m_spUOld, m_curTime);

        // increment check-back counter
        ++m_StepCheckBackCounter[m_stepLv];

        if (m_bSolverVerboseOutput)
            UG_LOGN("++++++ POINT IN TIME " << floor(m_curTime / m_dt + 0.5) * m_dt << "  END ++++++");
    }

    // end timeseries, produce gathering file if required
    if (m_bVTKOutput)
        m_spVTKOut->write_time_pvd(m_vtkFileName.c_str(), *m_spU);

}


template <typename TDomain>
void VDCC_BG_CN<TDomain>::update_potential(side_t* elem)
{
    // fill attachments with renewed values
    this->m_aaVm[elem] = m_spHNC->potential(elem);
}


template <typename TDomain>
const std::string VDCC_BG_CN<TDomain>::name() const
{
    return std::string("VDCC_BG_CN");
}


template <typename TDomain>
void VDCC_BG_CN<TDomain>::print_units() const
{
    std::string nm = name();
    size_t n = nm.size();
    UG_LOG(std::endl);
    UG_LOG("+------------------------------------------------------------------------------+"<< std::endl);
    UG_LOG("|  Units used in the implementation of " << nm << std::string(n>=40?0:40-n, ' ') << "|" << std::endl);
    UG_LOG("|------------------------------------------------------------------------------|"<< std::endl);
    UG_LOG("|    Input                                                                     |"<< std::endl);
    UG_LOG("|      [Ca_cyt]  mM (= mol/m^3)                                                |"<< std::endl);
    UG_LOG("|      [Ca_ext]  mM (= mol/m^3)                                                |"<< std::endl);
    UG_LOG("|      V_m       V                                                             |"<< std::endl);
    UG_LOG("|                                                                              |"<< std::endl);
    UG_LOG("|    Output                                                                    |"<< std::endl);
    UG_LOG("|      Ca flux   mol/s                                                         |"<< std::endl);
    UG_LOG("+------------------------------------------------------------------------------+"<< std::endl);
    UG_LOG(std::endl);
}


// explicit template specializations
#ifdef UG_DIM_1
    template class VDCC_BG_CN<Domain1d>;
#endif
#ifdef UG_DIM_2
    template class VDCC_BG_CN<Domain2d>;
#endif
#ifdef UG_DIM_3
    template class VDCC_BG_CN<Domain3d>;
#endif


} // namespace neuro_collection
} // namespace ug
