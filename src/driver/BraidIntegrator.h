//
// Created by parnet on 24.05.21.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_BRAIDINTEGRATOR_H
#define UG_PLUGIN_XBRAIDFORUG4_BRAIDINTEGRATOR_H


#include "BraidGridFunctionBase.h"


template<typename TDomain, typename TAlgebra>
class LevelIntegrator {
public: // todo set better modes

    typedef ug::ITimeIntegrator<TDomain, TAlgebra> TTimeIntegrator;
    typedef SmartPtr<TTimeIntegrator> SPTimeIntegrator;

    SPTimeIntegrator m_time_integrator;

    void set_time_integrator(SPTimeIntegrator integrator){
        m_time_integrator = integrator;
    }

    SPTimeIntegrator get_time_integrator(){
        return m_time_integrator;
    }
};




template<typename TDomain, typename TAlgebra>
class BraidIntegrator : public BraidGridFunctionBase<TDomain,TAlgebra> {
public: // todo set better modes
    /* *************************************************************************************
     * Type definitions
     ************************************************************************************* */
    typedef typename TAlgebra::matrix_type TMatrix;
    typedef typename TAlgebra::vector_type TVector;
    typedef ug::GridFunction<TDomain, TAlgebra> TGridFunction;
    typedef SmartPtr<TGridFunction> SPGridFunction;


    typedef ug::ITimeIntegrator<TDomain, TAlgebra> TTimeIntegrator;
    typedef SmartPtr<TTimeIntegrator> SPTimeIntegrator;
    typedef std::vector<LevelIntegrator<TDomain,TAlgebra>> LevelIntegratorList;

    typedef ug::INonlinearTimeIntegrator<TDomain, TAlgebra> TBaseType; // todo derivate ?
    typedef ug::INonlinearTimeIntegrator<TDomain, TAlgebra> TINonlinearTimeIntegrator;
    typedef typename TBaseType::solver_type TSolverType;
    typedef SmartPtr<TSolverType> SPSolverType;

    typedef ug::IDomainDiscretization<TAlgebra>	TDomainDisc;
    typedef SmartPtr<TDomainDisc> SPDomainDisc;

    typedef ug::SimpleTimeIntegrator<TDomain, TAlgebra> TTimeIntegratorType;
    typedef ug::LinearImplicitEuler<TAlgebra> TTimeStepType;
    typedef SmartPtr<TTimeStepType > SPTimeStepType;
    typedef ug::ISubDiagErrorEst<TVector> TErrorEstimatorType;
    typedef ug::AitkenNevilleTimex<TVector> TTimexType;



    typedef ug::StdConvCheck<typename TAlgebra::vector_type> TConv;
    typedef SmartPtr<TConv> SPConv;

    typedef SmartPtr <SpaceTimeCommunicator> SPCommunicator;

    typedef Paralog TParalog;
    typedef SmartPtr<TParalog> SPParalog;

private:
    /* -------------------------------------------------------------------------------------
     * Member variables
     * ---------------------------------------------------------------------------------- */

private:
    SPTimeIntegrator  m_default_integrator;
    LevelIntegratorList m_level_integrator;

    // bool m_verbose = true; // todo set to false
    //bool m_timing = true;
    //bool m_writeparam = true;

    //size_t m_max_levels; // todo set get
    //bool m_force_conv = true; // todo set get
    //bool m_adapt_conv_check = false; // todo set_get // todo outsource?

    //int m_levels = 15;

    //int m_progress = 0;


    //SPCommunicator m_comm;


public:
    SPTimeIntegrator m_time_integrator;
    /**
    * Note that this default constructor does not create a consistent object. The parameter t_comm (of type MPI_Comm)
    * for the temporal communication has to be set.
    */
    BraidIntegrator() : BraidGridFunctionBase<TDomain,TAlgebra>() {
        //this->m_log->o << "debug::BraidIntegrator()" << std::endl<<std::flush;
        this->m_name = "Braid Integrator";
        //this->m_log->o << "Warning Braid Integrator :: default constructor was used!" << std::endl;
    }
    BraidIntegrator(MPI_Comm mpi_temporal, double tstart, double tstop, int steps) : BraidGridFunctionBase<TDomain,TAlgebra>(mpi_temporal, tstart,tstop, steps) {
        //this->m_log->o << "debug::BraidIntegrator(args)" << std::endl<<std::flush;
        this->m_name = "Braid Integrator";
    }
    ~BraidIntegrator() = default;
    /* -------------------------------------------------------------------------------------
     * XBraid Method definitions
     * ---------------------------------------------------------------------------------- */
    braid_Int Step(braid_Vector u_,braid_Vector ustop_,braid_Vector fstop_,BraidStepStatus &pstatus) override {
        //this->m_log->o << "debug::BraidIntegrator::Step[[args]]" << std::endl<<std::flush;
        //print_status(this->m_log->o,pstatus);

#if TRACE_CONST == 1
        u->m_const = false;
//ustop->m_const = false;
#endif
        int l; // level
        pstatus.GetLevel(&l);
        StartLevelOperationTimer(LevelObserver::TL_STEP, l);
        double t_start, t_stop;
        pstatus.GetTstartTstop(&t_start, &t_stop);

        double current_dt = t_stop - t_start;

// this->m_log->o << "num stages " << this->m_timeDisc->num_stages() << std::endl;

//this->m_log->o << "message" << std::endl;
        if (this->m_verbose) {
            int tindex;
            pstatus.GetTIndex(&tindex);
            int iteration;
            pstatus.GetIter(&iteration);
#if TRACE_INDEX == 1


            if (fstop_ == nullptr) {
                this->m_log->o << "u_" << u_->index << " = step_"<<l<<"_n( u_" << u_->index << ", u_" << ustop_->index
                                  << ", null, " << t_start << ", " << current_dt << ", " << t_stop << ", " << l
                                  << ")"
                                  << "\t\t % " << tindex << std::endl<< std::flush;
            } else {
                this->m_log->o << "u_" << u_->index << " = step_"<<l<<"_r( u_" << u_->index << ", u_" << ustop_->index
                                  << ", u_"
                                  << fstop_->index << ", " << t_start << ", " << current_dt << ", " << t_stop << ", "
                                  << l << ")"
                                  << " % " << tindex << std::endl << std::flush;
            }

#else
            this->m_log->o << std::setw(13) << iteration << "step for level " << l << " at position " << tindex << " and iteration" << std::endl;
#endif
        }
//this->m_log->o << "preparation" << std::endl;
        auto *sp_u_approx_tstart = (SPGridFunction *) u_->value;
        auto *constsp_u_approx_tstop = (SPGridFunction *) ustop_->value;
        SPGridFunction sp_u_tstop_approx = constsp_u_approx_tstop->get()->clone();
        SPGridFunction lp = constsp_u_approx_tstop->get()->clone();
        SPGridFunction sp_rhs = this->m_u0->clone_without_values(); // for rhs
        //const ug::GridLevel gridlevel = sp_u_approx_tstart->get()->grid_level();
        // todo adapt conv check?
        if (fstop_ != nullptr) {
            //this->m_log->o << "Warning residul is ignored" << std::endl<< std::flush;
        } else {
            //this->m_log->o << "residual not used " << std::endl<< std::flush;
        }


//this->m_log->o << "fstop" << std::endl;


//this->m_log->o << "solve" << std::endl;
        StartLevelOperationTimer(LevelObserver::TL_SOLVE, l);
        // smartptr, const smartptr

        //this->m_log->o << "integrator application is comeing " << std::endl<< std::flush;

        SPTimeIntegrator loc_time_integrator = m_time_integrator;
        m_time_integrator->set_time_step(current_dt);
        bool success = m_time_integrator->apply(sp_u_tstop_approx, t_stop, sp_u_approx_tstart->cast_const(), t_start);
        //bool success = true;
        //this->m_log->o << std::flush;
        //this->m_log->o << "integrator applied "<< std::endl << std::flush;

        //this->m_log->o << "integrator applied "<< std::endl << std::flush;

        StopLevelOperationTimer(LevelObserver::TL_SOLVE, l);
        if (!success ){//&& m_force_conv) {
            this->m_log->o << "!!! Failure convergence not reached" << std::endl;
            exit(127);
        } else {
            //this->m_log->o << "convergence reached " << std::endl<< std::flush;
        }
            //this->m_log->o << "output" << std::endl;

        *sp_u_approx_tstart = sp_u_tstop_approx;
//this->m_log->o << "end" << std::endl;
        StopLevelOperationTimer(LevelObserver::TL_STEP, l);
        //this->m_log->o << "debug::BraidIntegrator::Step[[end]]" << std::endl<<std::flush;
        return 0;
    };



    braid_Int Residual(braid_Vector u_,braid_Vector r_,BraidStepStatus &pstatus)override {
        //this->m_log->o << "debug::BraidIntegrator::Residual[[args]]" << std::endl<<std::flush;
        //print_status(this->m_log->o,pstatus);
        return 0;
        //this->m_log->o << "ERROR: called residual method for non residual supported integrator";
        //this->m_log->o << "debug::BraidIntegrator::Residual[[end]]" << std::endl<<std::flush;
        //return 0;
    };
    /* -------------------------------------------------------------------------------------
     * Methods
     * ---------------------------------------------------------------------------------- */
    void release() {
        //this->m_log->o << "debug::BraidIntegrator::release()" << std::endl<<std::flush;
#if TRACE_TIMINGS
        for (int i = Observer::T_INIT; i != Observer::N_OBSERVER; i++) {
            BraidUsageTimer &tel = this->m_time_log_manager.get(static_cast<Observer>(i));
            this->m_log->o << std::setw(20) << ObserverNames[i]
                              << std::setw(5) << " ;"
                              << std::setw(12) << tel.getTime() << ";"
                              << std::setw(12) << tel.getUsage() << ";"
                              << std::setw(12) << (tel.getTime() / tel.getUsage()) << ";"
                              << std::endl;

        }

        for (int i = 0; i != cLevelObserver; i++) {
            for (int l = 0; l < this->m_levels; l++) {
                BraidUsageTimer &tel = this->m_time_log_manager.get(static_cast<LevelObserver >(i), l);
                this->m_log->o << std::setw(20) << LevelObserverNames[i]
                                  << std::setw(5) << l << ";"
                                  << std::setw(12) << tel.getTime() << ";"
                                  << std::setw(12) << tel.getUsage() << ";"
                                  << std::setw(12) << (tel.getTime() / tel.getUsage()) << ";"
                                  << std::endl;
            }
        }
#endif
    }
    void print_settings() {
        this->m_log->o << "debug::BraidIntegrator::print_settings()" << std::endl<<std::flush;
        this->m_log->o << "===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====" << std::endl;
        this->m_log->o << "Application: " << this->m_name << std::endl;
        this->m_log->o << "Verbose: " << this->m_verbose << std::endl;
        this->m_log->o << "Max number of level: " << this->m_levels << std::endl;
        //this->m_log->o << "Solver: " << this->m_linSolver->name() << std::endl;
        //this->m_log->o << "Force convergence: " << this->m_force_conv << std::endl;
        //this->m_log->o << "Use adaptive ConvCheck: " << this->m_adapt_conv_check << std::endl;
        this->m_log->o << "===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====" << std::endl;
    }
    void setAdaptConv(bool b_state) { // todo setConvControl
        //this->m_adapt_conv_check = b_state;
    }

/*


    void setMaxLevels(size_t levelcount) {
#if TRACE_TIMINGS == 1
        if(levelcount >= 10 ) {
            std::cerr << "Warning: Tracing times is currently limited to 15 level" << std::endl;
        }
#endif
        this->m_levels = levelcount;
        this->m_levelStep.reserve(this->m_levels);
    }

    void add_stage(SPTimeIntegrator integrator) {
        m_level_integrator.template emplace_back(integrator);
    }
*/
    void set_time_integrator(SPTimeIntegrator integrator){
        m_time_integrator = integrator;
    }

    SPTimeIntegrator get_time_integrator(){
        return m_time_integrator;
    }

};
#endif //UG_PLUGIN_XBRAIDFORUG4_BRAIDINTEGRATOR_H