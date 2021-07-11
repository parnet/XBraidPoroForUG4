//
// Created by parnet on 07.06.21.
//

//
// Created by parnet on 24.05.21.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_BRAIDINTEGRATORFACTORY_H
#define UG_PLUGIN_XBRAIDFORUG4_BRAIDINTEGRATORFACTORY_H


#include "BraidGridFunctionBase.h"
#include "../interface/integrator_factory.h"

template<typename TDomain, typename TAlgebra>
class LevelIntegratorFactory {
public: // todo set better modes

    typedef ug::ITimeIntegrator<TDomain, TAlgebra> TTimeIntegrator;
    typedef SmartPtr<TTimeIntegrator> SPTimeIntegrator;

    SPTimeIntegrator m_time_integrator_factory;

    void set_time_integrator(SPTimeIntegrator integrator){
        m_time_integrator_factory = integrator;
    }

    SPTimeIntegrator get_time_integrator(){
        return m_time_integrator_factory;
    }
};




template<typename TDomain, typename TAlgebra>
class BraidIntegratorFactory : public BraidGridFunctionBase<TDomain,TAlgebra> {
public: // todo set better modes
    /* *************************************************************************************
     * Type definitions
     ************************************************************************************* */
    typedef typename TAlgebra::matrix_type TMatrix;
    typedef typename TAlgebra::vector_type TVector;
    typedef ug::GridFunction<TDomain, TAlgebra> TGridFunction;
    typedef SmartPtr<TGridFunction> SPGridFunction;


    typedef IntegratorFactory<TDomain,TAlgebra> TIntegratorFactory;
    typedef SmartPtr<TIntegratorFactory> SPIntegratorFactory;

    typedef ug::ITimeIntegrator<TDomain, TAlgebra> TTimeIntegrator;
    typedef SmartPtr<TTimeIntegrator> SPTimeIntegrator;
    typedef std::vector<SPTimeIntegrator> LevelIntegratorList;

    typedef ug::INonlinearTimeIntegrator<TDomain, TAlgebra> TBaseType; // todo derivate ?
    typedef ug::INonlinearTimeIntegrator<TDomain, TAlgebra> TINonlinearTimeIntegrator;
    typedef typename TBaseType::solver_type TSolverType;
    typedef SmartPtr<TSolverType> SPSolverType;

    typedef ug::IDomainDiscretization<TAlgebra>	TDomainDisc;
    typedef SmartPtr<TDomainDisc> SPDomainDisc;

    //typedef ug::SimpleTimeIntegrator<TDomain, TAlgebra> TTimeIntegratorType;
    //typedef ug::LinearImplicitEuler<TAlgebra> TTimeStepType;
    //typedef SmartPtr<TTimeStepType > SPTimeStepType;
    //typedef ug::ISubDiagErrorEst<TVector> TErrorEstimatorType;
    //typedef ug::AitkenNevilleTimex<TVector> TTimexType;



    typedef ug::StdConvCheck<typename TAlgebra::vector_type> TConv;
    typedef SmartPtr<TConv> SPConv;


    typedef SmartPtr<MATLABScriptor<TDomain, TAlgebra>> SPMATLABScriptor;
    typedef SmartPtr <SpaceTimeCommunicator> SPCommunicator;

    typedef Paralog TParalog;
    typedef SmartPtr<TParalog> SPParalog;


private:
    /* -------------------------------------------------------------------------------------
     * Member variables
     * ---------------------------------------------------------------------------------- */

private:
    SPTimeIntegrator  m_default_integrator;
    SPTimeIntegrator  m_fine_integrator;
    SPTimeIntegrator  m_coarse_integrator;
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


    SPIntegratorFactory m_fine_time_integrator_factory;
    SPIntegratorFactory m_coarse_time_integrator_factory;
    /**
    * Note that this default constructor does not create a consistent object. The parameter t_comm (of type MPI_Comm)
    * for the temporal communication has to be set.
    */
    BraidIntegratorFactory() : BraidGridFunctionBase<TDomain,TAlgebra>() {
        //this->m_log->o << "debug::BraidIntegrator()" << std::endl<<std::flush;
        this->m_name = "Braid Integrator Factory";
        //this->m_log->o << "Warning Braid Integrator :: default constructor was used!" << std::endl;
    }

    BraidIntegratorFactory(MPI_Comm mpi_temporal, double tstart, double tstop, int steps) : BraidGridFunctionBase<TDomain,TAlgebra>(mpi_temporal, tstart,tstop, steps) {
        //this->m_log->o << "debug::BraidIntegrator(args)" << std::endl<<std::flush;
        this->m_name = "Braid Integrator Factory";
    }

    ~BraidIntegratorFactory() = default;
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
        int level; // level
        pstatus.GetLevel(&level);
        std::cout <<"XBRAID LEVEL = " << level << std::endl;
        StartLevelOperationTimer(LevelObserver::TL_STEP, level);
        double t_start, t_stop;
        pstatus.GetTstartTstop(&t_start, &t_stop);
        int idone;
        pstatus.GetDone(&idone);

        double current_dt = t_stop - t_start;
        if (this->m_verbose) {
            int tindex;
            pstatus.GetTIndex(&tindex);
            int iteration;
            pstatus.GetIter(&iteration);
#if TRACE_INDEX == 1


            if (fstop_ == nullptr) {
                this->m_log->o << "u_" << u_->index << " = step_" << level << "_n( u_" << u_->index << ", u_" << ustop_->index
                               << ", null, " << t_start << ", " << current_dt << ", " << t_stop << ", " << level
                               << ")"
                               << "\t\t % " << "idx("<< tindex<<") iter=" << iteration << std::endl<< std::flush;
            } else {
                this->m_log->o << "u_" << u_->index << " = step_" << level << "_r( u_" << u_->index << ", u_" << ustop_->index
                               << ", u_"
                               << fstop_->index << ", " << t_start << ", " << current_dt << ", " << t_stop << ", "
                               << level << ")"
                               << " % " << tindex << std::endl << std::flush;
            }

#else
            //this->m_log->o << std::setw(13) << iteration << "step for level " << level << " at position " << tindex << " and iteration" << std::endl;
#endif
        }
        auto *sp_u_approx_tstart = (SPGridFunction *) u_->value;
        auto *constsp_u_approx_tstop = (SPGridFunction *) ustop_->value;

        SPGridFunction sp_u_tstop_approx = constsp_u_approx_tstop->get()->clone();
        SPGridFunction lp = constsp_u_approx_tstop->get()->clone();
        SPGridFunction sp_rhs = this->m_u0->clone_without_values(); // for rhs


        int tindex;
        pstatus.GetTIndex(&tindex);
        int iteration;
        pstatus.GetIter(&iteration);
        {
            SPGridFunction tempobject_b = sp_u_approx_tstart->get()->clone(); // clone to ensure consistency
            this->vtkScriptor->set_filename("sp_u_approx_tstart_before");
            this->vtkScriptor->write(tempobject_b,tindex,t_start,iteration,level);
        }

        {
            SPGridFunction tempobject_b = sp_u_tstop_approx->clone(); // clone to ensure consistency
            this->vtkScriptor->set_filename("sp_u_tstop_approx_before");
            this->vtkScriptor->write(tempobject_b,tindex,t_stop,iteration,level);
        }

        //const ug::GridLevel gridlevel = sp_u_approx_tstart->get()->grid_level();
        // todo adapt conv check?
        if (fstop_ != nullptr) {
            //this->m_log->o << "Warning residul is ignored" << std::endl<< std::flush;
        } else {
            //this->m_log->o << "residual not used " << std::endl<< std::flush;
        }


//this->m_log->o << "fstop" << std::endl;

#if TRACE_INDEX == 1
        MATLAB(sp_rhs->clone().get(), u->index, t_stop);
#endif
//this->m_log->o << "solve" << std::endl;
        StartLevelOperationTimer(LevelObserver::TL_SOLVE, level);
        // smartptr, const smartptr
        //this->m_log->o << "integrator creation and application starts " << std::endl<< std::flush;
        SPTimeIntegrator loc_time_integrator;
        std::cout << "XBRAID Integrator solving level = " << level << std::endl;
        if(level <= 0) {
            std::cout << level<< "XBRAID Integrator using fine Integrator "  << m_fine_time_integrator_factory->m_name << std::endl;
            loc_time_integrator = m_fine_time_integrator_factory->create_level_time_integrator(current_dt, bool(idone), level);
//            this->m_log->o << m_fine_time_integrator_factory->get_name() << std::endl;
        } else {
            std::cout << level<<"XBRAID Integrator using coarse Integrator "  << m_coarse_time_integrator_factory->m_name << std::endl;
            loc_time_integrator = m_coarse_time_integrator_factory->create_level_time_integrator(current_dt, bool(idone), level);
//            this->m_log->o << m_coarse_time_integrator_factory->get_name() << std::endl;
        }


        //SPGridFunction tempobject_a = sp_u_approx_tstart->get()->clone(); // clone to ensure consistency
        loc_time_integrator->init(*sp_u_approx_tstart->get()->clone());
        loc_time_integrator->prepare(*sp_u_approx_tstart->get()->clone());





        bool success = loc_time_integrator->apply(sp_u_tstop_approx, t_stop, sp_u_approx_tstart->cast_const(), t_start);
        //this->m_log->o << "integrator applied "<< std::endl << std::flush;
        StopLevelOperationTimer(LevelObserver::TL_SOLVE, level);
        if (!success ){//&& m_force_conv) {
            this->m_log->o << "!!! Failure convergence not reached" << std::endl;
            exit(127);
        } else {
            //this->m_log->o << "convergence reached " << std::endl<< std::flush;
        }
#if TRACE_INDEX == 1
        MATLAB(sp_u_tstop_approx.get(), u->index, t_stop);
#endif
        *sp_u_approx_tstart = sp_u_tstop_approx; // u_tstart is the return value


        {
            SPGridFunction tempobject_b = sp_u_approx_tstart->get()->clone(); // clone to ensure consistency
            this->vtkScriptor->set_filename("sp_u_approx_tstart_after");
            this->vtkScriptor->write(tempobject_b,tindex,t_stop,iteration,level);
        }

        {
            SPGridFunction tempobject_b = sp_u_tstop_approx->clone(); // clone to ensure consistency
            this->vtkScriptor->set_filename("sp_u_tstop_approx_after");
            this->vtkScriptor->write(tempobject_b,tindex,t_stop,iteration,level);
        }
        this->vtkScriptor->set_filename("access");

        StopLevelOperationTimer(LevelObserver::TL_STEP, level);
        //this->m_log->o << "debug::BraidIntegrator::Step[[end]]" << std::endl<<std::flush;
        return 0;
    };



    braid_Int Residual(braid_Vector u_,braid_Vector r_,BraidStepStatus &pstatus)override {
        //this->m_log->o << "debug::BraidIntegrator::Residual[[args]]" << std::endl<<std::flush;
        //print_status(this->m_log->o,pstatus);
        this->m_log->o << "residual" << std::flush << std::endl;
        //this->m_log->o << "ERROR: called residual method for non residual supported integrator";
        //this->m_log->o << "debug::BraidIntegrator::Residual[[end]]" << std::endl<<std::flush;
        return 0;
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

    void set_max_levels(int num_level) {
        BraidGridFunctionBase<TDomain,TAlgebra>::set_max_levels(num_level);
        this->m_level_integrator.resize(num_level,SPNULL);
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
    void set_fine_time_integrator(SPIntegratorFactory integrator){
        this->m_fine_time_integrator_factory = integrator;
    }

    void set_coarse_time_integrator(SPIntegratorFactory integrator){
        this->m_coarse_time_integrator_factory = integrator;
    }

    SPIntegratorFactory get_fine_time_integrator(){
        return m_fine_time_integrator_factory;
    }

    SPIntegratorFactory get_coarse_time_integrator(){
        return m_coarse_time_integrator_factory;
    }

};
#endif //UG_PLUGIN_XBRAIDFORUG4_BRAIDINTEGRATORFACTORY_H