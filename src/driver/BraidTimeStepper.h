// zzzzz
// Created by parnet on 28.05.21.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_BRAIDTIMESTEPPER_H
#define UG_PLUGIN_XBRAIDFORUG4_BRAIDTIMESTEPPER_H

#include "BraidGridFunctionBase.h"
#include "lib_disc/operator/linear_operator/assembled_linear_operator.h"
#include "lib_disc/time_disc/theta_time_step.h"
#include "lib_algebra/operator/linear_solver/linear_solver.h"

template<typename TDomain, typename TAlgebra>
class BraidTimeStepper : public BraidGridFunctionBase<TDomain,TAlgebra> {
public:
    /* *************************************************************************************
     * Type definitions
     ************************************************************************************* */
    typedef typename TAlgebra::matrix_type TMatrix;
    typedef typename TAlgebra::vector_type TVector;

    typedef ug::GridFunction<TDomain, TAlgebra> TGridFunction;
    typedef SmartPtr<TGridFunction> SPGridFunction;


    typedef ug::ThetaTimeStep<TAlgebra> TTimeStep;
    typedef SmartPtr<TTimeStep> SPTimeStep;

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

    typedef ug::LinearSolver<typename TAlgebra::vector_type> TSolver;
    typedef SmartPtr<TSolver> SPSolver;
private:
    /* -------------------------------------------------------------------------------------
     * Member variables
     * ---------------------------------------------------------------------------------- */

private:
    SPTimeStep m_default_time_step;
    SPTimeStep m_fine_time_step;
    SPTimeStep m_coarse_time_step;
    SPSolver linSolver;
    SPDomainDisc domainDisc;

public:

    /**
    * Note that this default constructor does not create a consistent object. The parameter t_comm (of type MPI_Comm)
    * for the temporal communication has to be set.
    */
    BraidTimeStepper() : BraidGridFunctionBase<TDomain,TAlgebra>() {
        this->m_name = "Braid Time Stepper";
        this-> provide_residual = true;
        std::cout << "Braid Time Stepper :: default constructor was used!" << std::endl;
    }
    BraidTimeStepper(MPI_Comm mpi_temporal, double tstart, double tstop, int steps) : BraidGridFunctionBase<TDomain,TAlgebra>(mpi_temporal, tstart,tstop, steps) {
        this-> provide_residual = true;
        this->m_name = "Braid Time Stepper";
    }
    ~BraidTimeStepper() = default;

    /* -------------------------------------------------------------------------------------
     * XBraid Method definitions
     * ---------------------------------------------------------------------------------- */

    braid_Int Step(braid_Vector u_,braid_Vector ustop_,braid_Vector fstop_,BraidStepStatus &pstatus) override {
        this->m_log->o << "debug::BraidTimeStepper::Step[[args]]" << std::endl<<std::flush;
        print_status(this->m_log->o,pstatus);
#if TRACE_CONST == 1
        u->m_const = false;
//ustop->m_const = false;
#endif
        int l; // level
        pstatus.GetLevel(&l);
        StartLevelOperationTimer(LevelObserver::TL_STEP,l);
        double t_start, t_stop;
        pstatus.GetTstartTstop(&t_start, &t_stop);

        int idone;
        pstatus.GetDone(&idone);

        double current_dt = t_stop - t_start;

// std::cout << "num stages " << this->m_timeDisc->num_stages() << std::endl;

//this->o << "adaptConv" << std::endl;

        //if (this->adaptConv) {  } // todo
//this->m_log->o << "message" << std::endl;
        if (this->m_verbose) {
            int tindex;
            pstatus.GetTIndex(&tindex);
            int iteration;
            pstatus.GetIter(&iteration);
#if TRACE_INDEX == 1
            /*if (fstop_ == nullptr) {
                this->debugwriter << "u_" << u_->index << " = step_"<<l<<"_n( u_" << u_->index << ", u_" << ustop_->index
                                  << ", null, " << t_start << ", " << current_dt << ", " << t_stop << ", " << l
                                  << ")"
                                  << "\t\t % " << tindex << std::endl;
            } else {
                this->debugwriter << "u_" << u_->index << " = step_"<<l<<"_r( u_" << u_->index << ", u_" << ustop_->index
                                  << ", u_"
                                  << fstop_->index << ", " << t_start << ", " << current_dt << ", " << t_stop << ", "
                                  << l << ")"
                                  << " % " << tindex << std::endl;
            }*/

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


        // todo adapt conv check?
        if (fstop_ != nullptr) {
            this->m_log->o << "Warning residul is ignored" << std::endl<< std::flush;
        } else {
            this->m_log->o << "residual not used " << std::endl<< std::flush;
        }


//this->m_log->o << "fstop" << std::endl;


//this->m_log->o << "solve" << std::endl;
        //StartLevelOperationTimer(LevelObserver::TL_SOLVE,l);

        bool success;

        //this->m_log->o << "165"<<std::flush << std::endl;
        this->m_default_time_step = make_sp(new ug::ThetaTimeStep<TAlgebra>(domainDisc));

        //this->m_log->o << "168"<<std::flush << std::endl;
        this->m_default_time_step->set_theta(1.0); // implicit euler;

        //this->m_log->o << "171"<<std::flush << std::endl;
        auto solTimeSeries = make_sp(new ug::VectorTimeSeries<typename TAlgebra::vector_type>());

        //this->m_log->o << "174"<<std::flush << std::endl;
        solTimeSeries->push(sp_u_approx_tstart->get()->clone(), t_start);

        //this->m_log->o << "177"<<std::flush << std::endl;
        const ug::GridLevel gridlevel = sp_u_approx_tstart->get()->grid_level();

        //this->m_log->o << "180"<<std::flush << std::endl;
        auto Operator_A = make_sp(new ug::AssembledLinearOperator<TAlgebra>(this->m_default_time_step, gridlevel));

        //this->m_log->o << "183"<<std::flush << std::endl;
        auto rhs = sp_u_approx_tstart->get()->clone();

        //this->m_log->o << "186"<<std::flush << std::endl;
        this->m_default_time_step->prepare_step(solTimeSeries,current_dt);

        //this->m_log->o << "189"<<std::flush << std::endl;
        auto * ptr_Operator_A = Operator_A.get();

        //this->m_log->o << "191"<<std::flush << std::endl;
        this->m_default_time_step->assemble_jacobian(*ptr_Operator_A, *sp_u_approx_tstart->get(), gridlevel);

        //this->m_log->o << "195"<<std::flush << std::endl;
        this->m_default_time_step->assemble_rhs(*rhs.get(), gridlevel);

        //this->m_log->o << "198"<<std::flush << std::endl;
        linSolver->init(Operator_A, *sp_u_approx_tstart->get());
        //linSolver->prepare(*sp_u_approx_tstart->get());

        //this->m_log->o << "202"<<std::flush << std::endl;
        success = linSolver->apply(*sp_u_approx_tstart->get(),*rhs.get());

        //this->m_log->o << "205"<<std::flush << std::endl;

        //ptr_time_integrator->apply(u_end, t_end, u_start, t_start);
        //StopLevelOperationTimer(LevelObserver::TL_SOLVE,l);

//#if TRACE_DEFECT == 1
//        if(this->m_verbose){
            //this->m_log->o << std::setw(20) << "@conv Iterations: " << std::setw(12) << m_linSolver->step()
            //                  << std::setw(20) << "Reduction: " << std::setw(12) << m_linSolver->reduction()
            //                  << std::setw(20) << "Defect: " << std::setw(12) << m_linSolver->defect() << std::endl;
        //}
//#endif

        if (!success) {
            this->m_log->o << "!!! Failure convergence not reached" << std::endl;
            exit(127);
        } //else {
            //this->m_log->o << "converged " << std::flush << std::endl;
        //}

        //this->m_log->o << "output" << std::endl;
        //*sp_u_approx_tstart = sp_u_tstop_approx;
        //series->clear();
//this->m_log->o << "end" << std::endl;
        //StopLevelOperationTimer(LevelObserver::TL_STEP,l);
        return 0;
    };

    braid_Int Residual(braid_Vector u_,braid_Vector r_,BraidStepStatus &pstatus) override {
#if TRACE_CONST == 1
        r->m_const = false;
// u->m_const = false;
#endif
        int timegrid_level; // level;
        pstatus.GetLevel(&timegrid_level);
        StartLevelOperationTimer(LevelObserver::TL_RESIDUAL, timegrid_level);
        double t_start, t_stop;
        pstatus.GetTstartTstop(&t_start, &t_stop);
        double current_dt = t_stop - t_start;

        int tindex;
        pstatus.GetTIndex(&tindex);

#if TRACE_INDEX == 1
        /*if (this->m_verbose) {
            this->debugwriter << "u_" << r_->index << " =  residual( u_" << u_->index << " , u_" << r_->index
                              << ", "
                              << t_start << ", " << current_dt << ", " << t_stop << ", " << timegrid_level << ")"
                              << " % " << tindex << std::endl;
        }*/
#endif

        auto *const_u_approx_tstop = (SPGridFunction *) u_->value;
        auto *u_approx_tstart = (SPGridFunction *) r_->value;


        const ug::GridLevel grid_level = const_u_approx_tstop->get()->grid_level();

        // series->push(*u_approx_tstart, t_start);
        // m_timeDisc->prepare_step(series, current_dt);
        // auto sp_rhs = this->m_u0->clone_without_values();
        // if (fabs(current_dt - this->m_assembled_dt) > 1e-14) {
        //    StartLevelOperationTimer(LevelObserver::TL_ASSEMBLE_OP,timegrid_level);
        //    m_timeDisc->assemble_linear(*m_A, *sp_rhs.get(), gridlevel);

        //    m_linSolver->init(m_A, *const_u_approx_tstop->get());

        //    this->m_assembled_dt = current_dt;
        //    StopLevelOperationTimer(LevelObserver::TL_ASSEMBLE_OP,timegrid_level);
        //} else {
        //    StartLevelOperationTimer(LevelObserver::TL_ASSEMBLE_RHS,timegrid_level);
        //    m_timeDisc->assemble_rhs(*sp_rhs.get(), gridlevel);
        //    StopLevelOperationTimer(LevelObserver::TL_ASSEMBLE_RHS,timegrid_level);
        // }
        // m_linSolver->linear_operator()->apply_sub(
        //        *sp_rhs.get(), // f co domain function [in / out]
        //        *const_u_approx_tstop->get() // u domain function [in]
        // ); // calculates r = r - A * u

        // auto & levelStep = this->m_levelStep[timegrid_level];
        // levelStep.m_stepper->assemble_defect(*sp_rhs.get(), ,grid_level);
        // r = rhs  - L*u_stop
        // (*sp_rhs) *= -1; // r = -r + A*u
        // *u_approx_tstart = sp_rhs;
        // SPGridFunction  def = sp_rhs->clone(); // todo alternative?
        // m_timeDisc->assemble_defect(*def.get(),*const_u_approx_tstop->get(),gridlevel);
        // VecAdd(1,*def.get(),-1,*u_approx_tstart->get());

        // *u_approx_tstart = *const_u_approx_tstop;
        // series->clear();
        // StopLevelOperationTimer(LevelObserver::TL_RESIDUAL, timegrid_level);
        return 0;
    }

    void print_settings(){
    }
    void setAdaptConv(bool conv){
    }
    void setForceConv(bool force){

    }

    void set_domain(SPDomainDisc domainDisc) {
        this->domainDisc = domainDisc;
    }

    void set_solver(SPSolver solver) {
        this->linSolver = solver;
    }
};

#endif //UG_PLUGIN_XBRAIDFORUG4_BRAIDTIMESTEPPER_H
