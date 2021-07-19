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
class BraidTimeStepper : public BraidGridFunctionBase<TDomain, TAlgebra> {
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

    typedef ug::IDomainDiscretization<TAlgebra> TDomainDisc;
    typedef SmartPtr<TDomainDisc> SPDomainDisc;

    typedef ug::SimpleTimeIntegrator<TDomain, TAlgebra> TTimeIntegratorType;
    typedef ug::LinearImplicitEuler<TAlgebra> TTimeStepType;
    typedef SmartPtr<TTimeStepType> SPTimeStepType;
    typedef ug::ISubDiagErrorEst<TVector> TErrorEstimatorType;
    typedef ug::AitkenNevilleTimex<TVector> TTimexType;


    typedef ug::StdConvCheck<typename TAlgebra::vector_type> TConv;
    typedef SmartPtr<TConv> SPConv;


    typedef SmartPtr<SpaceTimeCommunicator> SPCommunicator;

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

public:

    /**
    * Note that this default constructor does not create a consistent object. The parameter t_comm (of type MPI_Comm)
    * for the temporal communication has to be set.
    */
    BraidTimeStepper() : BraidGridFunctionBase<TDomain, TAlgebra>() {
        this->m_name = "Braid Time Stepper";
        this->provide_residual = true;
        std::cout << "Braid Time Stepper :: default constructor was used!" << std::endl;
    }

    BraidTimeStepper(MPI_Comm mpi_temporal, double tstart, double tstop, int steps)
            : BraidGridFunctionBase<TDomain, TAlgebra>(mpi_temporal, tstart, tstop, steps) {
        this->provide_residual = true;
        this->m_name = "Braid Time Stepper";
    }

    ~BraidTimeStepper() = default;

    /* -------------------------------------------------------------------------------------
     * XBraid Method definitions
     * ---------------------------------------------------------------------------------- */

    braid_Int Step(braid_Vector u_, braid_Vector ustop_, braid_Vector fstop_, BraidStepStatus &pstatus) override {
        this->m_log->o << "debug::BraidTimeStepper::Step[[args]]" << std::endl<<std::flush;
        //print_status(this->m_log->o,pstatus);
        this->m_log->o << "debug::BraidTimeStepper::retrieve_variables" << std::endl<<std::flush;
        int level; // level
        pstatus.GetLevel(&level);
        // StartLevelOperationTimer(LevelObserver::TL_STEP, level);
        double t_start, t_stop;
        pstatus.GetTstartTstop(&t_start, &t_stop);
        int idone;
        pstatus.GetDone(&idone);
        double current_dt = t_stop - t_start;

        int tindex;
        pstatus.GetTIndex(&tindex);
        int iteration;
        pstatus.GetIter(&iteration);

        this->m_log->o << "debug::BraidTimeStepper::print_script" << std::endl<<std::flush;
        if (fstop_ == nullptr) {
            this->m_script_log->o << "u_" << u_->index << " = integrator_"<< level <<".apply(u_" << ustop_->index << ", "<< t_stop << ", u_" << u_->index << ", " << t_start << ")"
                                  << "\t\t % " << "level=" << level << " dt=" << current_dt << " t_index="<< tindex<<" iteration=" << iteration << " integrator: ThetaTimeStep"  << std::endl<< std::flush;
        } else {
            this->m_script_log->o << "u_" << u_->index << " = resintegrator_" << level
                                                << ".apply( u_" << ustop_->index <<", "<< t_stop
                                                << ", u_" << u_->index <<", " << t_start
                                                << ", r_" << fstop_->index << ")"
                    << "\t\t % " << "level=" << level << " dt=" << current_dt << " t_index="<< tindex<<" iteration=" << iteration << " integrator: ThetaTimeStep"  << std::endl<< std::flush;
        }

        //this->m_log->o << "debug::BraidTimeStepper::convert_inputs" << std::endl<<std::flush;
        auto *sp_u_approx_tstart = (SPGridFunction *) u_->value;
        auto *constsp_u_approx_tstop = (SPGridFunction *) ustop_->value;
        SPGridFunction sp_rhs = this->m_u0->clone_without_values(); // for rhs

        //SPGridFunction sp_u_tstop_approx = constsp_u_approx_tstop->get()->clone();
        //SPGridFunction lp = constsp_u_approx_tstop->get()->clone();

        //StartLevelOperationTimer(LevelObserver::TL_SOLVE,l);
        bool success;
        //this->m_log->o << "debug::BraidTimeStepper::create_timestep" << std::endl<<std::flush;
        this->m_default_time_step = make_sp(new ug::ThetaTimeStep<TAlgebra>(this->m_domain_disc));
        this->m_default_time_step->set_theta(1.0); // implicit euler;
        //this->m_log->o << "debug::BraidTimeStepper::create_time_series" << std::endl<<std::flush;
        auto solTimeSeries = make_sp(new ug::VectorTimeSeries<typename TAlgebra::vector_type>());
        solTimeSeries->push(sp_u_approx_tstart->get()->clone(), t_start);
        //this->m_log->o << "debug::BraidTimeStepper::create_operator" << std::endl<<std::flush;
        const ug::GridLevel gridlevel = sp_u_approx_tstart->get()->grid_level();
        auto Operator_A = make_sp(new ug::AssembledLinearOperator<TAlgebra>(this->m_default_time_step, gridlevel));

        //this->m_log->o << "debug::BraidTimeStepper::prepare_step" << std::endl<<std::flush;
        auto rhs = sp_u_approx_tstart->get()->clone();
        this->m_default_time_step->prepare_step(solTimeSeries, current_dt);

        //this->m_log->o << "debug::BraidTimeStepper::assemble_operator" << std::endl<<std::flush;
        auto *ptr_Operator_A = Operator_A.get();
        this->m_default_time_step->assemble_jacobian(*ptr_Operator_A, *sp_u_approx_tstart->get(), gridlevel);
        //this->m_log->o << "debug::BraidTimeStepper::assemble_rhs" << std::endl<<std::flush;
        this->m_default_time_step->assemble_rhs(*rhs.get(), gridlevel);

        //this->m_log->o << "debug::BraidTimeStepper::adapt_rhs" << std::endl<<std::flush;
        if (fstop_ != nullptr) {
            //this->m_log->o << "residual used" << std::endl << std::flush;
            auto fstop = *(SPGridFunction *) fstop_->value;
            VecAdd(1.0, *rhs.get(), 1.0, *fstop.get());
        } else {
            this->m_log->o << "residual not used " << std::endl << std::flush;
        }
        //this->m_log->o << "debug::BraidTimeStepper::init_linear_solver" << std::endl<<std::flush;
        linSolver->init(Operator_A, *sp_u_approx_tstart->get());
        //this->m_log->o << "debug::BraidTimeStepper::apply_linear_solver" << std::endl<<std::flush;
        success = linSolver->apply(*sp_u_approx_tstart->get(), *rhs.get());
        //StopLevelOperationTimer(LevelObserver::TL_SOLVE,l);
        //this->m_log->o << "debug::BraidTimeStepper::solver_finished" << std::endl<<std::flush;

        if (!success) {
            this->m_log->o << "!!! Failure convergence not reached" << std::endl;
            exit(127);
        } else {
            this->m_log->o <<"converged" << std::endl;
        }


        //StopLevelOperationTimer(LevelObserver::TL_STEP,l);
        this->m_log->o << "debug::BraidTimeStepper::Step[[end]]" << std::endl<<std::flush;
        return 0;
    };




    /** @brief Compute the residual at time @a tstop, given the approximate
      solutions at @a tstart and @a tstop. The values of @a tstart and @a tstop
      can be obtained from @a pstatus.

      @param[in]     u_ Input: approximate solution at time @a tstop.
      @param[in,out] r_ Input: approximate solution at time @a tstart.
                        Output: residual at time @a tstop.

      @see braid_PtFcnResidual.
  */
    braid_Int Residual(braid_Vector u_, braid_Vector r_, BraidStepStatus &pstatus) override {
        this->m_log->o << "debug::BraidTimeStepper::Residual[[args]]" << std::endl<<std::flush;
        //this->m_log->o << "debug::BraidTimeStepper::convert_inputs" << std::endl<<std::flush;
        auto sp_u_approx_tstop = *(SPGridFunction *) u_->value;
        auto sp_u_approx_tstart = *(SPGridFunction *) r_->value;

        //this->m_log->o << "debug::BraidTimeStepper::retrieve_variables" << std::endl<<std::flush;
        int level; // level;
        pstatus.GetLevel(&level);
        // StartLevelOperationTimer(LevelObserver::TL_RESIDUAL, level);
        double t_start, t_stop;
        pstatus.GetTstartTstop(&t_start, &t_stop);
        double current_dt = t_stop - t_start;

        int tindex;
        pstatus.GetTIndex(&tindex);
        int iteration;
        pstatus.GetIter(&iteration);

        //this->m_log->o << "debug::BraidTimeStepper::create_time_step" << std::endl<<std::flush;
        this->m_default_time_step = make_sp(new ug::ThetaTimeStep<TAlgebra>(this->m_domain_disc));
        this->m_default_time_step->set_theta(1.0); // implicit euler;
        //this->m_log->o << "debug::BraidTimeStepper::create_time_series" << std::endl<<std::flush;
        auto solTimeSeries = make_sp(new ug::VectorTimeSeries<typename TAlgebra::vector_type>());
        solTimeSeries->push(sp_u_approx_tstart->clone(), t_start);
        //this->m_log->o << "debug::BraidTimeStepper::prepare_step" << std::endl<<std::flush;
        const ug::GridLevel gridlevel = sp_u_approx_tstart->grid_level();
        this->m_default_time_step->prepare_step(solTimeSeries, current_dt);

//    StartLevelOperationTimer(LevelObserver::TL_ASSEMBLE_OP,timegrid_level);
//    StopLevelOperationTimer(LevelObserver::TL_ASSEMBLE_OP,timegrid_level);
//    StartLevelOperationTimer(LevelObserver::TL_ASSEMBLE_RHS,timegrid_level);
//    StopLevelOperationTimer(LevelObserver::TL_ASSEMBLE_RHS,timegrid_level);
        //* \param[out] d 	Defect d(u) to be filled
        //* \param[in] 	u 	Current iterate
        //* \param[in]	gl	Grid Level
        //this->m_log->o << "debug::BraidTimeStepper::assemble_defect" << std::endl<<std::flush;
        this->m_default_time_step->assemble_defect(*sp_u_approx_tstart.get(),
                                                   *sp_u_approx_tstop.cast_const().get(),
                                                   gridlevel);

        // multiply with -1.0
        //this->m_log->o << "debug::BraidTimeStepper::invert" << std::endl<<std::flush;
        *(sp_u_approx_tstart).get() *= -1.0;

        //this->m_log->o << "debug::BraidTimeStepper::print" << std::endl<<std::flush;
        this->m_script_log->o << "u_" << r_->index << " =  residual_"<<level<<"( u_" << u_->index << " , u_" << r_->index
                              << ", "
                              << t_start << "," << t_stop << ")"
                              << " % " << tindex << std::endl;



        // StopLevelOperationTimer(LevelObserver::TL_RESIDUAL, timegrid_level);
        this->m_log->o << "debug::BraidTimeStepper::Residual[[end]]" << std::endl<<std::flush;
        return 0;
    }

    void print_settings() {
    }

    void setAdaptConv(bool conv) {
    }

    void setForceConv(bool force) {

    }

    void set_solver(SPSolver solver) {
        this->linSolver = solver;
    }
};

#endif //UG_PLUGIN_XBRAIDFORUG4_BRAIDTIMESTEPPER_H
