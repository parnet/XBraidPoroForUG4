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
        //print_status(this->m_log->o,pstatus);
        int level;
        int idone;
        int tindex;
        int iteration;

        double t_start;
        double t_stop;

        pstatus.GetLevel(&level);
        pstatus.GetTstartTstop(&t_start, &t_stop);
        pstatus.GetDone(&idone);
        pstatus.GetTIndex(&tindex);
        pstatus.GetIter(&iteration);

        double current_dt = t_stop - t_start;


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


        auto *sp_u_approx_tstart = (SPGridFunction *) u_->value;
        SPGridFunction sp_rhs = this->m_u0->clone_without_values();

        //SPGridFunction sp_u_tstop_approx = constsp_u_approx_tstop->get()->clone();
        //SPGridFunction lp = constsp_u_approx_tstop->get()->clone();

        //StartLevelOperationTimer(LevelObserver::TL_SOLVE,l);
        bool success;
        this->m_default_time_step = make_sp(new ug::ThetaTimeStep<TAlgebra>(this->m_domain_disc));
        this->m_default_time_step->set_theta(1.0); // implicit euler;

        auto solTimeSeries = make_sp(new ug::VectorTimeSeries<typename TAlgebra::vector_type>());
        solTimeSeries->push(sp_u_approx_tstart->get()->clone(), t_start);
        this->m_default_time_step->prepare_step(solTimeSeries, current_dt);

        const ug::GridLevel gridlevel = sp_u_approx_tstart->get()->grid_level();
        auto Operator_A = make_sp(new ug::AssembledLinearOperator<TAlgebra>(this->m_default_time_step, gridlevel));
        auto *ptr_Operator_A = Operator_A.get();
        this->m_default_time_step->assemble_jacobian(*ptr_Operator_A, *sp_u_approx_tstart->get(), gridlevel);


        this->m_default_time_step->assemble_rhs(*sp_rhs.get(), gridlevel);

        if (fstop_ != nullptr) { // apply correction for rhs
            auto fstop = *(SPGridFunction *) fstop_->value;
            VecAdd(1.0, *sp_rhs.get(), 1.0, *fstop.get());
        }

        linSolver->init(Operator_A, *sp_u_approx_tstart->get());
        success = linSolver->apply(*sp_u_approx_tstart->get(), *sp_rhs.get());

        if (!success) {
            this->m_log->o << "!!! Failure convergence not reached" << std::endl;
            exit(127);
        }


        //StopLevelOperationTimer(LevelObserver::TL_STEP,l);
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


        int level;
        int tindex;
        int iteration;

        double t_start;
        double t_stop;

        pstatus.GetLevel(&level);
        pstatus.GetTstartTstop(&t_start, &t_stop);

        pstatus.GetTIndex(&tindex);

        double current_dt = t_stop - t_start;

        this->m_script_log->o << "u_" << r_->index << " =  residual_"<<level<<"( u_" << u_->index << "," << t_stop
                              << " , u_" << r_->index          << ", "  << t_start << ")"
                              << " % " << tindex << std::endl;


        auto *const_u_approx_tstop = (SPGridFunction *) u_->value;
        auto *u_approx_tstart = (SPGridFunction *) r_->value;

        const ug::GridLevel gridlevel = const_u_approx_tstop->get()->grid_level();

        auto solTimeSeries = make_sp(new ug::VectorTimeSeries<typename TAlgebra::vector_type>());
        solTimeSeries->push(*u_approx_tstart, t_start);

        this->m_default_time_step = make_sp(new ug::ThetaTimeStep<TAlgebra>(this->m_domain_disc));
        this->m_default_time_step->set_theta(1.0); // implicit euler;
        this->m_default_time_step->prepare_step(solTimeSeries, current_dt);

        auto sp_rhs = this->m_u0->clone_without_values();
        auto Operator_A = make_sp(new ug::AssembledLinearOperator<TAlgebra>(this->m_default_time_step, gridlevel));
        auto *ptr_Operator_A = Operator_A.get();

        this->m_default_time_step->assemble_linear(*ptr_Operator_A, *sp_rhs.get(), gridlevel);
        this->linSolver->init(Operator_A, *const_u_approx_tstop->get());

        //this->m_default_time_step->assemble_rhs(*sp_rhs.get(), gridlevel);

        Operator_A->apply_sub(
                *sp_rhs.get(), // f co domain function [in / out]
                *const_u_approx_tstop->get() // u domain function [in]
        );

        (*sp_rhs) *= -1;
        *u_approx_tstart = sp_rhs;
        /*this->m_log->o << "debug::BraidTimeStepper::Residual[[args]]" << std::endl<<std::flush;
        //this->m_log->o << "debug::BraidTimeStepper::convert_inputs" << std::endl<<std::flush;
        auto u_tstop = *(SPGridFunction *) u_->value;
        auto u_tstart = *(SPGridFunction *) r_->value;

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


        this->m_default_time_step = make_sp(new ug::ThetaTimeStep<TAlgebra>(this->m_domain_disc));
        this->m_default_time_step->set_theta(1.0); // implicit euler;
        auto solTimeSeries = make_sp(new ug::VectorTimeSeries<typename TAlgebra::vector_type>());

        solTimeSeries->push(u_tstart->clone(), t_start);
        solTimeSeries->push(u_tstop->clone(), t_stop);
        const ug::GridLevel gridlevel = u_tstart->grid_level();
        this->m_default_time_step->prepare_step(solTimeSeries, current_dt);

        auto Operator_A = make_sp(new ug::AssembledLinearOperator<TAlgebra>(this->m_default_time_step, gridlevel));
        auto *ptr_Operator_A = Operator_A.get();
        this->m_default_time_step->assemble_jacobian(*ptr_Operator_A, *u_tstop.get(), gridlevel);

        //* \param[out] d 	Defect d(u) to be filled
        //* \param[in] 	u 	Current iterate
        //* \param[in]	gl	Grid Level
        //this->m_log->o << "debug::BraidTimeStepper::assemble_defect" << std::endl<<std::flush;
        this->m_default_time_step->assemble_defect(*u_tstart.get(),
                                                   *u_tstop.cast_const().get(),
                                                   gridlevel);

        // multiply with -1.0
        //this->m_log->o << "debug::BraidTimeStepper::invert" << std::endl<<std::flush;
        *(u_tstart).get() *= -1.0;

        //this->m_log->o << "debug::BraidTimeStepper::print" << std::endl<<std::flush;
        this->m_script_log->o << "u_" << r_->index << " =  residual_"<<level<<"( u_" << u_->index << "," << t_stop
                                << " , u_" << r_->index          << ", "  << t_start << ")"
                              << " % " << tindex << std::endl;



        // StopLevelOperationTimer(LevelObserver::TL_RESIDUAL, timegrid_level);
        this->m_log->o << "debug::BraidTimeStepper::Residual[[end]]" << std::endl<<std::flush;*/
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
