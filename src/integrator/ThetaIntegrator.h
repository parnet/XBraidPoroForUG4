//
// Created by parnet on 28.06.21.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_THETAINTEGRATOR_H
#define UG_PLUGIN_XBRAIDFORUG4_THETAINTEGRATOR_H

#include "lib_disc/operator/linear_operator/assembled_linear_operator.h"
#include "lib_disc/time_disc/theta_time_step.h"
#include "lib_algebra/operator/linear_solver/linear_solver.h"

template<typename TDomain, typename TAlgebra>
class ThetaIntegrator : public ug::ITimeIntegrator<TDomain, TAlgebra>  {
public:
    typedef ug::GridFunction<TDomain, TAlgebra> grid_function_type;

    typedef ug::IDomainDiscretization <TAlgebra> TDomainDisc;
    typedef SmartPtr <TDomainDisc> SPDomainDisc;

    typedef ug::LinearSolver<typename TAlgebra::vector_type> TSolver;
    typedef SmartPtr<TSolver> SPSolver;

    SPDomainDisc m_domain_disc;
    SPSolver m_linear_solver;

    ThetaIntegrator() : ug::ITimeIntegrator<TDomain,TAlgebra>() {
        //this->name = "Theta Integrator";
    }

    ~ThetaIntegrator() = default;

    void init(grid_function_type const &u) override {

    };

    void prepare(grid_function_type &u) override {

    };

    bool apply(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0) override{

        double current_dt = t1 - t0;

        auto thetatimestep = make_sp(new ug::ThetaTimeStep<TAlgebra>(m_domain_disc));
        thetatimestep->set_theta(1.0);

        auto solTimeSeries = make_sp(new ug::VectorTimeSeries<typename TAlgebra::vector_type>());

        // \todo: avoid this hack, use smart ptr properly haha
        SmartPtr<grid_function_type> u0_nonconst = u0.cast_const();
        solTimeSeries->push(u0_nonconst, t0);

        auto gridlevel = u0_nonconst->grid_level();
        thetatimestep->prepare_step(solTimeSeries, current_dt);

        auto Operator_A = make_sp(new ug::AssembledLinearOperator<TAlgebra>(thetatimestep, gridlevel));

        auto rhs = u0_nonconst->clone();

        thetatimestep->assemble_jacobian(*Operator_A.get(), *u0_nonconst.get(), gridlevel);
        thetatimestep->assemble_rhs(*rhs.get(), gridlevel);

        m_linear_solver->init(Operator_A, *u1.get());
        bool success = m_linear_solver->apply(*u1.get(), *rhs.get());
        return success;
    };

    void set_time_step(double dt) { this->m_dt = dt; }

    void set_domain(SPDomainDisc domain){
        this->m_domain_disc = domain;
    }

    void set_solver(SPSolver solver){
        this->m_linear_solver = solver;
    }
};


#endif //UG_PLUGIN_XBRAIDFORUG4_THETAINTEGRATOR_H
