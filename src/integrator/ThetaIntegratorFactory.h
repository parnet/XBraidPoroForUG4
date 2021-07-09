//
// Created by parnet on 28.06.21.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_THETAINTEGRATORFACTORY_H
#define UG_PLUGIN_XBRAIDFORUG4_THETAINTEGRATORFACTORY_H

#include "../interface/integrator_factory.h"

template<typename TDomain, typename TAlgebra>
class ThetaIntegratorFactory : public IntegratorFactory<TDomain,TAlgebra> {
public:

    typedef ug::ITimeIntegrator<TDomain, TAlgebra> TTimeIntegrator;
    typedef SmartPtr<TTimeIntegrator> SPTimeIntegrator;

    typedef ug::IDomainDiscretization <TAlgebra> TDomainDisc;
    typedef SmartPtr <TDomainDisc> SPDomainDisc;

    typedef ug::LinearSolver<typename TAlgebra::vector_type> TSolver;
    typedef SmartPtr<TSolver> SPSolver;

    typedef ThetaIntegrator<TDomain,TAlgebra> TThetaIntegrator;
    typedef SmartPtr<TThetaIntegrator> SPThetaIntegrator;

    SPSolver m_linear_solver;
    SPDomainDisc m_domain_disc;

    ThetaIntegratorFactory() : IntegratorFactory<TDomain, TAlgebra>(){
        this->m_name = "Theta Integrator";
    }

    ~ThetaIntegratorFactory() = default;
    SPTimeIntegrator create_level_time_integrator(double current_dt, bool done,int level)override {
        return create_time_integrator(current_dt, done);
    }

    SPTimeIntegrator create_time_integrator(double current_dt, bool done)override{
        SPThetaIntegrator integrator = make_sp(new ThetaIntegrator<TDomain,TAlgebra>());
        integrator->set_domain(m_domain_disc);
        integrator->set_solver(m_linear_solver);
        return integrator;
    };


    void set_domain(SPDomainDisc domain){
        this->m_domain_disc = domain;
    }

    void set_solver(SPSolver solver){
        this->m_linear_solver = solver;
    }
};

#endif //UG_PLUGIN_XBRAIDFORUG4_THETAINTEGRATORFACTORY_H
