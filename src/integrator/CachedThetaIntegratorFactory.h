//
// Created by parnet on 28.06.21.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_CACHEDTHETAINTEGRATORFACTORY_H
#define UG_PLUGIN_XBRAIDFORUG4_CACHEDTHETAINTEGRATORFACTORY_H

#include "../interface/integrator_factory.h"

template<typename TDomain, typename TAlgebra>
class CachedThetaIntegratorFactory : public IntegratorFactory<TDomain,TAlgebra> {
public:

    typedef ug::ITimeIntegrator<TDomain, TAlgebra> TTimeIntegrator;
    typedef SmartPtr<TTimeIntegrator> SPTimeIntegrator;

    typedef ug::IDomainDiscretization <TAlgebra> TDomainDisc;
    typedef SmartPtr <TDomainDisc> SPDomainDisc;

    typedef ug::LinearSolver<typename TAlgebra::vector_type> TSolver;
    typedef SmartPtr<TSolver> SPSolver;

    typedef CachedThetaIntegrator<TDomain,TAlgebra> TThetaIntegrator;
    typedef SmartPtr<TThetaIntegrator> SPThetaIntegrator;

    SPSolver m_linear_solver;
    SPDomainDisc m_domain_disc;

    double last_dt = 1e10;
    int last_idx = 0;

    std::vector<SPTimeIntegrator> integrators;
    std::vector<double> dts;

    CachedThetaIntegratorFactory() : IntegratorFactory<TDomain, TAlgebra>(){
        this->m_name = "Theta Integrator (cached)";
    }

    ~CachedThetaIntegratorFactory() = default;

    SPTimeIntegrator create_level_time_integrator(double current_dt, bool done,int level) override{
        return create_time_integrator(current_dt,done);
    }

    SPTimeIntegrator create_time_integrator(double current_dt, bool done) override{
        if(current_dt < last_dt){
            //std::cout << "Talasma :: Stepsize decreased create integrator "<< std::endl;
            SPThetaIntegrator integrator = make_sp(new CachedThetaIntegrator<TDomain,TAlgebra>());
            integrator->set_domain(m_domain_disc);
            integrator->set_solver(m_linear_solver);
            //std::cout << "Talamsa :: Managing " << std::endl;
            integrators.push_back(integrator);
            dts.push_back(current_dt);
            // add to list
            last_dt =  current_dt;
            //std::cout << "Talamsa :: return " << std::endl;
            return integrator;
        } else {
            for(int i = 0; i < dts.size(); i++){
                if (fabs(dts[i] - current_dt < 1e-15)){
                    //std::cout << "Talasma :: Found integrator -> reuse it"<< std::endl;
                    return integrators[i];
                }
            }
            //std::cout << "Talasma :: Fallback create integrator -> reuse it"<< std::endl;

            // just for safety.
            SPThetaIntegrator integrator = make_sp(new CachedThetaIntegrator<TDomain,TAlgebra>());
            integrator->set_domain(m_domain_disc);
            integrator->set_solver(m_linear_solver);
            //std::cout << "Talamsa :: Managing :: prime" << std::endl;
            integrators.template emplace_back(integrator);
            dts.template emplace_back(current_dt);
            //std::cout << "Talamsa :: return prime" << std::endl;
            // add to list
            return integrator;

        }
    };


    void set_domain(SPDomainDisc domain){
        this->m_domain_disc = domain;
    }

    void set_solver(SPSolver solver){
        this->m_linear_solver = solver;
    }
};

#endif //UG_PLUGIN_XBRAIDFORUG4_CACHEDTHETAINTEGRATORFACTORY_H
