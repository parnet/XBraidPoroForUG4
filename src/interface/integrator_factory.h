//
// Created by parnet on 07.06.21.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_INTEGRATOR_FACTORY_H
#define UG_PLUGIN_XBRAIDFORUG4_INTEGRATOR_FACTORY_H

#include "../../../Limex/time_disc/time_integrator.hpp"

template<typename TDomain, typename TAlgebra>
class IntegratorFactory {
public:
    const char * m_name = "Interface";
    typedef ug::ITimeIntegrator<TDomain, TAlgebra> TTimeIntegrator;
    typedef SmartPtr<TTimeIntegrator> SPTimeIntegrator;

    virtual SPTimeIntegrator create_level_time_integrator(double current_dt, bool done,int level) = 0;

    virtual SPTimeIntegrator create_time_integrator(double current_dt, bool done) = 0;

    IntegratorFactory() = default;
    ~IntegratorFactory() = default;

    const char * get_name(){
        return this->m_name;
    }
};


#endif //UG_PLUGIN_XBRAIDFORUG4_INTEGRATOR_FACTORY_H
