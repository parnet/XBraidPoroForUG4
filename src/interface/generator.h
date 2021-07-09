//
// Created by parnet on 27.05.21.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_IGENERATOR_H
#define UG_PLUGIN_XBRAIDFORUG4_IGENERATOR_H

#include "../core/SpaceTimeCommunicator.h"
#include "../util/Scriptor.h"
#include "../util/MemoryObserver.h"
#include "../util/BraidTimer.h"
#include "../util/BraidUsageTimer.h"
#include <iomanip>
#include <sstream>

#include <ugbase.h>
/*
#include "common/math/math_vector_matrix/math_vector_functions.h"
#include "common/serialization.h"
#include "lib_disc/spatial_disc/domain_disc.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/time_disc/theta_time_step.h"
#include "lib_algebra/vector_interface/vec_functions.h"
#include "lib_algebra/operator/interface/linear_operator_inverse.h"
#include "../../../Limex/time_disc/time_integrator.hpp"
#include "../../../Limex/time_disc/linear_implicit_timestep.h"
#include "lib_disc/dof_manager/function_pattern.h"

template<typename TDomain, typename TAlgebra>
class IGenerator {
public: // todo set better modes
    typedef ug::GridFunction<TDomain, TAlgebra> TGridFunction;
    typedef SmartPtr <ug::UserData<double, TGridFunction::dim>> SPData;  // Funktion zum Interpolieren / Generieren von außen
    typedef SmartPtr <ug::IDomainDiscretization<TAlgebra>>  SPDomain;

    const char *m_cmp;
    SPData m_data;
    SPDomain m_domainDisc;


    void setGeneratorComponent(const char *cmp) {
        this->m_cmp = cmp;
    }

    void setVectorGenerator(SPData p_data) {
        this->m_data = p_data;
    }
    void setDomainDisc(SPDomain p_domainDisc) {
        m_domainDisc = p_domainDisc;
    }
#ifdef UG_FOR_LUA
    void setVectorGenerator(const char *fctName) {
        setVectorGenerator(ug::LuaUserDataFactory<double, TDomain::dim>::create(fctName));
    }

    void setVectorGenerator(ug::LuaFunctionHandle fct) {
        setVectorGenerator(make_sp(new ug::LuaUserData<double, TDomain::dim>(fct)));
    }
#endif

};
 */


#endif //UG_PLUGIN_XBRAIDFORUG4_IGENERATOR_H
