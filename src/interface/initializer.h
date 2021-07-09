//
// Created by parnet on 12.06.21.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_INITIALIZER_H
#define UG_PLUGIN_XBRAIDFORUG4_INITIALIZER_H
#include "lib_disc/function_spaces/grid_function.h"

template<typename TDomain, typename TAlgebra>
class BraidInitializer {
public: // todo set better modes
/* ---------------------------------------------------------------------------------------------------------------------
 * Type definitions
 * ------------------------------------------------------------------------------------------------------------------ */
    typedef ug::GridFunction<TDomain, TAlgebra> TGridFunction;
    typedef SmartPtr<TGridFunction> SPGridFunction;


/* ---------------------------------------------------------------------------------------------------------------------
 * Methods
 * ------------------------------------------------------------------------------------------------------------------ */
    virtual void init(SPGridFunction &u, number time ) = 0;
};


// *vec = this->m_u0->clone_without_values(); // todo
// Interpolate(this->m_data, *vec, this->m_cmp, NULL, t); // todo
// m_domainDisc->adjust_solution(*vec->get(), t); // todo


#endif //UG_PLUGIN_XBRAIDFORUG4_INITIALIZER_H
