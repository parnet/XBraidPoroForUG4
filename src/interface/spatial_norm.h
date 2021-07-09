//
// Created by parnet on 12.06.21.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_SPATIAL_NORM_H
#define UG_PLUGIN_XBRAIDFORUG4_SPATIAL_NORM_H


#include "lib_disc/function_spaces/grid_function.h"

template<typename TDomain, typename TAlgebra>
class BraidSpatialNorm {
public:
    typedef ug::GridFunction<TDomain, TAlgebra> TGridFunction;
    typedef SmartPtr<TGridFunction> SPGridFunction;

    virtual double norm(SPGridFunction u) = 0;
};


#endif //UG_PLUGIN_XBRAIDFORUG4_SPATIAL_NORM_H
