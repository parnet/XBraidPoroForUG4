//
// Created by parnet on 12.06.21.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_SCRIPTOR_H
#define UG_PLUGIN_XBRAIDFORUG4_SCRIPTOR_H

#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/function_spaces/interpolate.h"
#include "lib_disc/spatial_disc/domain_disc.h"


template<typename TDomain, typename TAlgebra>
class Scriptor {
public:
    typedef ug::GridFunction <TDomain, TAlgebra> TGridFunction;
    typedef SmartPtr <TGridFunction> SPGridFunction;

    Scriptor() = default;

    virtual ~Scriptor() = default;

    virtual bool write(SPGridFunction u, int index, double time) = 0;

    virtual bool write(SPGridFunction u, int index, double time, int iteration, int level) = 0 ;

    //virtual void print(const char *filename, TGridFunction &u, int index, double time) = 0;
};

#endif //UG_PLUGIN_XBRAIDFORUG4_SCRIPTOR_H
