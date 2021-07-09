//
// Created by parnet on 12.06.21.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_EUCLIDIAN_NORM_H
#define UG_PLUGIN_XBRAIDFORUG4_EUCLIDIAN_NORM_H

#include "../interface/spatial_norm.h"

template<typename TDomain, typename TAlgebra>
class BraidEuclidianNorm  : public BraidSpatialNorm<TDomain,TAlgebra>{
public:
    typedef ug::GridFunction<TDomain, TAlgebra> TGridFunction;
    typedef SmartPtr<TGridFunction> SPGridFunction;



    BraidEuclidianNorm() : BraidSpatialNorm<TDomain,TAlgebra>(){}

    ~BraidEuclidianNorm() = default;

    double norm(SPGridFunction u) override {
        return u->norm();
    }
};
#endif //UG_PLUGIN_XBRAIDFORUG4_EUCLIDIAN_NORM_H
