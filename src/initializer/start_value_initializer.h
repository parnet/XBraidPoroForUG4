//
// Created by parnet on 12.06.21.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_START_VALUE_INITIALIZER_H
#define UG_PLUGIN_XBRAIDFORUG4_START_VALUE_INITIALIZER_H
template<typename TDomain, typename TAlgebra>
class StartValueInitializer : public BraidInitializer<TDomain,TAlgebra> {
public:
    typedef ug::GridFunction<TDomain, TAlgebra> TGridFunction;
    typedef SmartPtr<TGridFunction> SPGridFunction;
/* ---------------------------------------------------------------------------------------------------------------------
 * Member Variable
 * ------------------------------------------------------------------------------------------------------------------ */
    SPGridFunction m_u0; // start vector
    double t_start = 0.0;

    StartValueInitializer() = default;
    ~StartValueInitializer() = default;

    void init(SPGridFunction &u, double time ) override{
        u = this->m_u0->clone();
    }

    void set_start_vector(SPGridFunction u0){ // todo const
        this->m_u0 = u0;
    }
};

#endif //UG_PLUGIN_XBRAIDFORUG4_START_VALUE_INITIALIZER_H
