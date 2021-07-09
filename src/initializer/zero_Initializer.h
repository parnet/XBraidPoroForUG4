//
// Created by parnet on 12.06.21.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_ZERO_INITIALIZER_H
#define UG_PLUGIN_XBRAIDFORUG4_ZERO_INITIALIZER_H

template<typename TDomain, typename TAlgebra>
class ZeroValueInitializer : public BraidInitializer<TDomain,TAlgebra>{
public: // todo set better modes
    typedef ug::GridFunction<TDomain, TAlgebra> TGridFunction;
    typedef SmartPtr<TGridFunction> SPGridFunction;

/* ---------------------------------------------------------------------------------------------------------------------
 * Member Variable
 -------------------------------------------------------------------------------------------------------------------- */
    SPGridFunction m_u0; // start vector
    number t_start;

    void init(SPGridFunction &u, number time ) override {
        if(time == this->t_start){
            u = this->m_u0->clone();
        } else {
            u = this->m_u0->clone_without_values();

        }
    }

    void set_start_time(number tstart){
        this->t_start = tstart;
    }

    void set_start_vector(SPGridFunction  u0){
        this->m_u0 = u0;
    }
};

#endif //UG_PLUGIN_XBRAIDFORUG4_ZERO_INITIALIZER_H
