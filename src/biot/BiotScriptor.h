//
// Created by parnet on 08.06.21.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_BIOTSCRIPTOR_H
#define UG_PLUGIN_XBRAIDFORUG4_BIOTSCRIPTOR_H

#include "../util/Scriptor.h"
#include "../../../Poroelasticity/barry_mercer.h"


template<typename TDomain, typename TAlgebra>
class EvalScriptor : public Scriptor<TDomain, TAlgebra> {
public:
    typedef ug::GridFunction <TDomain, TAlgebra> TGridFunction;
    typedef SmartPtr <TGridFunction> SPGridFunction;
    //typedef SmartPtr <ug::VTKOutput<TDomain::dim>> SPVTKOutput;
    //typedef SmartPtr <ug::IDomainDiscretization<TAlgebra>> SPDomainDisc;

    //typedef SmartPtr <ug::UserData<double, TGridFunction::dim>> SPData;
    //SPData m_data;

    //const char *m_cmp;

    //bool m_relative = false;
    std::ofstream *outfile;
    //SPDomainDisc m_domainDisc;

    EvalScriptor() : Scriptor<TDomain, TAlgebra>() {

    }

    ~EvalScriptor() {
        //outfile->close();
        //delete outfile;
    }

    void  set_problem(SmartPtr<BarryMercerProblem> problem){

    }

    void set_file(const char *filename) {
        outfile = new std::ofstream();
        outfile->open(filename);
        (*this->outfile) << "iteration; level; index; time; || u - v ||;" ;
        if(this->m_relative){
            (*this->outfile) << "  || u ||; relative;";
        }
        (*this->outfile) << std::endl;
    }

    void set_file(std::ofstream &file) {
        this->outfile = &file;
    }

    void write_time_pvd(TGridFunction *u) {};

    bool write(TGridFunction *u, int index, double time) override {
        return this->write(u,index,time,-1,-1);
    };

    bool write(TGridFunction *u, int index, double time, int iteration, int level) override  {
        SPGridFunction vec = u->clone_without_values();

        VecAdd(-1.0, *vec.get(), 1.0, *u);
        if(iteration == -1){
            (*this->outfile) << std::setw(4) << ""<<";"<< std::setw(4) << "" <<";";
        } else {
            (*this->outfile) << std::setw(4) << iteration << ";"<< std::setw(4)  << level << ";";
        }

        double vecnorm = vec->norm();
        (*this->outfile) << std::setw(7) << index <<";"
                         << std::setw(12) << time << ";"
                         << std::setw(12) << vecnorm << ";";
        if(this->m_relative){
            SPGridFunction uvec = u->clone();
            double unorm = uvec->norm();
            (*this->outfile) << std::setw(12) << unorm << ";"
                             << std::setw(12) << (vecnorm / unorm) << ";";
        }

        (*this->outfile) << std::endl;

        return true;
    };

    void print(const char *filename, TGridFunction &u, int index, double time) override {
        TGridFunction *uptr = &u;
        this->write(uptr, index, time);
    }

    void write_time_pvd_fn(const char *filename, TGridFunction &u) {
    }


    void set_generator_component(const char *cmp) {
        this->m_cmp = cmp;
    }


    void set_vector_generator(SPData p_data) {
        this->m_data = p_data;
    }

    void set_problem_parameter(SPDomainDisc p_domainDisc) {
        m_domainDisc = p_domainDisc;
    }
};

#endif //UG_PLUGIN_XBRAIDFORUG4_BIOTSCRIPTOR_H
