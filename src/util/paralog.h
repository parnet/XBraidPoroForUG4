//
// Created by parnet on 06.06.21.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_PARALOG_H
#define UG_PLUGIN_XBRAIDFORUG4_PARALOG_H

#include <iostream>
#include <fstream>
#include "common/util/smart_pointer.h"
#include "../core/SpaceTimeCommunicator.h"

class Paralog {

public:
    typedef SmartPtr<SpaceTimeCommunicator> SPSpaceTimeCommunicator;
    int m_init = 0;
    std::ofstream o;
    const char * m_filename = "job";
    SPSpaceTimeCommunicator m_comm;

    void set_comm(SPSpaceTimeCommunicator comm){
        this->m_comm = comm;
    }

    void set_file_name(const char * filename){
        this->m_filename = filename;
    }

    void init() {
        if (this->m_init == 0) {
            std::stringstream ss;
            ss << this->m_filename << "_" << this->m_comm->get_temporal_rank() << ".output";
            this->o.open(ss.str());
        }
        this->m_init++;
    }

    void release(){
        this->m_init--;
        if(this->m_init == 0){
            // todo release
            this->o << std::endl << std::endl << "finished" << std::endl << std::flush;

        }
    }

    void write(const char * content){
        this->o << content << std::endl;
    }
};


#endif //UG_PLUGIN_XBRAIDFORUG4_PARALOG_H
