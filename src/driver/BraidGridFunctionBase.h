//
// Created by parnet on 28.05.21.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_BRAIDGRIDFUNCTIONBASE_H
#define UG_PLUGIN_XBRAIDFORUG4_BRAIDGRIDFUNCTIONBASE_H

// global includes
#include <iomanip>
#include <sstream>
// project includes
#include <ugbase.h>
#include "common/math/math_vector_matrix/math_vector_functions.h"
#include "common/serialization.h"
#include "lib_disc/spatial_disc/domain_disc.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/dof_manager/function_pattern.h"
#include "lib_disc/time_disc/theta_time_step.h"
#include "lib_algebra/vector_interface/vec_functions.h"
#include "lib_algebra/operator/interface/linear_operator_inverse.h"
// plugin includes
#include "../../../Limex/time_disc/time_integrator.hpp"
#include "../../../Limex/time_disc/linear_implicit_timestep.h"
// library includes
#include "../../libs/braid/braid/braid.hpp"
// local includes - core
#include "../core/BraidVectorStruct.h"
#include "../core/SpaceTimeCommunicator.h"
// local includes - util
#include "../util/MemoryObserver.h"
#include "../util/BraidTimer.h"
#include "../util/BraidUsageTimer.h"
#include "../util/paralog.h"
#include "../util/Scriptor.h"
// local includes - interface
#include "../interface/scriptor.h"
#include "../interface/initializer.h"
#include "../interface/spatial_norm.h"


template<typename TDomain, typename TAlgebra>
class BraidGridFunctionBase : public BraidApp  {
public: // todo set better modes
/* ---------------------------------------------------------------------------------------------------------------------
 * Type definitions
 -------------------------------------------------------------------------------------------------------------------- */
    typedef ug::IDomainDiscretization <TAlgebra> TDomainDisc;
    typedef SmartPtr <TDomainDisc> SPDomainDisc;
    SPDomainDisc m_domain_disc;


    typedef SmartPtr<SpaceTimeCommunicator> SPSpaceTimeCommunicator;
    typedef ug::GridFunction<TDomain, TAlgebra> TGridFunction;
    typedef SmartPtr<TGridFunction> SPGridFunction;

    typedef typename TAlgebra::vector_type::value_type TVectorValueType;

    typedef ug::ITimeIntegrator<TDomain, TAlgebra> TTimeIntegrator;
    typedef SmartPtr<TTimeIntegrator> SPTimeIntegrator;
    typedef SmartPtr<BraidInitializer<TDomain,TAlgebra>> SPBraidInitializer;

    typedef Scriptor<TDomain, TAlgebra> TScriptor;
    typedef SmartPtr<TScriptor> SPScriptor;

    typedef Paralog TParalog;
    typedef SmartPtr<TParalog> SPParalog;

    typedef BraidSpatialNorm<TDomain,TAlgebra> TSpatialNorm;
    typedef SmartPtr<TSpatialNorm> SPSpatialNorm;

    SmartPtr<VTKScriptor<TDomain,TAlgebra>> vtk_scriptor;

    SmartPtr<VTKScriptor<TDomain,TAlgebra>> vtk_ustart_before;
    SmartPtr<VTKScriptor<TDomain,TAlgebra>> vtk_ustart_after;
    SmartPtr<VTKScriptor<TDomain,TAlgebra>> vtk_uend_after;
    SmartPtr<VTKScriptor<TDomain,TAlgebra>> vtk_uend_before;

    SmartPtr<VTKScriptor<TDomain,TAlgebra>> vtk_resu_before;
    SmartPtr<VTKScriptor<TDomain,TAlgebra>> vtk_resu_after;
    SmartPtr<VTKScriptor<TDomain,TAlgebra>> vtk_resr_before;
    SmartPtr<VTKScriptor<TDomain,TAlgebra>> vtk_resr_after;

    bool write_input_and_outputs = false;





    SmartPtr<VTKScriptor<TDomain,TAlgebra>> vtk_norm;
/* ---------------------------------------------------------------------------------------------------------------------
 * Member Variables
 -------------------------------------------------------------------------------------------------------------------- */

    void set_vtk_scriptor(SmartPtr<VTKScriptor<TDomain,TAlgebra>> sp_scriptor){
        this->vtk_scriptor = sp_scriptor;
    }

    void set_vtk_ustart_before(SmartPtr<VTKScriptor<TDomain,TAlgebra>> sp_scriptor){
        this->vtk_ustart_before = sp_scriptor;
    }
    void set_vtk_ustart_after(SmartPtr<VTKScriptor<TDomain,TAlgebra>> sp_scriptor){
        this->vtk_ustart_after = sp_scriptor;
    }

    void set_vtk_uend_before(SmartPtr<VTKScriptor<TDomain,TAlgebra>> sp_scriptor){
        this->vtk_uend_before = sp_scriptor;
    }

    void set_vtk_uend_after(SmartPtr<VTKScriptor<TDomain,TAlgebra>> sp_scriptor){
        this->vtk_uend_after = sp_scriptor;
    }



    void set_vtk_resu_before(SmartPtr<VTKScriptor<TDomain,TAlgebra>> sp_scriptor){
        this->vtk_resu_before = sp_scriptor;
    }

    void set_vtk_resu_after(SmartPtr<VTKScriptor<TDomain,TAlgebra>> sp_scriptor){
        this->vtk_resu_after = sp_scriptor;
    }

    void set_vtk_resr_before(SmartPtr<VTKScriptor<TDomain,TAlgebra>> sp_scriptor){
        this->vtk_resr_before = sp_scriptor;
    }

    void set_vtk_resr_after(SmartPtr<VTKScriptor<TDomain,TAlgebra>> sp_scriptor){
        this->vtk_resr_after = sp_scriptor;
    }

    void set_vtk_norm(SmartPtr<VTKScriptor<TDomain,TAlgebra>> sp_scriptor){
        this->vtk_norm =  sp_scriptor;
    }


    SPParalog m_log;
    SPParalog m_script_log;

    const char * m_name = nullptr;
    SPSpaceTimeCommunicator m_comm;
    bool m_verbose = true;
    SPGridFunction m_u0; // for t = tstart

    int m_levels = 15;
    bool provide_residual = false;

    SPScriptor m_out;
    SPBraidInitializer m_initializer;

    SPTimeIntegrator m_spIntegratorC;
    SPTimeIntegrator m_spIntegratorF;

    SPSpatialNorm m_norm;

    BraidTimer m_timer;

    BraidTimeLogManager m_time_log_manager;


    /**
    * Note that this default constructor does not create a consistent object. The parameter t_comm (of type MPI_Comm)
    * for the temporal communication has to be set.
    */
    BraidGridFunctionBase() : BraidApp(nullptr, 0, 10, 10) {
        //this->m_log->o << "debug::BraidGridFunctionBase()" << std::endl<<std::flush;
        this->m_name = "Braid Gridfunction Base";
        //this->m_log->o << "Warning Braid Function Base :: default constructor was used!" << std::endl;
    }
    BraidGridFunctionBase(MPI_Comm mpi_temporal, double tstart, double tstop, int steps) : BraidApp(mpi_temporal, tstart,tstop, steps) {
        //std::cout << "debug::BraidGridFunctionBase[[args]]" << std::endl<<std::flush;
        this->m_name = "Braid Grid Function Base";
        //std::cout << "debug::BraidGridFunctionBase[[end]]" << std::endl<<std::flush;
    }
    ~BraidGridFunctionBase() override = default;

/* ---------------------------------------------------------------------------------------------------------------------
* XBraid Method definitions
* ------------------------------------------------------------------------------------------------------------------- */
    braid_Int Init(braid_Real t, braid_Vector *u_ptr) override {
        //this->m_log->o << "debug::BraidGridFunctionBase::Init[[args]]" << std::endl<<std::flush;
        StartOperationTimer(Observer::T_INIT);
        auto *u = (BraidVector *) malloc(sizeof(BraidVector));
        auto *vec = new SPGridFunction();
        m_initializer->init(*vec, t);

        u->value = vec;
        u->time = t;


        this->m_script_log->o << "u_" << indexpool << " = init(" << t << ")" << std::endl;
        u->index = indexpool;
        indexpool++;

        *u_ptr = u;

        // ----> test init norm
        // SPGridFunction tempobject = vec->get()->clone();
        // auto norm_ptr = tempobject->norm();
        // this->m_log->o << "\t\tbraid init u norm: " << norm_ptr << std::endl; // todo chech for non 0 values
        // ----<


        StopOperationTimer(Observer::T_INIT);
        //this->m_log->o << "debug::BraidGridFunctionBase::Init[[end]]" << std::endl<<std::flush;
        return 0;
    };



    braid_Int Clone(braid_Vector u_, braid_Vector *v_ptr) override {
        //this->m_log->o << "debug::BraidGridFunctionBase::Clone[[args]]" << std::endl<<std::flush;
        StartOperationTimer(Observer::T_CLONE);
        //this->m_log->o << "debug::BraidGridFunctionBase::Clone[unpack->]" << std::endl<<std::flush;

        auto *v = (BraidVector *) malloc(sizeof(BraidVector));

        auto *uref = (SPGridFunction *) u_->value;
        auto *vref = new SPGridFunction();
        //this->m_log->o << "debug::BraidGridFunctionBase::Clone[clone]" << std::endl<<std::flush;
        *vref = uref->get()->clone();
        v->value = vref;
        v->time = u_->time;


        this->m_script_log->o << "u_" << indexpool << " = clone(u_" << u_->index << ")" << std::endl;

        v->index = indexpool;
        indexpool++;


        *v_ptr = v;

        // ----> test clone norm
        //SPGridFunction tempobject = vref->get()->clone();
        //auto norm_ptr = tempobject->norm();
        //this->m_log->o << "\t\tbraid clone u norm: " << norm_ptr << std::endl; // todo chech for non 0 values
        // ----<

        StopOperationTimer(Observer::T_CLONE);
        //this->m_log->o << "debug::BraidGridFunctionBase::Clone[[end]]" << std::endl<<std::flush;
        return 0;
    };



    braid_Int Free(braid_Vector u_) override {
        //this->m_log->o << "debug::BraidGridFunctionBase::Free[[args]]" << std::endl<<std::flush;
        StartOperationTimer(Observer::T_FREE);


        this->m_script_log->o  << "u_" << u_->index << " = null" << std::endl;

        auto *u_value = (SPGridFunction *) u_->value;
        delete u_value;
        free(u_);

        StopOperationTimer(Observer::T_FREE);
        //this->m_log->o << "debug::BraidGridFunctionBase::Free[[end]]" << std::endl<<std::flush;
        return 0;
    };



    braid_Int Sum(braid_Real alpha,braid_Vector x_,braid_Real beta,braid_Vector y_) override {
        //this->m_log->o << "debug::BraidGridFunctionBase::Sum[[args]]" << std::endl<<std::flush;
        StartOperationTimer(Observer::T_SUM);

            if (alpha == 0) {
                this->m_script_log->o  << "u_" << y_->index << " = " << beta << "* u_" << y_->index << " % Skalierung "
                                  << std::endl;
            } else if (beta == 0) {
                this->m_script_log->o  << "u_" << y_->index << " = " << alpha << "*u_" << x_->index << "  % Ersetzung "
                                  << std::endl;
            } else {
                this->m_script_log->o  << "u_" << y_->index << " = " << alpha << "* u_" << x_->index << "  + " << beta
                                  << "* u_"
                                  << y_->index << " % Summe " << std::endl;
            }


        auto *xref = (SPGridFunction *) x_->value;
        auto *yref = (SPGridFunction *) y_->value;
        //this->m_log->o << "vec add " << x_->time << "\t" << y_->time << std::endl;
        VecAdd(beta, *yref->get(), alpha, *xref->get());



        // this->m_log->o << "adjust - time = " << y_->time << std::endl;
        StopOperationTimer(Observer::T_SUM);


        //this->m_log->o << "debug::BraidGridFunctionBase::Sum[[end]]" << std::endl<<std::flush;
        return 0;
    };



    braid_Int SpatialNorm(braid_Vector u_, braid_Real *norm_ptr) override {
        *norm_ptr = 0;
        //this->m_log->o << "debug::BraidGridFunctionBase::SpatialNorm[[args]]" << std::endl<<std::flush;
        auto *uref = (SPGridFunction *) u_->value;
        SPGridFunction tempobject = uref->get()->clone(); // clone to ensure consistency
        *norm_ptr = m_norm->norm(tempobject);

        {
            SPGridFunction tempobject_n = uref->get()->clone();
            vtk_norm->write(tempobject_n, 0, 0, 0, 0);
        }

        this->m_script_log->o  << "norm( u_" << u_->index << ") % value ="  << *norm_ptr << std::endl;

        //this->m_log->o << "debug::BraidGridFunctionBase::SpatialNorm[[end]]" << std::endl<<std::flush;
        return 0;
    };



    braid_Int Access(braid_Vector u_, BraidAccessStatus &astatus) override {
        //this->m_log->o << "debug::BraidGridFunctionBase::Access[[args]]" << std::endl<<std::flush;
        //print_status(this->m_log->o,astatus);

        StartOperationTimer(Observer::T_ACCESS);




        int v;
        int index;
        astatus.GetTIndex(&index);
        double timestamp;
        astatus.GetT(&timestamp);

        auto ref = ((SPGridFunction *) u_->value)->get()->clone();


        int iter;
        int lvl;
        int done;

        astatus.GetIter(&iter);
        astatus.GetLevel(&lvl);
        astatus.GetDone(&done);


        if (done == 1) {
            this->m_script_log->o << "access( u_" << u_->index<< ")    @t="<<u_->time<< std::endl;
            v = this->m_out->write(ref, index, timestamp);
        } else {
            this->m_script_log->o << "access( u_" << u_->index<< ")    @t="<<u_->time<<"  filename=access_k"<<iter<<"_l"<<lvl<<"_t"<<std::setw(4) << std::setfill('0') <<index << ".vtu"<< std::endl;
            v = this->m_out->write(ref, index, timestamp, iter, lvl);
        }

        StopOperationTimer(Observer::T_ACCESS);
        //this->m_log->o << "debug::BraidGridFunctionBase::Access[[end]]" << std::endl<<std::flush;
        return 0;
    };



    braid_Int BufSize(braid_Int *size_ptr, BraidBufferStatus &bstatus) override{
        *size_ptr = 2048;
        //this->m_log->o << "debug::BraidGridFunctionBase::BufSize[[args]]" << std::endl<<std::flush;
        //print_status(this->m_log->o,bstatus);
        *size_ptr = (sizeof(TVectorValueType) * (*this->m_u0).size() // actual vector size
                     + sizeof(size_t)) // index
                    + 2 * sizeof(int) // todo find out what this is :)
                    + sizeof(double);
        //this->m_log->o << "debug::BraidGridFunctionBase::BufSize[[end]]" << std::endl<<std::flush;
        //this->m_script_log->o << "buffer_size = " << *size_ptr << std::flush << std::endl;
        return 0;
    };



    braid_Int BufPack(braid_Vector u_, void *buffer, BraidBufferStatus &bstatus) override {
        //this->m_log->o << "debug::BraidGridFunctionBase::BufPack[[args]]" << std::endl<<std::flush;
        //print_status(this->m_log->o,bstatus);
        StartOperationTimer(Observer::T_SEND);

        this->m_script_log->o << "send(u_" << u_->index << ")" << std::endl << std::flush;


        int bufferSize = 0; // startposition of gridfunction (will be written first) in buffer
        // delegate writing of gridfunction
        auto *u_ref = (SPGridFunction *) u_->value;
        pack(buffer, u_ref->get(), &bufferSize); // bufferSize returns position of bufferpointer after writing the gridfunction

        // write rank of source processor
        char *chBuffer = (char *) buffer;
        int temprank = this->m_comm->get_temporal_rank();
        memcpy(chBuffer + bufferSize, &temprank, sizeof(int));
        bufferSize += sizeof(int);

        // write index of vector

        memcpy(chBuffer + bufferSize, &u_->index, sizeof(int));
        bufferSize += sizeof(int);


        // write time for gridfunction
        memcpy(chBuffer + bufferSize, &u_->time, sizeof(double));
        bufferSize += sizeof(double );
        this->m_log->o  << "time written: " << u_->time << std::endl;
        bstatus.SetSize(bufferSize);

        StopOperationTimer(Observer::T_SEND);

        double diff, total;
        this->m_timer.now(total, diff);
        this->m_script_log->o << std::setw(10) << "@time:"
                          << std::setw(12) << total << " ; "
                          << std::setw(12) << diff << " Vector Send" << std::endl;

        //this->m_log->o << "debug::BraidGridFunctionBase::BufPack[[end]]" << std::endl<<std::flush;
        return 0;
    };



    braid_Int BufUnpack(void *buffer, braid_Vector *u_ptr, BraidBufferStatus &bstatus) override {
        //this->m_log->o << "debug::BraidGridFunctionBase::BufUnPack[[args]]" << std::endl<<std::flush;
        //print_status(this->m_log->o,bstatus);

        double diff, total;
        this->m_timer.now(total, diff);
        this->m_script_log->o << std::setw(10) << "@time:"
                          << std::setw(12) << total << " ; "
                          << std::setw(12) << diff << " Vector Received" << std::endl;

        StartOperationTimer(Observer::T_RECV);
        int pos = 0 ; // startposition of gridfunction (will be read first) in buffer
        //BufSize(&bufferSize, bstatus);
        auto *u = (BraidVector *) malloc(sizeof(BraidVector));
        auto *sp_u = new SPGridFunction(new TGridFunction(*this->m_u0)); // todo

        // delegate reading of gridfunction
        unpack(buffer, sp_u->get(), &pos);// pos returns position of bufferpointer after writing the gridfunction
        u->value = sp_u;

        // read rank of source processor
        char *chBuffer = (char *) buffer;
        int temprank;
        memcpy(&temprank, chBuffer + pos, sizeof(int));
        pos += sizeof(int);

        // read temporal rank of source processor

        if (this->m_verbose) {
            this->m_script_log->o << "u_" << indexpool << " = ";
        }
        int index;
        memcpy(&index, chBuffer + pos, sizeof(int));
        pos += sizeof(int);
        this->m_script_log->o << "rec( v_" << temprank << "_" << index << ")" << std::endl;


        double t;
        memcpy(&t, chBuffer + pos, sizeof(double));
        pos += sizeof(int);

        u->time = t;
        this->m_log->o  << "time read: " << u->time << std::endl;


        u->index = indexpool;
        indexpool++;

        *u_ptr = u;
        StopOperationTimer(Observer::T_RECV);
        //this->m_log->o << "debug::BraidGridFunctionBase::BufUnPack[[end]]" << std::endl<<std::flush;
        return 0;
    };



    braid_Int Sync(BraidSyncStatus &sstatus) override {
        //this->m_log->o << "debug::BraidGridFunctionBase::Sync[[args]]" << std::endl<<std::flush;
        //print_status(this->m_log->o,sstatus);
        // todo ?
        //this->m_log->o << "debug::BraidGridFunctionBase::Sync[[end]]" << std::endl<<std::flush;
        return 0;
    };



    braid_Int Coarsen(braid_Vector fu_,braid_Vector *cu_ptr,BraidCoarsenRefStatus &status) override {
        //this->m_log->o << "debug::BraidGridFunctionBase::Coarsen[[args]]" << std::endl<<std::flush;
        //print_status(this->m_log->o,status);
        return 0;
        this->Clone(fu_, cu_ptr);

        auto *sp_fu = (SPGridFunction *) fu_->value;
        auto *sp_cu = (SPGridFunction *) (*cu_ptr)->value;

        double t_upper;
        double t_lower;
        status.GetT(&t_lower);
        status.GetCTstop(&t_upper);
        //status.GetFTstop(&t_upper);

        m_spIntegratorC->apply(*sp_fu, t_upper, sp_cu->cast_const(), t_lower); // todo check fu, cu order
        //this->m_log->o << "debug::BraidGridFunctionBase::Coarsen[[end]]" << std::endl<<std::flush;
        return 0;
    };



    braid_Int Refine(braid_Vector cu_,braid_Vector *fu_ptr,BraidCoarsenRefStatus &status) override {
        //this->m_log->o << "debug::BraidGridFunctionBase::Refine[[args]]" << std::endl<<std::flush;
        //print_status(this->m_log->o,status);
        return 0;
        this->Clone(cu_, fu_ptr);

        auto *sp_fu = (SPGridFunction *) (*fu_ptr)->value;
        auto *sp_cu = (SPGridFunction *) cu_->value;

        double t_upper;
        double t_lower;
        status.GetT(&t_lower);
        status.GetCTstop(&t_upper);
        //status.GetFTstop(&t_upper);

        m_spIntegratorF->apply(*sp_cu, t_upper, sp_fu->cast_const(), t_lower); // todo check fu, cu order
        //this->m_log->o << "debug::BraidGridFunctionBase::Refine[[end]]" << std::endl<<std::flush;

        return 0;
    };
/* ---------------------------------------------------------------------------------------------------------------------
 * Methods
 * ------------------------------------------------------------------------------------------------------------------ */
    void set_paralog(SPParalog log){
        this->m_log = log;
    }

    void set_paralog_script(SPParalog log){
        this->m_script_log = log;
    }

    void set_script_log(SPParalog log){
        this->m_script_log = log;
    }

    void init() {
        this->m_log->init();
        //this->m_log->o << "debug::BraidGridFunctionBase::init[[args]]" << std::endl<<std::flush;
        this->m_timer.start(); // to trace send and receive times / MPI
        this->m_time_log_manager = BraidTimeLogManager(this->m_levels);


        //const ug::GridLevel gridlevel = this->m_u0->grid_level();
        //m_A = SPAssembledOperator(new TAssembledOperator(m_timeDisc, gridlevel));
        //this->m_log->o << "finished init" << std::endl << std::flush;
        //this->m_log->o << "debug::BraidGridFunctionBase::init[[end]]" << std::endl<<std::flush;
    }



    inline void pack(void *buffer, TGridFunction *u_ref, int *bufferSize) {
        //this->m_log->o << "debug::BraidGridFunctionBase::pack[[args]]" << std::endl<<std::flush;
        char *chBuffer = (char *) buffer;
        *bufferSize = 0;
        size_t szVector = u_ref->size();

        //write number of gridpoints
        memcpy(buffer, &szVector, sizeof(size_t));
        *bufferSize += sizeof(size_t);

        // write the value for each gridpoint
        for (size_t i = 0; i < szVector; i++) {
            memcpy(chBuffer + *bufferSize, &(*u_ref)[i], sizeof(TVectorValueType));
            *bufferSize += sizeof(TVectorValueType);
        }


        //this->m_log->o << "debug::BraidGridFunctionBase::pack[[end]]" << std::endl<<std::flush;
    }



    inline void unpack(void *buffer, TGridFunction *u_ref, int *pos) {
        //this->m_log->o << "debug::BraidGridFunctionBase::unpack[[args]]" << std::endl<<std::flush;
        char *chBuffer = (char *) buffer;
        size_t szVector = 0;

        // read number of gridpoints todo consistency check?
        memcpy(&szVector, chBuffer + *pos, sizeof(size_t));
        *pos = sizeof(size_t);

        // read the values for each gridpoint
        for (size_t i = 0; i < szVector; i++) {
            TVectorValueType val = 0;
            memcpy(&val, chBuffer + *pos, sizeof(TVectorValueType));
            *pos += sizeof(TVectorValueType);
            (*u_ref)[i] = val;
        }


        //this->m_log->o << "debug::BraidGridFunctionBase::unpack[[end]]" << std::endl<<std::flush;
    }



    void set_time_values(double startTime, double endTime, int n) {
        this->tstart = startTime;
        this->tstop = endTime;
        this->ntime = n;
    }



    void set_start_time(double startTime) {
        this->tstart = startTime;
    }



    void set_end_time(double endTime) {
        this->tstop = endTime;
    }



    void set_number_of_timesteps(int n) {
        this->ntime = n;
    }



    void set_verbose(bool verbose) {
        this->m_verbose = verbose;
    }



    void set_start_vector(SPGridFunction p_u0) {
        this->m_u0 = p_u0;
    }



    void set_scriptor(SPScriptor p_out) {
        this->m_out = p_out;
    }


    void set_norm_provider(SPSpatialNorm norm){
        this->m_norm = norm;
    }

    void set_max_levels(size_t levelcount) {
        if(levelcount > 10) {
            std::cerr << "Warning: Tracing times is currently limited to 15 level" << std::endl; // todo clean up
        }
        this->m_levels = levelcount;
    }



    void set_initializer(SPBraidInitializer initializer){
        this->m_initializer = initializer;
    }

    void set_domain(SPDomainDisc domain){
        this->m_domain_disc = domain;
    }
};


#endif //UG_PLUGIN_XBRAIDFORUG4_BRAIDGRIDFUNCTIONBASE_H
