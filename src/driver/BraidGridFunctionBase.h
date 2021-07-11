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
    SmartPtr<VTKScriptor<TDomain,TAlgebra>> vtkScriptor;
/* ---------------------------------------------------------------------------------------------------------------------
 * Member Variables
 -------------------------------------------------------------------------------------------------------------------- */

    void set_vtk_scriptor(SmartPtr<VTKScriptor<TDomain,TAlgebra>> sp_scriptor){
        this->vtkScriptor = sp_scriptor;
    }

    SPParalog m_log;
    const char * m_name = nullptr;
    SPSpaceTimeCommunicator m_comm;
    bool m_verbose = true;
    SPGridFunction m_u0; // for t = tstart
    std::ofstream debugwriter;
    int m_levels = 15;
    bool provide_residual = false;

    SPScriptor m_out;
    SPBraidInitializer m_initializer;

    SPTimeIntegrator m_spIntegratorC;
    SPTimeIntegrator m_spIntegratorF;

    SPSpatialNorm m_norm;

    BraidTimer m_timer;
#if TRACE_TIMINGS == 1
    BraidTimeLogManager m_time_log_manager;
#endif
#if TRACE_GRIDFUNCTION == 1
    SPMATLABScriptor matlab;
#endif
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
#if TRACE_INDEX == 1
        if (this->m_verbose) {
            this->debugwriter << "u_" << indexpool << " = init(" << t << ")" << std::endl;
        }
#endif
        StartOperationTimer(Observer::T_INIT);
        auto *u = (BraidVector *) malloc(sizeof(BraidVector));
        auto *vec = new SPGridFunction();
        m_initializer->init(*vec, t);
        u->value = vec;

#if TRACE_INDEX == 1
        u->index = indexpool;
        indexpool++;
        MATLAB(vec->get(), u->index, t);
#endif
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
#if TRACE_INDEX == 1
        if (this->m_verbose) {
            this->debugwriter << "u_" << indexpool << " = clone(u_" << u_->index << ")" << std::endl;
        }
#endif
        //this->m_log->o << "debug::BraidGridFunctionBase::Clone[unpack->]" << std::endl<<std::flush;

        auto *v = (BraidVector *) malloc(sizeof(BraidVector));

        auto *uref = (SPGridFunction *) u_->value;
        auto *vref = new SPGridFunction();
        //this->m_log->o << "debug::BraidGridFunctionBase::Clone[clone]" << std::endl<<std::flush;
        *vref = uref->get()->clone();
        v->value = vref;

#if TRACE_INDEX == 1
        v->index = indexpool;
        indexpool++;
        MATLAB(vref->get(), v->index, -1.0);
#endif

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

#if TRACE_INDEX == 1
        if (this->m_verbose) {
            this->debugwriter << "u_" << u_->index << " = null" << std::endl;
        }
#if TRACE_CONST == 1
        if (u->m_const) {
        this->o << "u_" << u_->index << " was const" << std::endl;
        const_free++;
    }
#endif // Trace Index
#endif // trace const
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
#if TRACE_INDEX == 1
        if (this->m_verbose) {
            if (alpha == 0) {
                this->debugwriter << "u_" << y_->index << " = " << beta << "* u_" << y_->index << " % Scale "
                                  << std::endl;
            } else if (beta == 0) {
                this->debugwriter << "u_" << y_->index << " = " << alpha << "*u_" << x_->index << "  % Replace "
                                  << std::endl;
            } else {
                this->debugwriter << "u_" << y_->index << " = " << alpha << "* u_" << x_->index << "  + " << beta
                                  << "* u_"
                                  << y_->index << " % Sum " << std::endl;
            }
        }
#endif
#if TRACE_CONST == 1
        y->m_const = false;
#endif
        auto *xref = (SPGridFunction *) x_->value;
        auto *yref = (SPGridFunction *) y_->value;
        VecAdd(beta, *yref->get(), alpha, *xref->get());
        StopOperationTimer(Observer::T_SUM);
#if TRACE_INDEX == 1
        MATLAB(yref->get(), y->index, -1.0);
#endif
        //this->m_log->o << "debug::BraidGridFunctionBase::Sum[[end]]" << std::endl<<std::flush;
        return 0;
    };



    braid_Int SpatialNorm(braid_Vector u_, braid_Real *norm_ptr) override {
        *norm_ptr = 0;
        //this->m_log->o << "debug::BraidGridFunctionBase::SpatialNorm[[args]]" << std::endl<<std::flush;
        auto *uref = (SPGridFunction *) u_->value;
        SPGridFunction tempobject = uref->get()->clone(); // clone to ensure consistency
        *norm_ptr = m_norm->norm(tempobject);
#if TRACE_INDEX == 1
        if (this->m_verbose) {
            this->debugwriter << "norm( u_" << u_->index << ") % " << *norm_ptr << std::endl;
        }
#endif
        //this->m_log->o << "debug::BraidGridFunctionBase::SpatialNorm[[end]]" << std::endl<<std::flush;
        return 0;
    };



    braid_Int Access(braid_Vector u_, BraidAccessStatus &astatus) override {
        //this->m_log->o << "debug::BraidGridFunctionBase::Access[[args]]" << std::endl<<std::flush;
        //print_status(this->m_log->o,astatus);

        StartOperationTimer(Observer::T_ACCESS);
#if TRACE_INDEX == 1
        if (this->m_verbose) {
            this->debugwriter << "% \t Access \t" << u_->index << std::endl;
        }
#endif

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
            v = this->m_out->write(ref, index, timestamp);
        } else {
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
                    + 2 * sizeof(int); // todo find out what this is :)
        //this->m_log->o << "debug::BraidGridFunctionBase::BufSize[[end]]" << std::endl<<std::flush;
        return 0;
    };



    braid_Int BufPack(braid_Vector u_, void *buffer, BraidBufferStatus &bstatus) override {
        //this->m_log->o << "debug::BraidGridFunctionBase::BufPack[[args]]" << std::endl<<std::flush;
        //print_status(this->m_log->o,bstatus);
        StartOperationTimer(Observer::T_SEND);

        int bufferSize;
#if TRACE_INDEX == 1
        if (this->m_verbose) {
            debugwriter << "send(u_" << u_->index << ")" << std::endl << std::flush;
        }
#endif

        auto *u_ref = (SPGridFunction *) u_->value;


        pack(buffer, u_ref->get(), &bufferSize);

#if TRACE_INDEX == 1
        char *chBuffer = (char *) buffer;
        memcpy(chBuffer + bufferSize, &u_->index, sizeof(int));
        bufferSize += sizeof(int);
#endif
        bstatus.SetSize(bufferSize);

        StopOperationTimer(Observer::T_SEND);
#if TRACE_RECVTIME == 1
        double diff, total;
        this->m_timer.now(total, diff);
        this->debugwriter << std::setw(10) << "@time:"
                          << std::setw(12) << total << " ; "
                          << std::setw(12) << diff << " Vector Send" << std::endl;
#endif
        //this->m_log->o << "debug::BraidGridFunctionBase::BufPack[[end]]" << std::endl<<std::flush;
        return 0;
    };



    braid_Int BufUnpack(void *buffer, braid_Vector *u_ptr, BraidBufferStatus &bstatus) override {
        //this->m_log->o << "debug::BraidGridFunctionBase::BufUnPack[[args]]" << std::endl<<std::flush;
        //print_status(this->m_log->o,bstatus);
#if TRACE_RECVTIME == 1
        double diff, total;
        this->m_timer.now(total, diff);
        this->debugwriter << std::setw(10) << "@time:"
                          << std::setw(12) << total << " ; "
                          << std::setw(12) << diff << " Vector Received" << std::endl;
#endif
        StartOperationTimer(Observer::T_RECV);
#if TRACE_INDEX == 1
        if (this->m_verbose) {
            this->debugwriter << "u_" << indexpool << " = ";
        }
#endif

        int bufferSize;
        BufSize(&bufferSize, bstatus);
        auto *u = (BraidVector *) malloc(sizeof(BraidVector));
        auto *sp_u = new SPGridFunction(new TGridFunction(*this->m_u0)); // todo
        unpack(buffer, sp_u->get(), &bufferSize);
        u->value = sp_u;
#if TRACE_INDEX == 1
        u->index = indexpool;
        indexpool++;
        MATLAB(sp_u->get(), u->index, -1.0);
#endif
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



    void init() {
        this->m_log->init();
        //this->m_log->o << "debug::BraidGridFunctionBase::init[[args]]" << std::endl<<std::flush;


        this->m_timer.start(); // to trace send and receive times / MPI
#if TRACE_TIMINGS == 1
        this->m_time_log_manager = BraidTimeLogManager(this->m_levels);
#endif
#if TRACE_GRIDFUNCTION == 1
        this->matlab = SmartPtr<MATLABScriptor<TDomain, TAlgebra>>(new MATLABScriptor<TDomain, TAlgebra>(this->o));
#endif
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

        memcpy(buffer, &szVector, sizeof(size_t));
        *bufferSize += sizeof(size_t);
        for (size_t i = 0; i < szVector; i++) {
            memcpy(chBuffer + *bufferSize, &(*u_ref)[i], sizeof(TVectorValueType));
            *bufferSize += sizeof(TVectorValueType);
        }

        int temprank = this->m_comm->get_temporal_rank();

        memcpy(chBuffer + *bufferSize, &temprank, sizeof(int));
        *bufferSize += sizeof(int);
        //this->m_log->o << "debug::BraidGridFunctionBase::pack[[end]]" << std::endl<<std::flush;
    }



    inline void unpack(void *buffer, TGridFunction *u_ref, int *bufferSize) {
        //this->m_log->o << "debug::BraidGridFunctionBase::unpack[[args]]" << std::endl<<std::flush;
        char *chBuffer = (char *) buffer;
        size_t szVector = 0;
        memcpy(&szVector, chBuffer, sizeof(size_t));
        int pos = sizeof(size_t);
        for (size_t i = 0; i < szVector; i++) {
            TVectorValueType val = 0;
            memcpy(&val, chBuffer + pos, sizeof(TVectorValueType));
            pos += sizeof(TVectorValueType);
            (*u_ref)[i] = val;
        }

        int temprank;
        memcpy(&temprank, chBuffer + pos, sizeof(int));
        pos += sizeof(int);

#if TRACE_INDEX == 1
        int index;
        memcpy(&index, chBuffer + pos, sizeof(int));
        pos += sizeof(int);
        if (this->m_verbose) {
            debugwriter << "rec( v_" << temprank << "_" << index << ")" << std::endl;
        }
#endif
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
#if TRACE_TIMINGS == 1
        if(levelcount > 10) {
            std::cerr << "Warning: Tracing times is currently limited to 15 level" << std::endl; // todo clean up
        }
#endif
        this->m_levels = levelcount;
    }



    void set_initializer(SPBraidInitializer initializer){
        this->m_initializer = initializer;
    }
};


#endif //UG_PLUGIN_XBRAIDFORUG4_BRAIDGRIDFUNCTIONBASE_H
