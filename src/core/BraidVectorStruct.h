//
// Created by parnet on 12.05.19.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_BRAIDVECTORSTRUCT_H
#define UG_PLUGIN_XBRAIDFORUG4_BRAIDVECTORSTRUCT_H

#include "../../libs/braid/braid/braid.hpp"


#define WRITE_SCRIPT 1
/**
 * To display the time points to which the BufferUnpack and BufferPack functions are requested to
 */
#define TRACE_RECVTIME 1
#define TRACE_TIMINGS 1

#if TRACE_TIMINGS == 1
    #define StartOperationTimer(opt) this->m_time_log_manager.get(opt).start()
    #define StopOperationTimer(opt) this->m_time_log_manager.get(opt).stop()
    #define StartLevelOperationTimer(opt, l) this->m_time_log_manager.get(opt,l).start()
    #define StopLevelOperationTimer(opt, l) this->m_time_log_manager.get(opt,l).stop()
#else
    #define StartOperationTimer(opt)
    #define StopOperationTimer(opt)
    #define StartLevelOperationTimer(opt, l)
    #define StopLevelOperationTimer(opt, l)
#endif


/*
 * this structure represents the grid functions that are used within xbraid as vector types. It is possible to create
 * a class which stores additional information in the void pointer to use this structure within different context without
 * changing it.
 */
typedef struct _braid_Vector_struct {
    void *value{};

    double time = 0.0;

#if WRITE_SCRIPT == 1
    size_t index = 0;
#endif

} BraidVector;

#if WRITE_SCRIPT == 1
    size_t indexpool = 0;
#endif



void print_status(std::ofstream & printer, BraidAccessStatus &status) {
    printer << " ===== Braid Access Status ==== [[ " << std::endl;
    {
        double t_ptr;
        int iter_ptr;
        int level_ptr;
        int done_ptr;
        status.GetTILD(&t_ptr,&iter_ptr,&level_ptr,&done_ptr);

        printer << "\tt: " << t_ptr << std::endl;
        printer << "\titer: " << iter_ptr << std::endl;
        printer << "\tlevel: " << level_ptr << std::endl;
        printer << "\tdone: " << done_ptr << std::endl;

    }
    /*{
        double t_ptr;
        status.GetT(&t_ptr);

        printer << "\tt: " << t_ptr << std::endl;
    }*/
    {
        int tindex_ptr;
        status.GetTIndex(&tindex_ptr);

        printer << "\tt_idx: " << tindex_ptr << std::endl;
    }
    /*{
        int done_ptr;
        status.GetDone(&done_ptr);
        printer << "\tdone: " << done_ptr << std::endl;
    }*/
    /*{
        int level_ptr;
        status.GetLevel(&level_ptr);
        printer << "\tlevel: " << level_ptr << std::endl;
    }*/
    {
        int nlevels_ptr;
        status.GetNLevels(&nlevels_ptr);
        printer << "\tn_level: " << nlevels_ptr << std::endl;
    }
    /*{
        int iter_ptr;
        status.GetIter(&iter_ptr);
        printer << "\titer: " << iter_ptr << std::endl;
    }*/
    {
        int wtest_ptr;
        status.GetWrapperTest(&wtest_ptr);
        printer << "\twtest: " << wtest_ptr << std::endl;
    }
    {
        double rnorm_ptr;
        status.GetResidual(&rnorm_ptr);
        printer << "\tresidual: " << rnorm_ptr << std::endl;
    }
    {
        int nrefine_ptr;
        status.GetNRefine(&nrefine_ptr);
        printer << "\tn refine: " << nrefine_ptr << std::endl;

    }
    {
        int ntpoints_ptr;
        status.GetNTPoints(&ntpoints_ptr);
        printer << "\tn timepoints: " << ntpoints_ptr << std::endl;
    }
    {
        double estimate_ptr;
        status.GetSingleErrorEstAccess(&estimate_ptr);
        printer << "\terror estimate: " << estimate_ptr << std::endl;
    }
    {
        int callingfcn_ptr;
        status.GetCallingFunction(&callingfcn_ptr);
        printer << "\tcalling function: " << callingfcn_ptr << std::endl;
    }
    printer << " ]] ===== Braid Access Status ==== " << std::endl;
}

void print_status(std::ofstream & printer,BraidSyncStatus &status) {
    printer << " ===== Braid Sync Status ==== [[" << std::endl;
    {
        int i_upper;
        int i_lower;
        int level = 0;
        status.GetTIUL(&i_upper, &i_lower, level);
        printer << "\ti lower " << i_lower << std::endl;
        printer << "\ti upper: " << i_upper << std::endl;
    }
    {
        int nlevels_ptr;
        status.GetNLevels(&nlevels_ptr);
        printer << "\tn_level: " << nlevels_ptr << std::endl;
    }
    {
        int iter_ptr;
        status.GetIter(&iter_ptr);
        printer << "\titer: " << iter_ptr << std::endl;
    }
    {
        int level_ptr;
        status.GetLevel(&level_ptr);
        printer << "\tlevel: " << level_ptr << std::endl;
    }
    {
        int nrefine_ptr;
        status.GetNRefine(&nrefine_ptr);
        printer << "\tnrefine : " << nrefine_ptr << std::endl;
    }
    {
        int ntpoints_pt;
        status.GetNTPoints(&ntpoints_pt);
        printer << "\tntpoints: " << ntpoints_pt << std::endl;
    }
    {
        int done_ptr;
        status.GetDone(&done_ptr);
        printer << "\tdone: " << done_ptr << std::endl;
    }
    {
        int npoints_ptr;
        status.GetNumErrorEst(&npoints_ptr);
        printer << "\tn_points error est" << npoints_ptr << std::endl;
    }
    {
        double error_est_ptr;
        status.GetAllErrorEst(&error_est_ptr);
        printer << "\terr est " << error_est_ptr << std::endl;
    }

    {
        MPI_Comm *temporal_comm;
        status.GetTComm(temporal_comm);
    }
    {
        int proc_ptr;
        int level = 3;
        int index = 2;
        status.GetProc(&proc_ptr, level, index);
        printer << "\tproc(l=" << level << ", idx=" << index << "): " << proc_ptr << std::endl;
    }
    {
        int callingfcn_ptr;
        status.GetCallingFunction(&callingfcn_ptr);
        printer << "\tcalling fnc: " << callingfcn_ptr << std::endl;
    }
    printer << " ===== Braid Sync Status ==== " << std::endl;
}

void print_status(std::ofstream & printer,BraidStepStatus &status) {
    printer << " ===== Braid Step Status ==== [[" << std::endl;
    {
        int nrequest_ptr = 100;
        auto *rnorms = new double[nrequest_ptr];

        status.GetRNorms(&nrequest_ptr, rnorms);

        for (int i = 0; i < nrequest_ptr; i++) {
            printer << "\trnorm_" << i << ": " << rnorms[i] << std::endl;
        }
        delete[] rnorms;
    }

    {
        double tstart_ptr;
        double tstop_ptr;

        status.GetTstartTstop(&tstart_ptr, &tstop_ptr);
        printer << "\tt_start: " << tstart_ptr << std::endl;
        printer << "\tt_stop: " << tstop_ptr << std::endl;
    }

    /*{
        double tstart_ptr;
        status.GetT(&tstart_ptr);
        printer << "\tt_start: " << tstart_ptr << std::endl;
    }*/
    /*{
        double tstop_ptr;
        status.GetTstop(&tstop_ptr);
        printer << "\tt_stop: " << tstop_ptr << std::endl;
    }*/
    {
        int done;
        status.GetDone(&done);
        printer << "\tdone: " << done << std::endl;
    }
    {
        int index_ptr;
        status.GetTIndex(&index_ptr);
        printer << "\tt_idx: " << index_ptr << std::endl;
    }
    {
        int level_ptr;
        status.GetLevel(&level_ptr);
        printer << "\tlevel: " << level_ptr << std::endl;
    }
    {
        int nlevels_ptr;
        status.GetNLevels(&nlevels_ptr);
        printer << "\tn_level: " << nlevels_ptr << std::endl;
    }
    {
        int nrefine_ptr;
        status.GetNRefine(&nrefine_ptr);
        printer << "\tn_refine: " << nrefine_ptr << std::endl;
    }
    {
        int ntpoints_ptr;
        status.GetNTPoints(&ntpoints_ptr);
        printer << "\tn_timepoints: " << ntpoints_ptr << std::endl;
    }
    { // set
        // int rfactor = 2;
        // status.SetRFactor(rfactor);
    }
    {// set
        // int rspace = 2;
        // status.SetRSpace(rspace);
    }
    {
        double tol_ptr;
        status.GetTol(&tol_ptr);
        printer << "\ttol: " << tol_ptr << std::endl;
    }
    {
        int iter_ptr;
        status.GetIter(&iter_ptr);
        printer << "\titer: " << iter_ptr << std::endl;
    }
    {
        double old_fine_tolx_ptr;
        status.GetOldFineTolx(&old_fine_tolx_ptr);
        printer << "\toftolx: " << old_fine_tolx_ptr << std::endl;
    }
    {
        // double old_fine_tolx = 1.0;
        // status.SetOldFineTolx(old_fine_tolx);
    }
    {
        // double tight_fine_tolx = 1.0;
        // status.SetTightFineTolx(tight_fine_tolx);
    }
    {
        double estimate_ptr;
        status.GetSingleErrorEstStep(&estimate_ptr);
        printer << "\terror estimate: "<< estimate_ptr <<std::endl;
    }
    {
        double loose_tol = 1e-3;
        double tight_tol = 1e-9;
        double tol_ptr;
        status.GetSpatialAccuracy(loose_tol, tight_tol, &tol_ptr);
        printer << "\tadapt tol: " << tol_ptr << std::endl;
    }
    printer << " ]]===== Braid Step Status ==== " << std::endl;
}

void print_status(std::ofstream & printer,BraidCoarsenRefStatus &status) {
    printer << " ===== Braid Coarsen Ref Status ==== [[" << std::endl;
    {
        double tstart_ptr;
        double f_tprior_ptr;
        double f_tstop_ptr;
        double c_tprior_ptr;
        double c_tstop_ptr;
        status.GetTpriorTstop(&tstart_ptr, &f_tprior_ptr, &f_tstop_ptr, &c_tprior_ptr, &c_tstop_ptr);
        printer << "\ttstart_ptr: " << tstart_ptr<< std::endl;
        printer << "\tf_tprior_ptr: " << f_tprior_ptr<< std::endl;
        printer << "\tf_tstop_ptr: " << f_tstop_ptr<< std::endl;
        printer << "\tc_tprior_ptr: " << c_tprior_ptr<< std::endl;
        printer << "\tc_tstop_ptr: " << c_tstop_ptr<< std::endl;
    }
    {
        double tstart_ptr;
        status.GetT(&tstart_ptr);
        printer << "\ttstart_ptr: " << tstart_ptr<< std::endl;
    }
    {
        int tindex_ptr;
        status.GetTIndex(&tindex_ptr);
        printer << "\ttindex_ptr: " << tindex_ptr<< std::endl;
    }
    {
        int iter_ptr;
        status.GetIter(&iter_ptr);
        printer << "\titer_ptr: " << iter_ptr<< std::endl;
    }
    {
        double f_tstop_ptr;
        status.GetFTstop(&f_tstop_ptr);
        printer << "\tf_tstop_ptr: " << f_tstop_ptr<< std::endl;
    }
    {
        double f_tprior_ptr;
        status.GetFTprior(&f_tprior_ptr);
        printer << "\tf_tprior_ptr: " << f_tprior_ptr<< std::endl;
    }
    {
        double c_tstop_ptr;
        status.GetCTstop(&c_tstop_ptr);
        printer << "\tc_tstop_ptr: " << c_tstop_ptr<< std::endl;
    }
    {
        double c_tprior_ptr;
        status.GetCTprior(&c_tprior_ptr);
        printer << "\tc_tprior_ptr: " << c_tprior_ptr<< std::endl;
    }
    {
        int level_ptr;
        status.GetLevel(&level_ptr);
        printer << "\tlevel_ptr: " << level_ptr<< std::endl;
    }
    {
        int nlevels_ptr;
        status.GetNLevels(&nlevels_ptr);
        printer << "\tnlevels_ptr: " << nlevels_ptr<< std::endl;
    }
    {
        int nrefine_ptr;
        status.GetNRefine(&nrefine_ptr);
        printer << "\tnrefine_ptr: " << nrefine_ptr<< std::endl;
    }
    {
        int ntpoints_ptr;
        status.GetNTPoints(&ntpoints_ptr);
        printer << "\tntpoints_ptr: " << ntpoints_ptr<< std::endl;
    }
    printer << " ]] ===== Braid Coarsen Ref Status ==== " << std::endl;
}

void print_status(std::ofstream & printer,BraidBufferStatus &status) {
    printer << " ===== Braid Buffer Status ==== [[" << std::endl;
    {
        int messagetype_ptr;
        status.GetMessageType(&messagetype_ptr);
        printer << "\tmessage type: " << messagetype_ptr << std::endl;
    }
    {
        //int size = 20;
        //status.SetSize(size);
    }
    printer << " ]] ===== Braid Buffer Status ==== " << std::endl;
}

void print_status(std::ofstream & printer,BraidObjectiveStatus &status) {
    printer << " ===== Braid Objective Status ==== [[" << std::endl;
    {
        double tstart_ptr;
        status.GetT(&tstart_ptr);
        printer << "\ttstart_ptr : " << tstart_ptr << std::endl;
    }
    {
        int tindex_ptr;
        status.GetTIndex(&tindex_ptr);
        printer << "\ttindex_ptr : " << tindex_ptr << std::endl;
    }
    {
        int iter_ptr;
        status.GetIter(&iter_ptr);
        printer << "\titer_ptr : " << iter_ptr << std::endl;
    }
    {
        int level_ptr;
        status.GetLevel(&level_ptr);
        printer << "\tlevel_ptr : " << level_ptr << std::endl;
    }
    {
        int nlevels_ptr;
        status.GetNLevels(&nlevels_ptr);
        printer << "\tnlevels_ptr : " << nlevels_ptr << std::endl;
    }
    {
        int nrefine_ptr;
        status.GetNRefine(&nrefine_ptr);
        printer << "\tnrefine_ptr : " << nrefine_ptr << std::endl;
    }
    {
        int ntpoints_ptr;
        status.GetNTPoints(&ntpoints_ptr);
        printer << "\tntpoints_ptr : " << ntpoints_ptr << std::endl;
    }
    {
        double tol_ptr;
        status.GetTol(&tol_ptr);
        printer << "\ttol_ptr: " << tol_ptr << std::endl;
    }
    printer << " ]] ===== Braid Objective Status ==== " << std::endl;
}



#endif //UG_PLUGIN_XBRAIDFORUG4_BRAIDVECTORSTRUCT_H
