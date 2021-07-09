//
// Created by parnet on 05.06.21.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_STATUS_H
#define UG_PLUGIN_XBRAIDFORUG4_STATUS_H


/*
void print_status(BraidAccessStatus &status) {
    {
        double t_ptr;
        int iter_ptr;
        int level_ptr;
        int done_ptr;
        status.GetTILD(&t_ptr,
                       &iter_ptr,
                       &level_ptr,
                       &done_ptr);
    }
    {
        double t_ptr;
        status.GetT(&t_ptr);
    }
    {
        int tindex_ptr;
        status.GetTIndex(&tindex_ptr);
    }
    {
        int done_ptr;
        status.GetDone(&done_ptr);
    }
    {
        int level_ptr;
        status.GetLevel(&level_ptr);
    }
    {
        int nlevels_ptr;
        status.GetNLevels(&nlevels_ptr);
    }
    {
        int iter_ptr;
        status.GetIter(&iter_ptr);
    }
    {
        int wtest_ptr;
        status.GetWrapperTest(&wtest_ptr);
    }
    {
        double rnorm_ptr;
        status.GetResidual(&rnorm_ptr);
    }
    {
        int nrefine_ptr;
        status.GetNRefine(&nrefine_ptr);
    }
    {
        int ntpoints_ptr;
        status.GetNTPoints(&ntpoints_ptr);
    }
    {
        double estimate_ptr;
        status.GetSingleErrorEstAccess(&estimate_ptr);
    }
    {
        int callingfcn_ptr;
        status.GetCallingFunction(&callingfcn_ptr);
    }
}

void print_status(BraidSyncStatus &status) {
    {
        int i_upper;
        int i_lower;
        int level = 0;
        status.GetTIUL(&i_upper, &i_lower, level);
    }
    {
        int nlevels_ptr;
        status.GetNLevels(&nlevels_ptr);
    }
    {
        int iter_ptr;
        status.GetIter(&iter_ptr);
    }
    {
        int level_ptr;
        status.GetLevel(&level_ptr);
    }
    {
        int nrefine_ptr;
        status.GetNRefine(&nrefine_ptr);
    }
    {
        int ntpoints_pt;
        status.GetNTPoints(&ntpoints_pt);
    }
    {
        int done_ptr;
        status.GetDone(&done_ptr);
    }
    {
        int npoints_ptr;
        status.GetNumErrorEst(&npoints_ptr);
    }
    {
        double error_est_ptr;
        status.GetAllErrorEst(&error_est_ptr);
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
    }
    {
        int callingfcn_ptr;
        status.GetCallingFunction(&callingfcn_ptr);
    }
}

void print_status(BraidStepStatus &status) {
    {
        int nrequest_ptr = 2;
        auto *rnorms = new double[nrequest_ptr];

        status.GetRNorms(&nrequest_ptr, rnorms);

        delete[] rnorms;
    }

    {
        double tstart_ptr;
        double tstop_ptr;

        status.GetTstartTstop(&tstart_ptr, &tstop_ptr);
    }

    {
        double tstart_ptr;
        status.GetT(&tstart_ptr);
    }
    {
        double tstop_ptr;
        status.GetTstop(&tstop_ptr);
    }
    {
        int done;
        status.GetDone(&done);
    }
    {
        int index_ptr;
        status.GetTIndex(&index_ptr);
    }
    {
        int level_ptr;
        status.GetLevel(&level_ptr);
    }
    {
        int nlevels_ptr;
        status.GetNLevels(&nlevels_ptr);
    }
    {
        int nrefine_ptr;
        status.GetNRefine(&nrefine_ptr);
    }
    {
        int ntpoints_ptr;
        status.GetNTPoints(&ntpoints_ptr);
    }
    { // set
        int rfactor = 2;
        status.SetRFactor(rfactor);
    }
    {// set
        int rspace = 2;
        status.SetRSpace(rspace);
    }
    {
        double tol_ptr;
        status.GetTol(&tol_ptr);
    }
    {
        int iter_ptr;
        status.GetIter(&iter_ptr);
    }
    {
        double old_fine_tolx_ptr;
        status.GetOldFineTolx(&old_fine_tolx_ptr);
    }
    {
        double old_fine_tolx = 1.0;
        status.SetOldFineTolx(old_fine_tolx);
    }
    {
        double tight_fine_tolx = 1.0;
        status.SetTightFineTolx(tight_fine_tolx);
    }
    {
        double estimate_ptr;
        status.GetSingleErrorEstStep(&estimate_ptr);
    }
    {
        double loose_tol = 1e-3;
        double tight_tol = 1e-9;
        double tol_ptr;
        status.GetSpatialAccuracy(loose_tol, tight_tol, &tol_ptr);
    }
}

void print_status(BraidCoarsenRefStatus &status) {
    {
        double tstart_ptr;
        double f_tprior_ptr;
        double f_tstop_ptr;
        double c_tprior_ptr;
        double c_tstop_ptr;
        status.GetTpriorTstop(&tstart_ptr, &f_tprior_ptr, &f_tstop_ptr, &c_tprior_ptr, &c_tstop_ptr);
    }
    {
        double tstart_ptr;
        status.GetT(&tstart_ptr);
    }
    {
        int tindex_ptr;
        status.GetTIndex(&tindex_ptr);
    }
    {
        int iter_ptr;
        status.GetIter(&iter_ptr);
    }
    {
        double f_tstop_ptr;
        status.GetFTstop(&f_tstop_ptr);
    }
    {
        double f_tprior_ptr;
        status.GetFTprior(&f_tprior_ptr);
    }
    {
        double c_tstop_ptr;
        status.GetCTstop(&c_tstop_ptr);
    }
    {
        double c_tprior_ptr;
        status.GetCTprior(&c_tprior_ptr);
    }
    {
        int level_ptr;
        status.GetLevel(&level_ptr);
    }
    {
        int nlevels_ptr;
        status.GetNLevels(&nlevels_ptr);
    }
    {
        int nrefine_ptr;
        status.GetNRefine(&nrefine_ptr);
    }
    {
        int ntpoints_ptr;
        status.GetNTPoints(&ntpoints_ptr);
    }

}

void print_status(BraidBufferStatus &status) {
    {
        int messagetype_ptr;
        status.GetMessageType(&messagetype_ptr);
    }
    {
        int size = 20;
        status.SetSize(size);
    }
}

void print_status(BraidObjectiveStatus &status) {
    {
        double tstart_ptr;
        status.GetT(&tstart_ptr);
    }
    {
        int tindex_ptr;
        status.GetTIndex(&tindex_ptr);
    }
    {
        int iter_ptr;
        status.GetIter(&iter_ptr);
    }
    {
        int level_ptr;
        status.GetLevel(&level_ptr);
    }
    {
        int nlevels_ptr;
        status.GetNLevels(&nlevels_ptr);
    }
    {
        int nrefine_ptr;
        status.GetNRefine(&nrefine_ptr);
    }
    {
        int ntpoints_ptr;
        status.GetNTPoints(&ntpoints_ptr);
    }
    {
        double tol_ptr;
        status.GetTol(&tol_ptr);
    }
} */


#endif //UG_PLUGIN_XBRAIDFORUG4_STATUS_H
