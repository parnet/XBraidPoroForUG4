// zzzzzz
// Created by parnet on 28.05.21.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_BRAIDCONV_H
#define UG_PLUGIN_XBRAIDFORUG4_BRAIDCONV_H

#include <iostream>
#include <lib_algebra/operator/convergence_check.h>

#include "../util/BraidTimer.h"
#include "../util/BraidUsageTimer.h"
#include "../core/BraidVectorStruct.h"
/*
template <typename TAlgebra>
class IBraidConvControl{
    typedef ug::IConvergenceCheck<typename TAlgebra::vector_type> TConv;
    typedef SmartPtr<TConv> SPConv;

    virtual SPConv getAccuracy(int time_grid_level, int iteration) = 0;

    virtual void  print_settings() = 0;
};

template <typename TAlgebra>
class BraidStaticConv : public IBraidConvControl<TAlgebra> {

    typedef ug::IConvergenceCheck<typename TAlgebra::vector_type> TConv;
    typedef SmartPtr<TConv> SPConv;

    double relReduction = 1e-3;
    int maxSteps = 100;
    double minDefect = 1-9;

    SPConv getAccuracy(int time_grid_level, int iteration) override{
        return StdConvCheck(maxSteps, minDefect, relReduction);
    };

    void print_settings(){
        std::cout << "===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====" << std::endl;
        std::cout << "===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====" << std::endl;
    }
};

template <typename TAlgebra>
class BraidStaticLeveldependend : public IBraidConvControl<TAlgebra> {
    typedef ug::IConvergenceCheck<typename TAlgebra::vector_type> TConv;
    typedef SmartPtr<TConv> SPConv;

    double tol = 1e-3;
    SPConv getAccuracy(int time_grid_level, int iteration) override{

    };
};
/*
template <typename TAlgebra>
class BraidConvControl {
    typedef ug::IConvergenceCheck<typename TAlgebra::vector_type> TConv;
    typedef SmartPtr<TConv> SPConv;

public: // todo set better modes
    number m_loose_tol = 1e-2;
    number m_tight_tol = 1e-9;

    double norm_request[2];

    bool m_exact_first_iteration = false;
    bool m_exact_recurring_iteration = false;
    size_t m_recurring_interval = 4;

    int lastiter[15] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

    SPConv getAccuracy(int time_grid_level) {

        //StartRedoranLevel(LevelObserver::T_ADAPTIVE_TOL, l);
        ///int iteration;
        //pstatus.GetIter(&iteration);

        //int iter;
        //double tol;
        //if (iteration != this->lastiter[l]) {
        //    this->lastiter[l] = iteration;
        //    if (this->m_timing) {
        //        double diff, total;
        //        this->timer.now(total, diff);
        //        this->debugwriter << std::setw(10) << "@time:"
        //                          << std::setw(12) << total << " ; "
        //                          << std::setw(12) << diff << " Begin iteration for level" << l << std::endl;
        //    }
            //if (this->m_verbose) {
            //    this->o << "========== ==========<< " << l << " >>========== ==========" << std::endl;
            //}
        //    if(this->strongFirstIteration && iteration == 0) {
        //        tol = tight_tol;
        //        iter = 10;
        //        this->debugwriter << "strong first cycle    ";
        //    } else if (this->strongRecurring && iteration % this->recurringInterval == 0){
        //        tol = tight_tol;
        //        iter = 10;
        //        this->debugwriter << "strong recurring    ";
        //    } else if (l == 0) {
        //        this->GetSpatialAccuracy(pstatus, l, loose_tol, tight_tol, &tol);
        //        iter = 100;
        //        this->debugwriter << "level 0     ";
        //    } else {
        //        tol = this->loose_tol;
        //        iter = 2;
        //        this->debugwriter << "level != 0    ";
        //    }

        //    SPConv conv = SPConv(new TConv(iter, 1e-32, tol, false));
        //    this->debugwriter << "% new Tolerance \t " << tol << std::endl;
        //    //this->m_linSolver->set_convergence_check(conv);
        //}
        //StopRedoranLevel(LevelObserver::T_ADAPTIVE_TOL, l);
    }

    braid_Int GetSpatialAccuracy(BraidStepStatus &sstatus,
                                 int level,
                                 double loose_tol,
                                 double tight_tol,
                                 double *tol_ptr) {

        //X void GetRNorms(braid_Int *nrequest_ptr, braid_Real *norm_request) // residual history, negative -> last n norms, positive -> first  n norms
        void GetTstartTstop(braid_Real *tstart_ptr, braid_Real *tstop_ptr)
        void GetT(braid_Real *tstart_ptr)
        void GetTstop(braid_Real *tstop_ptr)
        void GetDone(braid_Int *done)
        void GetTIndex(braid_Int *tindex_ptr)
        void GetLevel(braid_Int *level_ptr)
        void GetNLevels(braid_Int *nlevels_ptr)
        void GetNRefine(braid_Int *nrefine_ptr)
        void GetNTPoints(braid_Int *ntpoints_ptr)
        void SetRFactor(braid_Int rfactor)
        void SetRSpace(braid_Int rspace)
        //X void GetTol(braid_Real *tol_ptr)
        void GetIter(braid_Int *iter_ptr)
        void GetOldFineTolx(braid_Real *old_fine_tolx_ptr)
        // X void SetOldFineTolx(braid_Real old_fine_tolx)
        // X void SetTightFineTolx(braid_Int tight_fine_tolx)
        void GetSingleErrorEstStep(braid_Real *estimate_ptr)
        void GetSpatialAccuracy( braid_Real loose_tol, braid_Real tight_tol, braid_Real *tol_ptr)/

        /*
        braid_Real stol, tol, rnorm, rnorm0, old_fine_tolx;
        braid_Real l_rnorm, l_ltol, l_ttol, l_tol;

        sstatus.GetTol(&tol);
        sstatus.GetLevel(&level)
        //xbraid_updated: pstatus.StepStatusGetTol(&tol);
        sstatus.GetOldFineTolx(old_fine_tolx);


        /* Get the first and then the current residual norms /
        int nrequest = 2; // positive -> first  n norms
        norm_request[0] = -1.0;
        norm_request[1] = -1.0;
        sstatus.GetRNorms(&nrequest, norm_request);
        if ((norm_request[0] == -1.0) && (norm_request[1] != -1.0)) {
            rnorm0 = norm_request[1];
        } else {
            rnorm0 = norm_request[0];
        }


        nrequest = -2; // negative -> last n norms
        sstatus.GetRNorms(&nrequest, norm_request);
        if ((norm_request[1] == -1.0) && (norm_request[0] != -1.0)) {
            rnorm = norm_request[0];
        } else {
            rnorm = norm_request[1];
        }


        if ((level > 0) || (nrequest == 0) || (rnorm0 == -1.0)) {
            /* Always return the loose tolerance, if
             * (1) On a coarse grid computation
             * (2) There is no residual history yet (this is the first Braid iteration with skip turned on) /
            *tol_ptr = loose_tol;
        } else {
            /* Else, do a variable tolerance for the fine grid /
            l_rnorm = -log10(rnorm / rnorm0);
            l_tol = -log10(tol / rnorm0);
            l_ltol = -log10(loose_tol);
            l_ttol = -log10(tight_tol);

            if (l_rnorm >= (7.0 / 8.0) * l_tol) {
                /* Close to convergence, return tight_tol /
                *tol_ptr = tight_tol;
            } else {
                /* linear interpolation between loose_tol and tight_tol/
                stol = (l_rnorm / l_tol) * (l_ttol - l_ltol) + l_ltol;
                *tol_ptr = pow(10, -stol);
                  if (((*tol_ptr) > old_fine_tolx) && (old_fine_tolx > 0)) {
                    *tol_ptr = old_fine_tolx; // take old tolerance for fine grid if calculated value is more loose.
                }
            }
        }

        if (level == 0) {
            /* Store this fine grid tolerance /
            sstatus.SetOldFineTolx(*tol_ptr);
            std::cout << "new tolerance " << *tol_ptr << std::endl;
            /* If we've reached the "tight tolerance", then indicate to Braid that we can halt /
            if (*tol_ptr == tight_tol) {
                std::cout << "tolerance reached" << std::endl;
                sstatus.SetTightFineTolx(1);
            } else {
                sstatus.SetTightFineTolx(0);
            }

        }

        return 0;
    }

    void print_summary(){
        if (this->m_adapt_conv_check) {
            std::cout << "\t loose: " << this->m_loose_tol << std::endl;
            std::cout << "\t tight: " << this->m_tight_tol << std::endl;
            std::cout << "\t exact first iteration: " << this->m_exact_first_iteration << std::endl;
            std::cout << "\t exact recurring iteration: " << this->m_exact_recurring_iteration << std::endl;
            if (this->m_exact_recurring_iteration) {
                std::cout << "\t \tinterval: " << this->recurringInterval << std::endl;
            }
        }
    }

    void setStrongFirstIteration(bool b_state) {
        this->m_exact_first_iteration = b_state;
    }


    void setRecurringStrongIteration(int interval) {
        this->m_exact_recurring_iteration = (interval != 0);
        this->strongFirstIteration = interval;
    }


    void setLooseTol(double pLooseTol) {
        this->m_loose_tol = pLooseTol;
    }

    void setTightTol(double pTightTol) {
        this->m_tight_tol = pTightTol;
    }
    void print_settings(){
        std::cout << "===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====" << std::endl;
        std::cout << "===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====" << std::endl;
    }
};
*/

#endif //UG_PLUGIN_XBRAIDFORUG4_BRAIDCONV_H
