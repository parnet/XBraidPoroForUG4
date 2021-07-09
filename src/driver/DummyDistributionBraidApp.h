//
// Created by parnet on 19.05.21.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_DUMMYDISTRIBUTIONBRAIDAPP_H
#define UG_PLUGIN_XBRAIDFORUG4_DUMMYDISTRIBUTIONBRAIDAPP_H
/*
#include "BraidVectorStruct.h"
#include "SpaceTimeCommunicator.h"
#include "MemoryObserver.h"
#include "BraidTimer.h"
#include "BraidUsageTimer.h"
#include <iomanip>
#include <sstream>

#include "lib_disc/io/vtkoutput.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/function_spaces/interpolate.h"
#include "lib_disc/spatial_disc/domain_disc.h"
#include "lib_algebra/vector_interface/vec_functions.h"

class DummyDistributionBraidApp : public BraidApp {
     public: // todo set better modes
    const char *name;
    SmartPtr<SpaceTimeCommunicator> m_comm;
    std::ofstream o;

    DummyDistributionBraidApp(MPI_Comm _comm_t, braid_Real _tstart = 0.0, braid_Real _tstop = 6.0,
                              braid_Int _ntime = 6)
            : BraidApp(_comm_t, _tstart, _tstop, _ntime) {
        std::cout << "Braid App constructor" << std::endl;
    }

    virtual ~DummyDistributionBraidApp() {
        std::cout << "Braid App destroyed" << std::endl;
    }



    virtual braid_Int Step(braid_Vector u_,
                           braid_Vector ustop_,
                           braid_Vector fstop_,
                           BraidStepStatus &pstatus) override {
        std::cout << "Braid App step" << std::endl;
        return 0;
    };

    virtual braid_Int Residual(braid_Vector u_,
                               braid_Vector r_,
                               BraidStepStatus &pstatus) override {
        std::cout << "Braid App residual" << std::endl;
        return 0;
    }

    virtual braid_Int Clone(braid_Vector u_,
                            braid_Vector *v_ptr) override {
        std::cout << "Braid App clone" << std::endl;
        return 0;
    }

    virtual braid_Int Init(braid_Real t, braid_Vector *u_ptr) override {
        std::cout << "Braid App init" << std::endl;
        return 0;
    }

    virtual braid_Int Free(braid_Vector u_) override {
        std::cout << "Braid App free" << std::endl;
        return 0;
    };

    virtual braid_Int Sum(braid_Real alpha,
                          braid_Vector x_,
                          braid_Real beta,
                          braid_Vector y_) override {
        std::cout << "Braid App sum" << std::endl;
        return 0;
    };

    virtual braid_Int SpatialNorm(braid_Vector u_,
                                  braid_Real *norm_ptr) override {
        std::cout << "Braid App norm" << std::endl;
        return 0;
    };

    /// @see braid_PtFcnAccess.
    virtual braid_Int Access(braid_Vector u_,
                             BraidAccessStatus &astatus) override {
        std::cout << "Braid App access" << std::endl;
        return 0;
    };

    /// @see braid_PtFcnBufSize.
    virtual braid_Int BufSize(braid_Int *size_ptr,
                              BraidBufferStatus &bstatus) override {
        std::cout << "Braid App buffer size" << std::endl;
        return 0;
    };

    /// @see braid_PtFcnBufPack.
    virtual braid_Int BufPack(braid_Vector u_,
                              void *buffer,
                              BraidBufferStatus &bstatus) override {
        std::cout << "Braid App buffer pack" << std::endl;
        return 0;
    };

    /// @see braid_PtFcnBufUnpack.
    virtual braid_Int BufUnpack(void *buffer,
                                braid_Vector *u_ptr,
                                BraidBufferStatus &bstatus) override {
        std::cout << "Braid App buffer unpack" << std::endl;
        return 0;
    };

    virtual braid_Int Coarsen(braid_Vector fu_,
                              braid_Vector *cu_ptr,
                              BraidCoarsenRefStatus &status) override {
        std::cout << "Braid App coarsen" << std::endl;
        Clone(fu_, cu_ptr);
        return 0;
    }

    /// @see braid_PtFcnSRefine.
    virtual braid_Int Refine(braid_Vector cu_,
                             braid_Vector *fu_ptr,
                             BraidCoarsenRefStatus &status) override {
        std::cout << "Braid App refine" << std::endl;
        Clone(cu_, fu_ptr);
        return 0;
    }

    virtual braid_Int Sync(BraidSyncStatus &sstatus) override {
        std::cout << "Braid App sync" << std::endl;
        return 0;
    }
};


#endif //UG_PLUGIN_XBRAIDFORUG4_DUMMYDISTRIBUTIONBRAIDAPP_H
*/
#endif