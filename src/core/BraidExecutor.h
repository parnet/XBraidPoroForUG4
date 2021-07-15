// zzzzz
// Created by parnet on 12.05.19.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_BRAIDEXECUTOR_H
#define UG_PLUGIN_XBRAIDFORUG4_BRAIDEXECUTOR_H

#include "common/assert.h"
// #include "pcl/pcl_comm_world.h"

#include "bindings/lua/lua_function_handle.h"
#include "bindings/lua/lua_user_data.h"

// #include "lib_disc/time_disc/theta_time_step.h"
// #include "lib_disc/time_disc/solution_time_series.h"

// #include "PBraidApp.h"

#include "SpaceTimeCommunicator.h"
#include "../driver/BraidIntegrator.h"
#include "../driver/BraidTimeStepper.h"

#include "../util/Scriptor.h"

class BraidSettings {
public:
    // https://github.com/XBraid/xbraid/blob/master/braid/braid.c
    bool m_increase_max_level = false;
    bool m_skip = false;
    bool m_cycle_fmg = false;
    bool m_residual = false;
    bool m_sync = false;
    bool m_final_fc_relax = false;
    bool finished = false;
    bool m_file_io_level = true;
    bool m_relax_only_cg = false;
    bool m_est_error = false; // richardson-based error estimation
    bool m_richardson = false; //  richardson-based extrapolation for finest grid
    bool m_sequential = false;
    bool m_full_r_norm = false;
    bool m_coarsen_and_refine = false;
    bool m_refine = false;

    int eval_level = 0;
    int m_max_level = 30;
    int m_min_coarse = -1;
    int m_max_iter = 100;
    int m_temp_norm = 2;
    int m_store_values = -1;
    int m_t_points_cutoff = -1;
    int m_max_refinements = 200;
    int m_access_level = 1;
    int m_print_level = 2;
    int m_cycle_nfmg = -1 ;
    int m_cycle_nfmgv = 1;
    int m_reverted_ranks = 0;
    int m_local_order = -1; // 2 for bachward euler
    int m_periodic = -1;




    double m_abs_tol = -1.0;
    double m_rel_tol = -1.0;

    const char * m_print_file = "braid_runtime.out";

    int m_n_relax_default = 1;
    std::map<int, int> m_n_relax;

    int m_c_factor_default = 2 ;
    std::map<int, int> m_c_factor;

    double m_c_relax_weight_default = 1;
    std::map<int, double> m_c_relax_weight;

};

template<typename TDomain, typename TAlgebra>
class BraidExecutor {
public: // todo set better modes
    typedef ug::ThetaTimeStep<TAlgebra> TTimeStep;
    //typedef ug::ILinearOperatorInverse<typename TAlgebra::vector_type> TSolver;
    typedef typename TAlgebra::vector_type vector_type;

    typedef ug::GridFunction<TDomain, TAlgebra> TGridFunction;
    typedef SmartPtr<TGridFunction> SPGridFunction;
    typedef ConstSmartPtr<TGridFunction> CSPGridFunction;


    typedef BraidGridFunctionBase<TDomain, TAlgebra> TBraidBaseApp;
    typedef SmartPtr<TBraidBaseApp> SPBraidBaseApp;

    typedef BraidIntegrator <TDomain, TAlgebra> TBraidIntegrator;
    typedef SmartPtr<TBraidIntegrator> SPBraidIntegrator;

    typedef BraidTimeStepper <TDomain, TAlgebra> TBraidTimeStepper;
    typedef SmartPtr<TBraidTimeStepper> SPBraidTimeStepper;

    typedef SmartPtr<ug::DomainDiscretization<TDomain, TAlgebra>> SPDomainDisc;

    typedef SmartPtr<SpaceTimeCommunicator> SPXCommunicator;
    typedef SmartPtr<Scriptor<TDomain, TAlgebra>> SPScriptor;

    typedef SmartPtr<BraidInitializer<TDomain,TAlgebra>> SPBraidInitializer;

    typedef BraidSettings TBraidSettings;
    typedef SmartPtr<TBraidSettings> SPBraidSettings;

    typedef Paralog TParalog;
    typedef SmartPtr<TParalog> SPParalog;

    typedef BraidSpatialNorm<TDomain,TAlgebra> TSpatialNorm;
    typedef SmartPtr<TSpatialNorm> SPSpatialNorm;

    // -----------------------------------------------------------------------------------------------------------------
    // member variables
    // -----------------------------------------------------------------------------------------------------------------
    SPParalog m_log;
    BraidCore *m_braid_core = nullptr;
    SPBraidBaseApp m_app;
    SmartPtr<SpaceTimeCommunicator> m_comm;
    const char *m_filename;
    BraidSettings m_braid_settings;

    // -----------------------------------------------------------------------------------------------------------------
    // constructor and destructor
    // -----------------------------------------------------------------------------------------------------------------
    BraidExecutor(SPXCommunicator& p_comm, SPBraidBaseApp& p_app) { // todo ctrlx
        this->m_comm = p_comm;
        this->set_app(p_app);
    }

    ~BraidExecutor() { // todo ctrlx
        delete m_braid_core;
    }




    // -----------------------------------------------------------------------------------------------------------------
    //
    // -----------------------------------------------------------------------------------------------------------------


    void set_residual(bool residual) { // todo ctrlx
        if (residual && this->m_app->provide_residual) {
            this->m_braid_settings.m_residual = residual;
            this->m_braid_core->SetResidual();
        }
    }

    void set_n_relax(int level, int number){
        if(level == -1){
            this->m_braid_settings.m_n_relax_default = number;
        } else {
            this->m_braid_settings.m_n_relax[level] = number;
            if (level >  this->m_braid_settings.eval_level){
                this->m_braid_settings.eval_level = level;
            }
        }
        this->m_braid_core->SetNRelax(level, number);
    }

    void set_c_factor(int level, int factor) {
        if(level == -1){
            this->m_braid_settings.m_c_factor_default = factor;
        } else {
            this->m_braid_settings.m_c_factor[level] = factor;
            if (level >  this->m_braid_settings.eval_level){
                this->m_braid_settings.eval_level = level;
            }
        }
        this->m_braid_core->SetCFactor(level, factor);
    }

    void set_max_levels(int maxLevel) {
        this->m_braid_settings.m_max_level = maxLevel;
        this->m_braid_core->SetMaxLevels(maxLevel);
        this->m_app->set_max_levels(maxLevel);
    }

    void set_skip_downcycle_work(bool skip){
        if(skip) {
            this->m_braid_settings.m_skip = skip;
            this->m_braid_core->SetSkip(skip);
        }
    }


    void set_min_coarse(int minCoarse) { // todo what about
        this->m_braid_settings.m_min_coarse = minCoarse;
        this->m_braid_core->SetMinCoarse(minCoarse);
    }


    // todo set ConvCheck


    void set_max_iterations(int max_iter) {
        this->m_braid_settings.m_max_iter = max_iter;
        this->m_braid_core->SetMaxIter(max_iter);
    }

    void set_absolute_tol(double tol) {
        this->m_braid_settings.m_abs_tol = tol;
        this->m_braid_core->SetAbsTol(tol);
    }

    void set_relative_tol(double tol) {
        this->m_braid_settings.m_rel_tol = tol;
        this->m_braid_core->SetRelTol(tol);
    }

    void set_temporal_norm(int nrm) {
        this->m_braid_settings.m_temp_norm = nrm;
        this->m_braid_core->SetTemporalNorm(nrm);
    }

    void set_sequential(bool sequential) {
        this->m_braid_settings.m_sequential = sequential;
        if (sequential) {
            this->m_braid_core->SetSeqSoln(1);
        } else {
            this->m_braid_core->SetSeqSoln(0);
        }
    }

    void set_store_values(int level) {
        this->m_braid_settings.m_store_values = level;
        this->m_braid_core->SetStorage(level);
    }

    void set_spatial_coarsen_and_refine(bool cnr) {
        if (cnr) {
            this->m_braid_settings.m_coarsen_and_refine = cnr;
            this->m_braid_core->SetSpatialCoarsenAndRefine();
        }
    }

    void set_refine(bool ref) { // todo ctrlx
        this->m_braid_settings.m_refine = ref;
        if(ref) {
            this->m_braid_core->SetRefine(1);
        } else {
            this->m_braid_core->SetRefine(0);
        }
    }

    void set_max_refinements(int number) { // todo ctrlx
        this->m_braid_settings.m_max_refinements = number;
        this->m_braid_core->SetMaxRefinements(number);
    }

    void set_access_level(int level) {
        this->m_braid_settings.m_access_level = level;
        this->m_braid_core->SetAccessLevel(level);
    }


    void set_print_level(int level) {
        this->m_braid_settings.m_print_level = level;
        this->m_braid_core->SetPrintLevel(level);
    }

    void set_default_print_file(){
        const char *file = "braid_runtime.out";
        this->set_print_file(file);
    }

    void set_print_file(const char *file) { // todo ctrlx
        this->m_braid_settings.m_print_file = file;
        this->m_braid_core->SetPrintFile(file);
    }

    void set_cycle_fmg() {
        this->m_braid_settings.m_cycle_fmg = true;
        this->m_braid_core->SetFMG();
    }

    void set_cycle_nfmg(int vu) {
        this->m_braid_settings.m_cycle_nfmg = vu;
        this->m_braid_core->SetNFMG(vu);
    }

    void set_cycle_nfmgv(int mu) {
        this->m_braid_settings.m_cycle_nfmgv = mu;
        this->m_braid_core->SetNFMGVcyc(mu);
    }

    void set_cycle_type(const char *ctyp) { // todo ctrlx
        std::cout << "debug::BraidExecutor::setCycleType(args)" << std::endl<<std::flush;
        if (strcmp(ctyp, "V FCF") == 0) {
            this->m_braid_core->SetNRelax(-1, 1);
            this->m_braid_core->SetNRelax(0, 1);

        } else if (strcmp(ctyp, "V F") == 0) {
            this->m_braid_core->SetNRelax(-1, 0);
            this->m_braid_core->SetNRelax(0, 0);

        } else if (strcmp(ctyp, "V F-FCF") == 0) {
            this->m_braid_core->SetNRelax(-1, 1);
            this->m_braid_core->SetNRelax(0, 0);

        } else if (strcmp(ctyp, "F FCF") == 0) {
            this->m_braid_core->SetNRelax(-1, 1);
            this->m_braid_core->SetNRelax(0, 1);
            this->m_braid_settings.m_cycle_fmg = true;
            this->m_braid_core->SetFMG();


        } else if (strcmp(ctyp, "F F") == 0) {
            this->m_braid_core->SetNRelax(-1, 0);
            this->m_braid_core->SetNRelax(0, 0);
            this->m_braid_settings.m_cycle_fmg = true;
            this->m_braid_core->SetFMG();


        } else if (strcmp(ctyp, "F F-FCF") == 0) {
            this->m_braid_core->SetNRelax(-1, 1);
            this->m_braid_core->SetNRelax(0, 0);
            this->m_braid_settings.m_cycle_fmg = true;
            this->m_braid_core->SetFMG();
        } else {
            std::cout << "invalid cycle type" << std::endl;
        }
    }

    void set_sync(){ // todo what about?
        this->m_braid_settings.m_sync = true;
        this->m_braid_core->SetSync();
    }

    void set_increase_max_levels(){ // todo what about?
        this->m_braid_settings.m_increase_max_level = true;
        this->m_braid_core->SetIncrMaxLevels();
    }

    void set_relax_only_cg(bool setting){ // todo what about?
        this->m_braid_settings.m_relax_only_cg = setting;
        if ( setting) {
            this->m_braid_core->SetRelaxOnlyCG(1);
        } else{
            this->m_braid_core->SetRelaxOnlyCG(0);
        }
    }

    void set_agg_c_factor(int cfactor0){
        // todo controll braid implementation
        this->m_braid_core->SetAggCFactor(cfactor0);
    }

    void set_periodic(int periodic){ // todo what about
        this->m_braid_settings.m_periodic = periodic;
        this->m_braid_core->SetPeriodic(periodic);
    }

    void set_final_fc_relax(){// todo what about?
        this->m_braid_settings.m_final_fc_relax = true;
        this->m_braid_core->SetFinalFCRelax();
    }

    void set_reverted_ranks(int ranks){ // todo what about? -> adjoint problems
        this->m_braid_settings.m_reverted_ranks = ranks;
        this->m_braid_core->SetRevertedRanks(ranks);
    }

    void set_richardson_estimation(bool use_richardson, bool use_extrapolation, int  local_order){
        this->m_braid_settings.m_est_error = use_richardson;
        this->m_braid_settings.m_richardson = use_extrapolation;
        this->m_braid_settings.m_local_order = local_order;
        this->m_braid_core->SetRichardsonEstimation(use_richardson, use_extrapolation, local_order);
    }

    void set_file_io_level(int level){
        this->m_braid_settings.m_file_io_level = level;
        this->m_braid_core->SetFileIOLevel(level);
    }

    void set_c_relax_weight(int level, double weight){
        //todo structure
        if(level == -1){
            this->m_braid_settings.m_c_relax_weight_default = weight;
        } else {
            this->m_braid_settings.m_c_relax_weight[level] = weight;
            if (level >  this->m_braid_settings.eval_level){
                this->m_braid_settings.eval_level = level;
            }
        }
        this->m_braid_core->SetCRelaxWt(level, weight);
    }

    void set_t_points_cutoff(int cutoff){ // todo what about
        this->m_braid_settings.m_t_points_cutoff = cutoff;
        this->m_braid_core->SetTPointsCutoff(cutoff);
    }


    void set_full_residual_norm(){
        this->m_braid_settings.m_full_r_norm = true;
        this->m_braid_core->SetFullRNormRes(_BraidAppResidual);
    }

    void set_time_grid(){ // todo what about
        // todo print summary
        // todo print settings
        //todo this->m_braid_core->SetTimeGrid(braid_PtFcnTimeGrid tgrid)
    }

    int get_num_iteration(){
        int iter = 0;
        this->m_braid_core->GetNumIter(&iter);
        return iter;
    }

    void get_c_factor(){
        // todo check braid implementation for level == 0
        // this->m_braid_core->GetCFactor(braid_Int *cfactor_ptr) // todo
    }

    void get_residual_norms(){ // todo ctrlx
        //todo
        //todo this->m_braid_core->GetRNorms(braid_Int *nrequest_ptr, braid_Real *rnorms)
    }



    int get_num_level(){
        int number_of_level = 0;
        this->m_braid_core->GetNLevels(&number_of_level);
        return number_of_level;
    }



    int get_warm_restart(){
       return this->m_braid_core->GetWarmRestart();
    }

    int get_distribution_lower(){
        int lower;
        int upper;
        this->m_braid_core->GetDistribution(&lower,&upper);
        return lower;
    }

    int get_distribution_upper(){
        int lower;
        int upper;
        this->m_braid_core->GetDistribution(&lower,&upper);
        return upper;
    }

    void get_distribution(int &lower, int &upper){ // todo ctrlx
        this->m_braid_core->GetDistribution(&lower,&upper);
    }

    int get_id(){ // todo ctrlx
        int id = 0;
        this->m_braid_core->GetMyID(&id);
        return id;
    }





    int run() { // todo ctrlx
        this->m_app->init();
        this->m_log->o << "debug::BraidExecutor::run()" << std::endl<<std::flush;
        this->m_braid_core->Drive();

        // GetRNorms
        // GetNLevels
        // GetNumIter
        return 0;
    }
    braid_Core get_real_braid_core(){ // todo ctrlx
        return this->m_braid_core->GetCore();
    }

    // -----------------------------------------------------------------------------------------------------------------
    // Braid Exectution Control
    // -----------------------------------------------------------------------------------------------------------------

    BraidCore * get_braid_core(){
        return this->m_braid_core;
    }




    void reset_core(SPBraidBaseApp p_app){ // todo ctrlx
        delete m_braid_core;
        m_braid_core = new BraidCore(this->m_comm->GLOBAL, this->m_app.get());
    }

    void set_app(SPBraidBaseApp p_app) { // todo ctrlx
        this->m_app = p_app;
        this->m_app->m_comm = this->m_comm;
        this->m_app->comm_t = this->m_comm->TEMPORAL;
        this->reset_core(p_app);
    }

    SPBraidBaseApp get_app() { // todo ctrlx
        return this->m_app;
    }

    void set_settings(SPBraidSettings settings) { // todo ctrlx
    }

    SPBraidSettings get_settings(){ // todo ctrlx
    }
    void apply_settings(TBraidSettings settings){ // todo ctrlx

    }

    void create_access(SmartPtr<ug::VTKOutput<TDomain::dim>> out) { // todo ctrlx
        // this->m_out = out;-
        // mgf->m_out = this->m_out;
    }

    void set_filename(const char *filename) {  // todo ctrlx
        // mgf->m_filename = filename;
        this->m_filename = filename;
    }


    void print_settings(){ // todo ctrlx
        // this->m_log->o << "debug::BraidExecutor::print_settings()" << std::endl<<std::flush;
        this->m_log->o << "===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====" << std::endl;
        this->m_log->o << "Module Name:                                             Braid Executor" << std::endl;
        this->m_log->o << "----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----" << std::endl;
        this->m_log->o << "Max Level: " <<  this->m_braid_settings.m_max_level << std::endl;
        this->m_log->o << "Increase max level: " <<  this->m_braid_settings.m_increase_max_level << std::endl;
        this->m_log->o << "Skip Down Cylce Work: " << this->m_braid_settings.m_skip << std::endl;

        this->m_log->o << "Absolute Tol: " << this->m_braid_settings.m_abs_tol << std::endl;
        this->m_log->o << "Relative Tol: " << this->m_braid_settings.m_rel_tol << std::endl;
        this->m_log->o << "Temporal Norm: " << this->m_braid_settings.m_temp_norm << std::endl;
        this->m_log->o << "Max Iteration: " << this->m_braid_settings.m_max_iter << std::endl;

        this->m_log->o << "Sync: " << this->m_braid_settings.m_sync << std::endl;
        this->m_log->o << "Residual: " << this->m_braid_settings.m_residual << std::endl;
        this->m_log->o << "Cycle FMG: " << this->m_braid_settings.m_cycle_fmg << std::endl;
        this->m_log->o << "Cycle NFMG: " << this->m_braid_settings.m_cycle_nfmg << std::endl;
        this->m_log->o << "Cycle NFMGV: " << this->m_braid_settings.m_cycle_nfmgv << std::endl;
        this->m_log->o << "Store Values: " << this->m_braid_settings.m_store_values << std::endl;
        this->m_log->o << "Print File: " << this->m_braid_settings.m_print_file << std::endl;
        this->m_log->o << "Print Level: " << this->m_braid_settings.m_print_level << std::endl;
        this->m_log->o << "File IO Level: " << this->m_braid_settings.m_file_io_level << std::endl;
        this->m_log->o << "Access Level: " << this->m_braid_settings.m_access_level << std::endl;
        this->m_log->o << "Sequential: " << this->m_braid_settings.m_sequential << std::endl;

        this->m_log->o << "Coarsen & Refine: " << this->m_braid_settings.m_coarsen_and_refine << std::endl;
        this->m_log->o << "Min Coarse: " << this->m_braid_settings.m_min_coarse << std::endl;

        this->m_log->o << "Final FC Relax: " << this->m_braid_settings.m_final_fc_relax << std::endl;
        this->m_log->o << "Refine: " << this->m_braid_settings.m_refine << std::endl;
        this->m_log->o << "Max Refinements: " << this->m_braid_settings.m_max_refinements << std::endl;
        this->m_log->o << "Relax only CG: " << this->m_braid_settings.m_relax_only_cg << std::endl;
        this->m_log->o << "t-point cutoff: " << this->m_braid_settings.m_t_points_cutoff  << std::endl;

        this->m_log->o << "Braid ID: " << this->get_id() << std::endl;
        this->m_log->o << "Distribution: " << this->get_distribution_lower() << " , " << this->get_distribution_upper() << std::endl;

        this->m_log->o << "Number of Level: " << this->get_num_level()  << std::endl;
        this->m_log->o << "Number of Iteration: " << this->get_num_iteration()  << std::endl;

        this->m_log->o << "Richardson-based Error Estimation: " << this->m_braid_settings.m_est_error << std::endl;
        this->m_log->o << "Richardson-based Extrapolation: " << this->m_braid_settings.m_richardson << std::endl;
        this->m_log->o << "Local order: " << this->m_braid_settings.m_local_order << std::endl;
        this->m_log->o << "Reverted Ranks: " << this->m_braid_settings.m_reverted_ranks << std::endl;

        this->m_log->o << "Periodic: " << this->m_braid_settings.m_periodic << std::endl;
        this->m_log->o << "Full residual norm: " << this->m_braid_settings.m_full_r_norm << std::endl;
        this->m_log->o << "----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----" << std::endl;

        this->m_log->o << "default" << m_braid_settings.m_n_relax_default << "    "
                               << m_braid_settings.m_c_relax_weight_default << "    "
                               << m_braid_settings.m_c_factor_default << "    "
                               <<std::endl;
        for(int i = 0; i <= this->m_braid_settings.eval_level; i++){
            this->m_log->o << i << "    ";
            if(this->m_braid_settings.m_n_relax.find(i) != this->m_braid_settings.m_n_relax.end()){
                this->m_log->o <<" "<< this->m_braid_settings.m_n_relax[i] << " ";
            } else {
                this->m_log->o <<"("<< m_braid_settings.m_n_relax_default<<")";
            }
            this->m_log->o <<  "    ";

            if(this->m_braid_settings.m_c_relax_weight.find(i) != this->m_braid_settings.m_c_relax_weight.end()){
                this->m_log->o << " " <<this->m_braid_settings.m_c_relax_weight[i] << " ";
            } else {
                this->m_log->o<< "(" << m_braid_settings.m_c_relax_weight_default << ")";
            }
            this->m_log->o << "    ";
            if(this->m_braid_settings.m_c_factor.find(i) != this->m_braid_settings.m_c_factor.end()){
                this->m_log->o << " " << this->m_braid_settings.m_c_factor[i] << " ";
            } else {
                this->m_log->o << "(" << m_braid_settings.m_c_factor_default << ")";
            }
            this->m_log->o << std::endl;

        }
        this->m_log->o << "===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====" << std::endl;
    }

    void print_summary(){
        this->m_log->o << "----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----" << std::endl;
        this->m_log->o << "Number of Level: " << this->get_num_level()  << std::endl;
        this->m_log->o << "Number of Iteration: " << this->get_num_iteration()  << std::endl;

        //   void GetRNorms(braid_Int *nrequest_ptr, braid_Real *rnorms) { braid_GetRNorms(core, nrequest_ptr, rnorms); }
        int request = this->m_braid_settings.m_max_iter;
        auto * norms = new double[this->m_braid_settings.m_max_iter];
        this->m_braid_core->GetRNorms(&request, norms);
        for(int i = 0; i < request ; i++){
            this->m_log->o << i  << ": " << norms[i] << std::endl;
        }
        // todo c factors afterwards
        this->m_log->o << "===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====" << std::endl;
    }

    int test(SPGridFunction u0,
             const char *generator,
             const char *cmp,
             SPScriptor output) {// todo ctrlx
        this->m_log->o << "debug::BraidExecutor::test(args)" << std::endl<<std::flush;
        //this->m_app->setStartVector(u0);
        // this->m_app->setVectorGenerator(BraidInitializer);
        //this->m_app->setGeneratorComponent(cmp);
        this->m_app->set_scriptor(output);
        this->m_app->init();

        FILE *file;
        file = fopen("myfile", "w");

        BraidUtil bu = BraidUtil();
        bu.TestResidual(this->m_app.get(),
                        this->m_comm->SPATIAL,
                        file,
                        0,
                        0.1);

        //this->m_app->release();
        return 1; // todo return value ?
    }

    void set_output(SPScriptor output){ // todo ctrlx
        this->m_app->set_scriptor(output);
    }

    void set_initializer(SPBraidInitializer initializer){
        this->m_app->m_initializer = initializer;
    }

    void set_norm_provider(SPSpatialNorm norm){
        this->m_app->m_norm = norm;
    }

    bool apply(SPGridFunction u_stop,
            number t_stop,
            SPGridFunction u0,
            number t0) { // todo ctrlx
        this->m_log->o << "debug::BraidExecutor::apply(args)" << std::endl<<std::flush;
        this->m_log->o << "BraidExecutor run" << std::endl;

        this->m_app->set_start_time(t0);
        this->m_app->set_end_time(t_stop);
        this->m_app->set_start_vector(u0);
        //this->m_app->init();
        //this->m_log->o << "debug::BraidExecutor::apply()#d" << std::endl<<std::flush;
        this->m_braid_core->Drive();
        this->m_log->o << "debug::BraidExecutor::apply()#e" << std::endl<<std::flush;
        //this->m_app->release();
        return true;
    }

    void set_paralog(SPParalog log) {
        this->m_log = log;
        this->m_app->set_paralog(log);
    }


    void set_paralog_script(SPParalog log){
        this->m_app->set_paralog_script(log);
    }
};



#endif //UG_PLUGIN_XBRAIDFORUG4_BRAIDEXECUTOR_H

