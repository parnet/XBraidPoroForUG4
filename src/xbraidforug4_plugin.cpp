//  zzzzz
// Created by parnet on 12.05.19.
//


#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"


#include "control/BraidConv.h"

#include "core/BraidExecutor.h"
#include "core/BraidVectorStruct.h"
#include "core/SpaceTimeCommunicator.h"

#include "driver/BraidGridFunctionBase.h"
#include "driver/BraidIntegrator.h"
#include "driver/BraidIntegratorFactory.h"
#include "driver/BraidTimeStepper.h"
//#include "driver/DummyDistribution.h"

#include "initializer/start_value_initializer.h"
#include "initializer/zero_Initializer.h"

#include "interface/generator.h"
#include "interface/initializer.h"
#include "interface/integrator_factory.h"
#include "interface/scriptor.h"
#include "interface/spatial_norm.h"
#include "spatial_norm/euclidian_norm.h"

#include "integrator/ThetaIntegrator.h"
#include "integrator/ThetaIntegratorFactory.h"
#include "integrator/CachedThetaIntegrator.h"
#include "integrator/CachedThetaIntegratorFactory.h"

#include "util/BraidUsageTimer.h"
#include "util/MemoryObserver.h"
#include "util/paralog.h"
#include "util/Scriptor.h"
#include "util/BraidTimer.h"
#include "util/status.h"

using namespace std;
using namespace ug::bridge;

namespace ug {

    namespace XBraidForUG4 {

        struct Functionality {


            template<typename TDomain, typename TAlgebra>
            static void DomainAlgebra(Registry &reg, string grp) {

                string suffix = GetDomainAlgebraSuffix<TDomain, TAlgebra>();
                string tag = GetDomainAlgebraTag<TDomain, TAlgebra>();



                // Initializer Classes
                {
                    // BraidInitializer
                    {
                        typedef BraidInitializer<TDomain, TAlgebra> TBraidInitializer;
                        typedef ug::GridFunction<TDomain, TAlgebra> TGridFunction;
                        typedef SmartPtr<TGridFunction> SPGridFunction;

                        string name = string("BraidInitializer").append(suffix);
                        reg.add_class_<TBraidInitializer>(name, grp)
                                .set_construct_as_smart_pointer(true);
                        reg.add_class_to_group(name, "BraidInitializer", tag);
                    }

                    //start_value_initializer
                    {
                        typedef StartValueInitializer<TDomain, TAlgebra> TStartValueInitializer;
                        typedef BraidInitializer<TDomain, TAlgebra> TBraidInitializer;

                        typedef ug::GridFunction<TDomain, TAlgebra> TGridFunction;
                        typedef SmartPtr<TGridFunction> SPGridFunction;

                        string name = string("StartValueInitializer").append(suffix);
                        reg.add_class_<TStartValueInitializer, TBraidInitializer>(name, grp)
                                .add_constructor()
                                        //.add_method("init", static_cast<void (TStartValueInitializer::*)(SPGridFunction ,number)>(&TStartValueInitializer::init),"None", "init", "adds a Scriptor to the storage")
                                .add_method("set_start_vector", &TStartValueInitializer::set_start_vector, "None",
                                            "Gridfunction u0", "set the vector for t=t0")
                                .set_construct_as_smart_pointer(true);
                        reg.add_class_to_group(name, "StartValueInitializer", tag);
                    }
                    //ZeroValueInitializer
                    {
                        typedef ZeroValueInitializer<TDomain, TAlgebra> TZeroValueInitializer;
                        typedef BraidInitializer<TDomain, TAlgebra> TBraidInitializer;

                        typedef ug::GridFunction<TDomain, TAlgebra> TGridFunction;
                        typedef SmartPtr<TGridFunction> SPGridFunction;

                        string name = string("ZeroValueInitializer").append(suffix);
                        reg.add_class_<TZeroValueInitializer, TBraidInitializer>(name, grp)
                                .add_constructor()
                                        //.add_method("init", static_cast<void (TZeroValueInitializer::*)(SPGridFunction ,double)>(&TZeroValueInitializer::init),"None", "a derived class of Scriptor", "adds a Scriptor to the storage")
                                .add_method("set_start_time", &TZeroValueInitializer::set_start_time, "None",
                                            "a derived class of Scriptor", "adds a Scriptor to the storage")
                                .add_method("set_start_vector", &TZeroValueInitializer::set_start_vector, "None",
                                            "Gridfunction u0", "set the vector for t=t0")
                                .set_construct_as_smart_pointer(true);
                        reg.add_class_to_group(name, "ZeroValueInitializer", tag);
                    }
                    // todo add Generator / Initializer
                    //.add_method("setGeneratorComponent", &TBraidGridFunctionBase::setGeneratorComponent, "None","function component", "set the function component for interpolation")
                    //.add_method("setVectorGenerator",static_cast<void (TBraidGridFunctionBase::*)(const char *)>(&TBraidGridFunctionBase::setVectorGenerator),"None", "LuaFunction name","set vector for time t != t0 as guess to speed up iterations")
                    //.add_method("setVectorGenerator",static_cast<void (TBraidGridFunctionBase::*)(SPData)>(&TBraidGridFunctionBase::setVectorGenerator),"None", "UserData", "set vector for time t != t0 as guess to speed up iterations")
                    //.add_method("setVectorGenerator", static_cast<void (TBraidGridFunctionBase::*)(ug::LuaFunctionHandle)>(&TBraidGridFunctionBase::setVectorGenerator), "None", "LuaFunctionHandle","set vector for time t != t0 as guess to speed up iterations");
                }

                // Scriptor Classes
                {
                    // Scriptor
                    {
                        typedef Scriptor<TDomain, TAlgebra> TScriptor;
                        string name = string("Scriptor").append(suffix);
                        reg.add_class_<TScriptor>(name, grp)
                                .add_method("lua_write",&TScriptor::lua_write,"","","")
                                .set_construct_as_smart_pointer(true);
                        reg.add_class_to_group(name, "Scriptor", tag);
                    }
                    // VTKScriptor
                    {
                        typedef Scriptor<TDomain, TAlgebra> TScriptor;
                        typedef SmartPtr<ug::VTKOutput<TDomain::dim>> SPVTKOutput;
                        typedef VTKScriptor<TDomain, TAlgebra> TVTKScriptor;
                        string name_vtk = string("VTKScriptor").append(suffix);
                        reg.add_class_<TVTKScriptor, TScriptor>(name_vtk, grp)
                                .template add_constructor<void (*)(SPVTKOutput, const char *)>("")
                                //.add_method("print", &TVTKScriptor::print, "None", "Gridfunction u# index # timestamp",
                                //            "writes a VTK file for the grid function")
                                .add_method("write_time_pvd", &TVTKScriptor::write_time_pvd, "None", "Gridfunction u",
                                            "creates the PVD for the series")
                                .set_construct_as_smart_pointer(true);;
                        reg.add_class_to_group(name_vtk, "VTKScriptor", tag);
                    }
                    //VTKIterationScriptor
                    {
                        typedef Scriptor<TDomain, TAlgebra> TScriptor;
                        typedef SmartPtr<ug::VTKOutput<TDomain::dim>> SPVTKOutput;
                        typedef VTKIterationScriptor<TDomain, TAlgebra> TVTKIScriptor;
                        string name_vtk = string("VTKIterationScriptor").append(suffix);
                        reg.add_class_<TVTKIScriptor, TScriptor>(name_vtk, grp)
                                .template add_constructor<void (*)(SPVTKOutput, const char *)>("")
                                //.add_method("print", &TVTKIScriptor::print, "None", "Gridfunction u# index # timestamp",
                                //            "writes a VTK file for the grid function")
                                .add_method("write_time_pvd", &TVTKIScriptor::write_time_pvd, "None", "Gridfunction u",
                                            "creates the PVD for the series")
                                .set_construct_as_smart_pointer(true);;
                        reg.add_class_to_group(name_vtk, "VTKIterationScriptor", tag);
                    }

                    // VTKModScriptor
                    {
                        typedef Scriptor<TDomain, TAlgebra> TScriptor;
                        typedef VTKModScriptor<TDomain, TAlgebra> TVTKModScriptor;
                        typedef SmartPtr<ug::VTKOutput<TDomain::dim>> SPVTKOutput;
                        string name_vtkmod = string("VTKModScriptor").append(suffix);
                        reg.add_class_<TVTKModScriptor, TScriptor>(name_vtkmod, grp)
                                .template add_constructor<void (*)(SPVTKOutput, const char *)>("")
                                //.add_method("print", &TVTKModScriptor::print, "None",
                                //            "Gridfunction u# index # timestamp",
                                //            "writes a VTK file for the grid function")
                                .add_method("setModal", &TVTKModScriptor::set_modal, "None", "n",
                                            "write only every n'th value")
                                .add_method("write_time_pvd", &TVTKModScriptor::write_time_pvd, "None",
                                            "Gridfunction u", "creates the PVD for the series")
                                .set_construct_as_smart_pointer(true);
                        reg.add_class_to_group(name_vtkmod, "VTKModScriptor", tag);
                    }

                    // EvalScriptor
                    /*{
                        typedef Scriptor<TDomain, TAlgebra> TScriptor;
                        typedef EvalScriptor<TDomain, TAlgebra> TEvalScriptor;
                        string name_eval = string("EvalScriptor").append(suffix);
                        reg.add_class_<TEvalScriptor, TScriptor>(name_eval, grp)
                                .add_constructor()
                                .add_method("setFile",static_cast<void (TEvalScriptor::*)(const char *)>(&TEvalScriptor::set_file),"void","filename", "Creates a file")
                                .add_method("setGeneratorComponent", &TEvalScriptor::set_generator_component, "None","function component", "Set the function component for interpolation")
                                .add_method("setVectorGenerator", static_cast<void (TEvalScriptor::*)(const char *)>(&TEvalScriptor::set_vector_generator), "None", "lua function name","Set a LUA function to generate Gridfunctions to different timestamps")
                                //.add_method("setVectorGenerator",static_cast<void (TEvalScriptor::*)(SPData)>(&TEvalScriptor::setVectorGenerator),"None", "userdata", "Set user data to generate Gridfunctions")
                                .add_method("setDomain", &TEvalScriptor::set_domain, "None", "domain discretization","Set the Domain in which the Gridfunctions are generated")
                                .add_method("setRelative", &TEvalScriptor::set_relative, "None","display norm of u and norm of norm(u-v)/nrom(u)", "")
                                .add_method("print", &TEvalScriptor::print, "None", "Gridfunction u# index # timestamp","generates a Gridfunction and compare this with the given Gridfunction")
                                .add_method("write_time_pvd", &TEvalScriptor::write_time_pvd, "None", "","finishes the access")
                                .set_construct_as_smart_pointer(true);
                        reg.add_class_to_group(name_eval, "EvalScriptor", tag);
                    }*/

                    // MATLAB Scriptor
                    {
                    }
                    // MultiScriptor
                    {
                        typedef Scriptor<TDomain, TAlgebra> TScriptor;
                        typedef MultiScriptor<TDomain, TAlgebra> TMultiScriptor;
                        string name_multi = string("MultiScriptor").append(suffix);
                        reg.add_class_<TMultiScriptor, TScriptor>(name_multi, grp)
                                .add_constructor()
                                .add_method("addScriptor", &TMultiScriptor::add_scriptor, "None",
                                            "a derived class of Scriptor", "adds a Scriptor to the storage")
                                //.add_method("print", &TMultiScriptor::print, "None", "Gridfunction u#index#timestamp",
                                //            "accesses a given gridfunction and delegate it to the stored Scriptors")
                                .add_method("write_time_pvd", &TMultiScriptor::write_time_pvd, "None", "Gridfunction u",
                                            "finishes the access and write a pvd")
                                .set_construct_as_smart_pointer(true);
                        reg.add_class_to_group(name_multi, "MultiScriptor", tag);
                    }
                }
                // Braid Convergence Control
                {
                    // Interface
                    {

                    }
                    // Static
                    {

                    }
                    // Adaptive
                    {
                        /*
                        .add_method("setTightTol", &TRGFBraidApp::setTightTol, "None", "tight tol","sets the tight tolerance for adaptive convergence check")
                                .add_method("setLooseTol", &TRGFBraidApp::setLooseTol, "None", "loose tol","sets the loose tolerance for adaptive convergence check")
                                .add_method("setAdaptConv", &TRGFBraidApp::setAdaptConv, "None", "adaptivity","set adaptive convergence check")
                                .add_method("setForceConv", &TRGFBraidApp::setForceConv, "None", "force convergence","set cancel on convergence check fail")
                                .add_method("setStrongFirstIteration", &TRGFBraidApp::setStrongFirstIteration, "None", "value","apply the first iteration with tight tolerance")
                                .add_method("setRecurringStrongIteration", &TRGFBraidApp::setRecurringStrongIteration, "None","interval", "apply every n'th iteration with tight tolerance")*/
                    }

                }
                // Norm Provider
                {
                    //Norm Provider
                    {
                        typedef BraidSpatialNorm<TDomain, TAlgebra> TBraidSpatialNorm;
                        string name_gf = string("BraidSpatialNorm").append(suffix);
                        reg.add_class_<TBraidSpatialNorm>(name_gf, grp);
                        reg.add_class_to_group(name_gf, "BraidSpatialNorm", tag);
                    }

                    //
                    {
                        typedef BraidEuclidianNorm<TDomain, TAlgebra> TBraidEuclidianNorm;
                        typedef BraidSpatialNorm<TDomain, TAlgebra> TBraidSpatialNorm;
                        string name_gf = string("BraidEuclidianNorm").append(suffix);
                        reg.add_class_<TBraidEuclidianNorm,TBraidSpatialNorm>(name_gf, grp)
                                .add_constructor()
                                .add_method("norm", &TBraidEuclidianNorm::norm, "None", "verbose",
                                            "set the level of verbose (true / false)")
                                .set_construct_as_smart_pointer(true);;
                        reg.add_class_to_group(name_gf, "BraidEuclidianNorm", tag);
                    }
                }
                // Integrator
                {
                    {
                        typedef ThetaIntegrator<TDomain, TAlgebra> TThetaIntegrator;
                        typedef ug::ITimeIntegrator<TDomain, TAlgebra> TBase;
                        string name_gf = string("ThetaIntegrator").append(suffix);
                        reg.add_class_<TThetaIntegrator, TBase>(name_gf, grp) // todo constructor set by executor
                                .add_constructor()
                                .add_method("set_domain", &TThetaIntegrator::set_domain, "None", "verbose",
                                            "set the level of verbose (true / false)")
                                .add_method("set_solver", &TThetaIntegrator::set_domain, "None", "verbose",
                                            "set the level of verbose (true / false)")
                                .set_construct_as_smart_pointer(true);
                        reg.add_class_to_group(name_gf, "ThetaIntegrator", tag);
                    }
                    {
                        typedef CachedThetaIntegrator<TDomain, TAlgebra> TThetaIntegrator;
                        typedef ug::ITimeIntegrator<TDomain, TAlgebra> TBase;
                        string name_gf = string("CachedThetaIntegrator").append(suffix);
                        reg.add_class_<TThetaIntegrator, TBase>(name_gf, grp) // todo constructor set by executor
                                .add_constructor()
                                .add_method("set_domain", &TThetaIntegrator::set_domain, "None", "verbose",
                                            "set the level of verbose (true / false)")
                                .add_method("set_solver", &TThetaIntegrator::set_domain, "None", "verbose",
                                            "set the level of verbose (true / false)")
                                .set_construct_as_smart_pointer(true);
                        reg.add_class_to_group(name_gf, "CachedThetaIntegrator", tag);
                    }
                }
                // Integrator Factory
                {
                    //SimpleLimexFactory
                    {
                        typedef IntegratorFactory<TDomain, TAlgebra> TIntegratorFactory;
                        string name_gf = string("IntegratorFactory").append(suffix);
                        reg.add_class_<TIntegratorFactory>(name_gf, grp); // todo constructor set by executor
                        reg.add_class_to_group(name_gf, "IntegratorFactory", tag);
                    }

                    // Theta Factory
                    {
                        typedef ThetaIntegratorFactory<TDomain, TAlgebra> TThetaIntegrator;
                        typedef IntegratorFactory<TDomain, TAlgebra> TBase;
                        string name_gf = string("ThetaIntegratorFactory").append(suffix);
                        reg.add_class_<TThetaIntegrator, TBase>(name_gf, grp) // todo constructor set by executor
                                .add_constructor()
                                .add_method("set_domain", &TThetaIntegrator::set_domain, "None", "verbose",
                                            "set the level of verbose (true / false)")
                                .add_method("set_solver", &TThetaIntegrator::set_solver, "None", "verbose","set the level of verbose (true / false)")
                                .add_method("create_time_integrator", &TThetaIntegrator::create_time_integrator, "None", "Gridfunction u0","set the vector for t=t0")
                                .add_method("create_level_time_integrator", &TThetaIntegrator::create_level_time_integrator, "None", "Gridfunction u0","set the vector for t=t0")
                                .set_construct_as_smart_pointer(true);
                        reg.add_class_to_group(name_gf, "ThetaIntegratorFactory", tag);
                    }

                    {
                        typedef CachedThetaIntegratorFactory<TDomain, TAlgebra> TThetaIntegrator;
                        typedef IntegratorFactory<TDomain, TAlgebra> TBase;
                        string name_gf = string("CachedThetaIntegratorFactory").append(suffix);
                        reg.add_class_<TThetaIntegrator, TBase>(name_gf, grp) // todo constructor set by executor
                                .add_constructor()
                                .add_method("set_domain", &TThetaIntegrator::set_domain, "None", "verbose",
                                            "set the level of verbose (true / false)")
                                .add_method("set_solver", &TThetaIntegrator::set_solver, "None", "verbose",
                                            "set the level of verbose (true / false)")
                                .add_method("create_time_integrator", &TThetaIntegrator::create_time_integrator, "None", "Gridfunction u0","set the vector for t=t0")
                                .add_method("create_level_time_integrator", &TThetaIntegrator::create_level_time_integrator, "None", "Gridfunction u0","set the vector for t=t0")
                                .set_construct_as_smart_pointer(true);
                        reg.add_class_to_group(name_gf, "CachedThetaIntegratorFactory", tag);
                    }
                    /*{

                        typedef IntegratorFactory<TDomain,TAlgebra> TIntegratorFactory;
                        typedef SimpleLimexFactory<TDomain, TAlgebra> TLimexFactory;
                        string name_gf = string("LimexFactory").append(suffix);
                        reg.add_class_<TLimexFactory,TIntegratorFactory>(name_gf, grp) // todo constructor set by executor
                                .add_constructor()
                                .add_method("set_domain_disc", &TLimexFactory::set_domain_disc, "None", "verbose","set the level of verbose (true / false)")
                                .add_method("set_solver", &TLimexFactory::set_solver, "None", "verbose","set the level of verbose (true / false)")
                                .add_method("set_error_estimator", &TLimexFactory::set_error_estimator, "None", "initial time","set t0 as initial time")
                                .add_method("attach_observer", &TLimexFactory::attach_observer, "None", "end time", "set tN as endtime")
                                .add_method("create_time_integrator", &TLimexFactory::create_time_integrator, "","number of timesteps", "set N as number of timesteps")
                                .set_construct_as_smart_pointer(true);
                        reg.add_class_to_group(name_gf, "LimexFactory", tag);
                    }


                    // SimpleIntegratorFactory
                    {
                        typedef IntegratorFactory<TDomain,TAlgebra> TIntegratorFactory;
                        typedef SimpleIntegratorFactory<TDomain, TAlgebra> TBraidGridFunctionBase;
                        string name_gf = string("SimpleIntegratorFactory").append(suffix);
                        reg.add_class_<TBraidGridFunctionBase,TIntegratorFactory>(name_gf, grp) // todo constructor set by executor
                                .add_constructor()
                                .add_method("set_time_stepper", &TBraidGridFunctionBase::set_time_stepper, "None", "verbose","set the level of verbose (true / false)")
                                .add_method("set_solver", &TBraidGridFunctionBase::set_solver, "None", "initial time","set t0 as initial time")
                                .add_method("set_dt_min", &TBraidGridFunctionBase::set_dt_min, "None", "end time", "set tN as endtime")
                                .add_method("set_dt_max", &TBraidGridFunctionBase::set_dt_max, "","number of timesteps", "set N as number of timesteps")
                                .add_method("set_reduction_factor", &TBraidGridFunctionBase::set_reduction_factor,"None", "t0#tN#N", "sets tstart, tstop, number of timesteps")
                                        //.add_method("setDomainDisc", &TBraidGridFunctionBase::setDomainDisc, "None", "domain discretization","set the domain")
                                .add_method("create_time_integrator", &TBraidGridFunctionBase::create_time_integrator, "None", "Gridfunction u0","set the vector for t=t0")

                                .set_construct_as_smart_pointer(true);
                        reg.add_class_to_group(name_gf, "SimpleIntegratorFactory", tag);
                    }
                }*/

                    // Braid App Classes
                    {
                        // Braid Grid Function Base
                        {
                            typedef BraidGridFunctionBase<TDomain, TAlgebra> TBraidGridFunctionBase;
                            string name_gf = string("BraidGridFunctionBase").append(suffix);
                            reg.add_class_<TBraidGridFunctionBase>(name_gf, grp) // todo constructor set by executor
                                    .add_method("init", &TBraidGridFunctionBase::init, "None", "verbose",
                                                "set the level of verbose (true / false)")
                                    .add_method("set_verbose", &TBraidGridFunctionBase::set_verbose, "None", "verbose",
                                                "set the level of verbose (true / false)")
                                    .add_method("set_start_time", &TBraidGridFunctionBase::set_start_time, "None",
                                                "initial time", "set t0 as initial time")
                                    .add_method("set_end_time", &TBraidGridFunctionBase::set_end_time, "None",
                                                "end time", "set tN as endtime")
                                    .add_method("set_number_of_timesteps",
                                                &TBraidGridFunctionBase::set_number_of_timesteps, "",
                                                "number of timesteps", "set N as number of timesteps")
                                    .add_method("set_time_values",
                                                static_cast<void (TBraidGridFunctionBase::*)(double, double,
                                                                                             int)>(&TBraidGridFunctionBase::set_time_values),
                                                "None", "t0#tN#N", "sets tstart, tstop, number of timesteps")
                                            //.add_method("setDomainDisc", &TBraidGridFunctionBase::setDomainDisc, "None", "domain discretization","set the domain")
                                    .add_method("set_start_vector", &TBraidGridFunctionBase::set_start_vector, "None",                                                "Gridfunction u0", "set the vector for t=t0")
                                    .add_method("set_norm_provider", &TBraidGridFunctionBase::set_norm_provider, "None",                                                "Gridfunction u0", "set the vector for t=t0")
                                    .add_method("set_scriptor", &TBraidGridFunctionBase::set_scriptor, "None",
                                                "Gridfunction u0", "set the vector for t=t0")
                                    .add_method("set_max_levels", &TBraidGridFunctionBase::set_max_levels, "None",
                                                "number", "set maximum number of level");
                            reg.add_class_to_group(name_gf, "BraidGridFunctionBase", tag);
                        }

                        // Braid Time Integrator
                        {
                            typedef BraidIntegrator<TDomain, TAlgebra> TBraidIntegrator;
                            typedef BraidGridFunctionBase<TDomain, TAlgebra> TBraidGridFunctionBase;
                            string name_gf = string("BraidIntegrator").append(suffix);
                            reg.add_class_<TBraidIntegrator, TBraidGridFunctionBase>(name_gf,
                                                                                     grp) // todo constructor set by executor // todo release?
                                    .add_constructor()
                                    .add_method("print_settings", &TBraidIntegrator::print_settings, "None", "verbose",
                                                "set the level of verbose (true / false)")
                                    .add_method("set_time_integrator", &TBraidIntegrator::set_time_integrator, "None",
                                                "verbose", "set the level of verbose (true / false)")

                                            //.add_method("set_adapt_convergence", &TBraidIntegrator::set_adapt_conv, "None", "initial time","set t0 as initial time")
                                    .set_construct_as_smart_pointer(true);
                            //.add_method("set_force_convergence", &TBraidIntegrator::setForceConv, "None", "end time", "set tN as endtime");
                            // todo add Generator / Initializer
                            //.add_method("setGeneratorComponent", &TBraidGridFunctionBase::setGeneratorComponent, "None","function component", "set the function component for interpolation")
                            //.add_method("setVectorGenerator",static_cast<void (TBraidGridFunctionBase::*)(const char *)>(&TBraidGridFunctionBase::setVectorGenerator),"None", "LuaFunction name","set vector for time t != t0 as guess to speed up iterations")
                            //.add_method("setVectorGenerator",static_cast<void (TBraidGridFunctionBase::*)(SPData)>(&TBraidGridFunctionBase::setVectorGenerator),"None", "UserData", "set vector for time t != t0 as guess to speed up iterations")
                            //.add_method("setVectorGenerator", static_cast<void (TBraidGridFunctionBase::*)(ug::LuaFunctionHandle)>(&TBraidGridFunctionBase::setVectorGenerator), "None", "LuaFunctionHandle","set vector for time t != t0 as guess to speed up iterations");
                            reg.add_class_to_group(name_gf, "BraidIntegrator", tag);
                        }
                        // Braid Time Integrator Factory
                        {
                            typedef BraidIntegratorFactory<TDomain, TAlgebra> TBraidIntegratorFactory;
                            typedef BraidGridFunctionBase<TDomain, TAlgebra> TBraidGridFunctionBase;
                            string name_gf = string("BraidIntegratorFactory").append(suffix);
                            reg.add_class_<TBraidIntegratorFactory, TBraidGridFunctionBase>(name_gf,
                                                                                            grp) // todo constructor set by executor // todo release?
                                    .add_constructor()
                                    .add_method("print_settings", &TBraidIntegratorFactory::print_settings, "None",
                                                "verbose", "set the level of verbose (true / false)")
                                    .add_method("set_fine_time_integrator",
                                                &TBraidIntegratorFactory::set_fine_time_integrator, "None", "verbose",
                                                "set the level of verbose (true / false)")
                                    .add_method("set_coarse_time_integrator",
                                                &TBraidIntegratorFactory::set_coarse_time_integrator, "None", "verbose",
                                                "set the level of verbose (true / false)")

                                            //.add_method("set_adapt_convergence", &TBraidIntegrator::set_adapt_conv, "None", "initial time","set t0 as initial time")
                                    .set_construct_as_smart_pointer(true);
                            //.add_method("set_force_convergence", &TBraidIntegrator::setForceConv, "None", "end time", "set tN as endtime");
                            // todo add Generator / Initializer
                            //.add_method("setGeneratorComponent", &TBraidGridFunctionBase::setGeneratorComponent, "None","function component", "set the function component for interpolation")
                            //.add_method("setVectorGenerator",static_cast<void (TBraidGridFunctionBase::*)(const char *)>(&TBraidGridFunctionBase::setVectorGenerator),"None", "LuaFunction name","set vector for time t != t0 as guess to speed up iterations")
                            //.add_method("setVectorGenerator",static_cast<void (TBraidGridFunctionBase::*)(SPData)>(&TBraidGridFunctionBase::setVectorGenerator),"None", "UserData", "set vector for time t != t0 as guess to speed up iterations")
                            //.add_method("setVectorGenerator", static_cast<void (TBraidGridFunctionBase::*)(ug::LuaFunctionHandle)>(&TBraidGridFunctionBase::setVectorGenerator), "None", "LuaFunctionHandle","set vector for time t != t0 as guess to speed up iterations");
                            reg.add_class_to_group(name_gf, "BraidIntegratorFactory", tag);
                        }
                        // Braid Time Stepper
                        {
                            typedef BraidTimeStepper<TDomain, TAlgebra> TBraidTimeStepper;
                            typedef BraidGridFunctionBase<TDomain, TAlgebra> TBraidGridFunctionBase;
                            string name_gf = string("BraidTimeStepper").append(suffix);
                            reg.add_class_<TBraidTimeStepper, TBraidGridFunctionBase>(name_gf,
                                                                                      grp) // todo constructor set by executor // todo release?
                                    .add_constructor()
                                    .add_method("print_settings", &TBraidTimeStepper::print_settings, "None", "verbose",
                                                "set the level of verbose (true / false)")
                                    .add_method("set_adapt_convergence", &TBraidTimeStepper::setAdaptConv, "None",
                                                "initial time", "set t0 as initial time")
                                    .add_method("set_domain", &TBraidTimeStepper::set_domain, "None",
                                                "initial time", "set t0 as initial time")
                                    .add_method("set_solver", &TBraidTimeStepper::set_solver, "None",
                                                "initial time", "set t0 as initial time")
                                    .add_method("set_force_convergence", &TBraidTimeStepper::setForceConv, "None",
                                                "end time", "set tN as endtime")
                                    .set_construct_as_smart_pointer(true);
                            reg.add_class_to_group(name_gf, "BraidTimeStepper", tag);
                        }
/*
.add_method("setStoreOperator", &TITSGFBraidApp::setStoreOperator, "None", "store","store operator for level l")
.add_method("setResidual", &TITSGFBraidApp::setResidual, "None", "method","set the residual method")
*/

                        // Dummy Braid App
                        {

                        }
                    }



                    // Braid Executor
                    {
                        // Grid Function Executor
                        {
                            typedef BraidExecutor<TDomain, TAlgebra> TBraidExecutor;

                            typedef SmartPtr<SpaceTimeCommunicator> SPXCommunicator;
                            typedef SmartPtr<BraidGridFunctionBase<TDomain, TAlgebra>> SPBraidGridFunctionBase;
                            typedef SmartPtr<Scriptor<TDomain, TAlgebra>> SPScriptor;

                            typedef ug::GridFunction<TDomain, TAlgebra> TGridFunction;
                            typedef SmartPtr<TGridFunction> SPGridFunction;
                            typedef ConstSmartPtr<TGridFunction> CSPGridFunction;

                            string name = string("BraidExecutor").append(suffix);
                            reg.add_class_<TBraidExecutor>(name, grp)
                                    .template add_constructor<void (*)(SPXCommunicator, SPBraidGridFunctionBase)>(
                                            "Communicator, BraidApp")
#ifdef UG_FOR_LUA
                                    .add_method("apply", static_cast<bool (TBraidExecutor::*)(SPGridFunction, number,
                                                                                              SPGridFunction,
                                                                                              number)>(&TBraidExecutor::apply),
                                                "", "", "starts XBraid iteration")
                                    .add_method("test",
                                                static_cast<int (TBraidExecutor::*)(SPGridFunction, const char *,
                                                                                    const char *,
                                                                                    SPScriptor)>(&TBraidExecutor::test),
                                                "", "",
                                                "tests XBraid core functionality, grid has to be defined at time t=0.0")
#endif

                                    .add_method("set_residual", &TBraidExecutor::set_residual, "", "",
                                                "activates the usage of the residual, if supported by the app")
                                    .add_method("set_n_relax", &TBraidExecutor::set_n_relax, "", "",
                                                "set the number of relaxation for the provided level, both values must be > 0")
                                    .add_method("set_c_factor", &TBraidExecutor::set_c_factor, "", "",
                                                "set the coarsening factor for the given level, both values must be > 0")
                                    .add_method("set_max_levels", &TBraidExecutor::set_max_levels, "", "",
                                                "set the maximum number of level for the time grid")
                                    .add_method("set_skip_downcycle_work", &TBraidExecutor::set_skip_downcycle_work, "",
                                                "", "skips the first downcycle if set to true")
                                    .add_method("set_min_coarse", &TBraidExecutor::set_min_coarse, "", "",
                                                "tests XBraid core functionality, grid has to be defined at time t=0.0")

                                    .add_method("set_max_iterations", &TBraidExecutor::set_max_iterations, "", "",
                                                "set the maximum number of allowed iterations must be > 0")
                                    .add_method("set_absolute_tol", &TBraidExecutor::set_absolute_tol, "", "",
                                                "sets the absolute tolerance > 0.0")
                                    .add_method("set_relative_tol", &TBraidExecutor::set_relative_tol, "", "",
                                                "set the relative tolerance")
                                    .add_method("set_temporal_norm", &TBraidExecutor::set_temporal_norm, "", "",
                                                "set the temporal norm which should be used by braid ")

                                    .add_method("set_sequential", &TBraidExecutor::set_sequential, "", "",
                                                "for testing, uses sequential time stepping")
                                    .add_method("set_store_values", &TBraidExecutor::set_store_values, "", "",
                                                "activates the storage of values up to the given level")
                                    .add_method("set_spatial_coarsen_and_refine",
                                                &TBraidExecutor::set_spatial_coarsen_and_refine, "", "",
                                                "activates coarsening for space and time")
                                    .add_method("set_refine", &TBraidExecutor::set_refine, "", "", "todo") // todo
                                    .add_method("set_max_refinements", &TBraidExecutor::set_max_refinements, "", "",
                                                "todo") // todo
                                    .add_method("set_access_level", &TBraidExecutor::set_access_level, "", "",
                                                "0: never access files, 1: access after iteration finished (fine grid only),  2: after each iteration on each level")
                                    .add_method("set_print_level", &TBraidExecutor::set_print_level, "", "",
                                                "0: no output, 1: residual history, 2: plust post run statics (default), 3: debug output")
                                    .add_method("set_print_file", &TBraidExecutor::set_print_file, "", "",
                                                "set a file to create the output for xbraid prints")
                                    .add_method("set_default_print_file", &TBraidExecutor::set_default_print_file, "",
                                                "", "set the printfile to the default print file ")

                                    .add_method("set_cycle_fmg", &TBraidExecutor::set_cycle_fmg, "", "",
                                                "activates f cycle instead of v cylce")
                                    .add_method("set_cycle_nfmg", &TBraidExecutor::set_cycle_nfmg, "", "",
                                                "number of f cycles before switching to v cycle ")
                                    .add_method("set_cycle_nfmgv", &TBraidExecutor::set_cycle_nfmgv, "", "",
                                                "v cycle on each level ( default is 1) ")

                                    .add_method("set_cycle_type", &TBraidExecutor::set_cycle_type, "", "",
                                                " cycletype 'V FCF', 'V F', 'V F-FCF', 'F FCF', 'F F', 'F F-FCF'")
                                    .add_method("set_sync", &TBraidExecutor::set_sync, "", "",
                                                "activates syncronisation")
                                    .add_method("set_increase_max_levels", &TBraidExecutor::set_increase_max_levels, "",
                                                "", "todo") // todo
                                    .add_method("set_relax_only_cg", &TBraidExecutor::set_relax_only_cg, "", "",
                                                "todo") // todo
                                    .add_method("set_agg_c_factor", &TBraidExecutor::set_agg_c_factor, "", "",
                                                "todo") // todo
                                    .add_method("set_periodic", &TBraidExecutor::set_periodic, "", "", "todo") // todo
                                    .add_method("set_final_fc_relax", &TBraidExecutor::set_final_fc_relax, "", "",
                                                "todo") // todo
                                    .add_method("set_reverted_ranks", &TBraidExecutor::set_reverted_ranks, "", "",
                                                "todo") // todo
                                    .add_method("set_richardson_estimation", &TBraidExecutor::set_richardson_estimation,
                                                "", "", "todo") // todo

                                    .add_method("set_file_io_level", &TBraidExecutor::set_file_io_level, "", "",
                                                "0: no output, 1: save output for braid.out.cycle")
                                    .add_method("set_c_relax_weight", &TBraidExecutor::set_c_relax_weight, "", "",
                                                "todo") // todo
                                    .add_method("set_t_points_cutoff", &TBraidExecutor::set_t_points_cutoff, "", "",
                                                "todo") // todo
                                    .add_method("set_full_residual_norm", &TBraidExecutor::set_full_residual_norm, "",
                                                "", "todo") // todo
                                    .add_method("set_time_grid", &TBraidExecutor::set_time_grid, "", "", "todo") // todo
                                    .add_method("get_num_iteration", &TBraidExecutor::get_num_iteration, "", "",
                                                "todo") // todo
                                    .add_method("get_c_factor", &TBraidExecutor::get_c_factor, "", "", "todo") // todo
                                    .add_method("get_residual_norms", &TBraidExecutor::get_residual_norms, "", "",
                                                "todo") // todo
                                    .add_method("get_num_level", &TBraidExecutor::get_num_level, "", "", "todo") // todo
                                    .add_method("get_warm_restart", &TBraidExecutor::get_warm_restart, "", "",
                                                "todo") // todo
                                    .add_method("get_distribution_lower", &TBraidExecutor::get_distribution_lower, "",
                                                "", "todo") // todo
                                    .add_method("get_distribution_upper", &TBraidExecutor::get_distribution_upper, "",
                                                "", "todo") // todo
                                    .add_method("get_id", &TBraidExecutor::get_id, "", "", "todo") // todo



                                    .add_method("set_app", &TBraidExecutor::set_app, "", "",
                                                "set vector for time t = t0")
                                    .add_method("set_paralog", &TBraidExecutor::set_paralog, "", "",
                                                "set vector for time t = t0")
                                    .add_method("get_app", &TBraidExecutor::get_app, "braid application", "",
                                                "set vector for time t = t0")


                                    .add_method("set_initializer", &TBraidExecutor::set_initializer, "", "","tests XBraid core functionality, grid has to be defined at time t=0.0")
                                    .add_method("set_norm_provider", &TBraidExecutor::set_norm_provider, "", "","tests XBraid core functionality, grid has to be defined at time t=0.0")

                                    .add_method("set_output", &TBraidExecutor::create_access, "", "",
                                                "set a VTK Output")
                                    .add_method("set_filename", &TBraidExecutor::set_filename, "", "",
                                                "set a filename for VTK Output") // todo delete
                                                .add_method("print_settings", &TBraidExecutor::print_settings, "", "","set a filename for VTK Output") // todo delete
                                    .add_method("print_summary", &TBraidExecutor::print_summary, "", "",
                                                "set a filename for VTK Output") // todo delete
                                    .set_construct_as_smart_pointer(true);
                            reg.add_class_to_group(name, "BraidExecutor", tag);
                        }

                        // Dummy Executor
                        {

                        }
                    }
                }
            }


            template<typename TDomain>
            static void Domain(Registry &reg, string grp) {
                /*
                string suffix = GetDomainSuffix<TDomain>();
                string tag = GetDomainTag<TDomain>();
                 */

            }

            template<int dim>
            static void Dimension(Registry &reg, string grp) {
                /*
                string suffix = GetDimensionSuffix<dim>();
                string tag = GetDimensionTag<dim>();
                 */
            }

            template<typename TAlgebra>
            static void Algebra(Registry &reg, string grp) {
                string suffix = GetAlgebraSuffix<TAlgebra>();
                string tag = GetAlgebraTag<TAlgebra>();
            }

            static void Common(Registry &reg, string grp) { }

        };
    } // end namespace XBraidForUG4

    extern "C" void
    InitUGPlugin_XBraidForUG4(Registry *reg, string grp) {
        using namespace XBraidForUG4;
        grp.append("XBraidForUG4");
        // Space Time Communicator
        {
            typedef SpaceTimeCommunicator TSpaceTimeCommunicator;
            string name = "SpaceTimeCommunicator";
            reg->add_class_<TSpaceTimeCommunicator>(name, "XBraid", "")
                    .add_method("split", &TSpaceTimeCommunicator::split)
                    .add_method("unsplit", &TSpaceTimeCommunicator::unsplit)

                    .add_method("get_global_rank", &TSpaceTimeCommunicator::get_global_rank)
                    .add_method("get_spatial_rank", &TSpaceTimeCommunicator::get_spatial_rank)
                    .add_method("get_temporal_rank", &TSpaceTimeCommunicator::get_temporal_rank)

                    .add_method("get_global_size", &TSpaceTimeCommunicator::get_global_size)
                    .add_method("get_spatial_size", &TSpaceTimeCommunicator::get_spatial_size)
                    .add_method("get_temporal_size", &TSpaceTimeCommunicator::get_temporal_size)
                    .add_method("sleep", &TSpaceTimeCommunicator::sleep)

                    .add_constructor()
                    .set_construct_as_smart_pointer(true);
        }
        // Paralog
        {
            typedef Paralog TParalog;
            string name = "Paralog";
            reg->add_class_<TParalog>(name, "XBraid","")
                    .add_constructor()
                    .add_method("set_file_name", &TParalog::set_file_name, "None", "Gridfunction u0",
                                "set the vector for t=t0")
                    .add_method("set_comm", &TParalog::set_comm, "None", "Gridfunction u0",
                                "set the vector for t=t0")
                    .add_method("init", &TParalog::init, "None", "Gridfunction u0", "set the vector for t=t0")
                    .add_method("release", &TParalog::release, "None", "Gridfunction u0", "set the vector for t=t0")
                    .add_method("write", &TParalog::write, "None", "Gridfunction u0", "set the vector for t=t0")
                    .set_construct_as_smart_pointer(true);
            //reg.add_class_to_group(name, "Paralog", tag);
        }

        // ReplaceStandardStream
        {
            typedef ReplaceStandardStream TReplaceStandardStream;
            string name = "ReplaceStandardStream";
            reg->add_class_<TReplaceStandardStream>(name, "XBraid","")
                    .add_constructor()
                    .add_method("apply", &TReplaceStandardStream::apply, "None", "Gridfunction u0",
                                "set the vector for t=t0")
                    .add_method("undo", &TReplaceStandardStream::undo, "None", "Gridfunction u0","set the vector for t=t0")
                    .add_method("set_space_time_comm", &TReplaceStandardStream::set_space_time_comm, "None", "Gridfunction u0","set the vector for t=t0")
                    .set_construct_as_smart_pointer(true);
            //reg.add_class_to_group(name, "Paralog", tag);
        }

        reg->add_class_<BraidTimer>("BraidTimer", grp, "Class to measure time differences")
                .add_method("start", &BraidTimer::start)
                .add_method("stop", &BraidTimer::stop)
                .add_method("get", &BraidTimer::get)
                .add_constructor()
                .set_construct_as_smart_pointer(true);

        try {
            RegisterCommon<Functionality>(*reg, grp);
            RegisterDimensionDependent<Functionality>(*reg, grp);
            RegisterDomainDependent<Functionality>(*reg, grp);
            RegisterAlgebraDependent<Functionality>(*reg, grp);
            RegisterDomainAlgebraDependent<Functionality>(*reg, grp);
        }
        UG_REGISTRY_CATCH_THROW(grp);

    }
}