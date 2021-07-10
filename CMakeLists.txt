cmake_minimum_required(VERSION 2.8)

set(pluginName XBraidForUG4)
project(UG_PLUGIN_${pluginName})

set(UG4_XBraid_plugin_VERSION_MAJOR 0.1)
set(UG4_XBraid_plugin_VERSION_MINOR 0.0)

# ----------  ----------  ---------- ---------- Compiling libraries ---------- ----------  ----------  ----------
if (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/libs/braid/braid/libbraid.a")
	MESSAGE("Compiling XBRAID")
    add_custom_command(OUTPUT braidlib
            COMMAND make braid
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/libs/braid
            )
else ()
	MESSAGE("XBRAID Exists")
     add_custom_command(
            OUTPUT braidlib
            COMMAND touch hello.txt
        )
endif ()
add_custom_target(run ALL DEPENDS braidlib)



set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS_INIT} -fPIC")



set(SOURCES
        src/xbraidforug4_plugin.cpp
        src/util/BraidUsageTimer.h
        src/util/BraidTimer.h
        src/core/BraidExecutor.h
        src/core/BraidVectorStruct.h

        src/PBraidApp.h
        src/core/SpaceTimeCommunicator.h
        src/util/Scriptor.h

        src/util/MemoryObserver.h
        src/driver/DummyDistributionBraidApp.h
        src/driver/BraidIntegrator.h
        src/interface/generator.h
        src/driver/BraidGridFunctionBase.h
        src/driver/BraidTimeStepper.h
        src/control/BraidConv.h
        src/util/status.cpp
        src/util/status.h
        src/util/paralog.cpp
        src/util/paralog.h
        src/interface/integrator_factory.h
        src/driver/BraidIntegratorFactory.h
        src/biot/BiotScriptor.h
        src/interface/initializer.h src/initializer/start_value_initializer.h src/initializer/zero_Initializer.h src/interface/scriptor.h src/interface/spatial_norm.h src/spatial_norm/euclidian_norm.h src/integrator/ThetaIntegrator.h src/integrator/ThetaIntegratorFactory.h)


if(UG_ROOT_CMAKE_PATH)
    MESSAGE("XBraid UG Root was found")
else(UG_ROOT_CMAKE_PATH)
    set(UG_ROOT_CMAKE_PATH /pro/ug/ugcore/cmake)
endif(UG_ROOT_CMAKE_PATH)

include(${UG_ROOT_CMAKE_PATH}/ug_plugin_includes.cmake)

if (buildEmbeddedPlugins)
    EXPORTSOURCES(${CMAKE_CURRENT_SOURCE_DIR} ${SOURCES})
else (buildEmbeddedPlugins)
    add_library(${pluginName} SHARED ${SOURCES})
    #add_library(${pluginName} SHARED ${CMAKE_CURRENT_SOURCE_DIR}/libs/braid/braid/libbraid.a)
    target_link_libraries(${pluginName} ${CMAKE_CURRENT_SOURCE_DIR}/libs/braid/braid/libbraid.a)
    target_link_libraries(${pluginName} ug4)
    #add_dependencies(${pluginName} UG_PLUGIN_ConvectionDiffusion)
    LINK_DIRECTORIES(${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
    MESSAGE(${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
    # target_link_libraries(${pluginName} ConvectionDiffusion)
    # target_link_libraries(${pluginName} Limex)
endif (buildEmbeddedPlugins)




