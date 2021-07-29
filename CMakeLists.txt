cmake_minimum_required(VERSION 3.0)
project(charm_3d)

if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release")
   message(STATUS "Build type not specified: Use Release by default")
endif(NOT CMAKE_BUILD_TYPE)


set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_CXX_COMPILER "mpicxx")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O2") # for production
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g -DCHARM_DEBUG") # for developing
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2") # for production
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -DCHARM_DEBUG") # for developing

if (DEFINED POGGI)
    message("\nWARNING: Building for Poggi's task for RMI\n")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DPOGGI") # for RMI modeling
endif()

set(CONTRIB_DIR "${CMAKE_SOURCE_DIR}/contrib")

set(P4EST_INSTALL "${CONTRIB_DIR}/p4est/local")
set(P4EST_INCLUDE "${P4EST_INSTALL}/include" "${P4EST_INSTALL}/include/sc_builtin")
set(P4EST_LIB "${P4EST_INSTALL}/lib")

set(YAML_INCLUDE "${CONTRIB_DIR}/yaml-cpp/local/include")

set(YAML_LIB "${CONTRIB_DIR}/yaml-cpp/local/lib")
set(YAML_LIB64 "${CONTRIB_DIR}/yaml-cpp/local/lib64")

set(CHARM_INCLUDE "${CMAKE_SOURCE_DIR}/include")


include_directories(${CHARM_INCLUDE} ${P4EST_INCLUDE} ${YAML_INCLUDE})
link_directories(${P4EST_LIB} ${YAML_LIB} ${YAML_LIB64})

add_executable(	charm_3d

        include/charm_def.h
        include/charm_amr.h
        include/charm_base_func.h
        include/charm_bnd_cond.h
        include/charm_eos.h
        include/charm_fluxes.h
        include/charm_globals.h
        include/charm_limiter.h
        include/charm_models.h
        include/charm_vtk.h
        include/charm_version.h

        src/charm_3d.c
        src/charm_connectivity.cpp
        src/charm_timesteps.c
        src/charm_fluxes.c
        src/charm_amr.c
        src/charm_init.c
        src/charm_init_yaml.cpp
        src/charm_output.c
        src/charm_globals.c
        src/charm_bnd_cond.c
        src/charm_geom.c
        src/charm_base_func.c
        src/limiters/charm_limiter.c
        src/charm_eos.c
        src/charm_vtk.c
        src/limiters/charm_limiter_bj.c

        src/models/euler/charm_model.c
        src/models/euler/charm_model_euler_init.cpp

        src/models/advection/charm_model.c
        src/models/advection/charm_model_adv_init.cpp
        src/models/advection/charm_model_adv_amr.c


        src/models/ns/charm_model.c
        src/models/ns/charm_model_ns_init.cpp
        src/models/ns/charm_model_ns_conv.c
        src/models/ns/charm_model_ns_diff.c
        src/models/ns/charm_model_ns_diff_grad.c
        src/models/ns/charm_model_ns_diff_integrals.c
        src/models/ns/charm_model_ns_chem.c
        src/models/ns/charm_model_ns_turb.c
        src/models/ns/turb/charm_model_ns_turb_sst.c
        src/models/ns/charm_model_ns_geom.c src/models/ns/turb/charm_model_ns_turb_sa.c src/models/ns/turb/charm_model_ns_turb_sa_grad.c src/models/ns/charm_model_ns_init_cond.c)

target_link_libraries(charm_3d p4est sc m yaml-cpp)
set_property(TARGET charm_3d PROPERTY C_STANDARD 99)
set_property(TARGET charm_3d PROPERTY CXX_STANDARD 14)


add_executable(	tet2hex
        src/utils/tet2hex/tet2hex.cpp
)