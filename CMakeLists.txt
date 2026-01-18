cmake_minimum_required(VERSION 3.5)
include(cmake/git.cmake)

# project(P4EST
# LANGUAGES C
# HOMEPAGE_URL https://www.p4est.org/
# DESCRIPTION "p4est manages a collection---a forest---of octrees in parallel."
# VERSION ${PROJECT_VERSION})

project(charm_3d 
VERSION ${PROJECT_VERSION})

message(STATUS "")
message(STATUS "*************************************************")
message(STATUS "*  Building CHARM_3D version " ${CHARM_VERSION})
message(STATUS "*************************************************")
message(STATUS "")


if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release")
   message(STATUS "Build type not specified: Use Release by default")
endif(NOT CMAKE_BUILD_TYPE)

find_package(MPI)

set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -DCHARM_PROJECT_VERSION=${CHARM_VERSION}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCHARM_PROJECT_VERSION=${CHARM_VERSION}")

# set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -DP4EST_WITH_METIS")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DP4EST_WITH_METIS")

# set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -DCHARM_PROJECT_VERSION=${CMAKE_PROJECT_VERSION}")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCHARM_PROJECT_VERSION=${CMAKE_PROJECT_VERSION}")

set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O2") # for production
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g -DCHARM_DEBUG") # for developing
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2") # for production
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -DCHARM_DEBUG") # for developing

if (DEFINED POGGI)
    message("\nWARNING: Building for Poggi's task for RMI\n")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DPOGGI") # for RMI modeling
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")
find_package(METIS REQUIRED)

set(CHARM_CONTRIB_DIR "${CMAKE_SOURCE_DIR}/contrib")

include(cmake/AddGitSubmodule.cmake)
include(cmake/p4est.cmake)


add_git_submodule(${CHARM_CONTRIB_DIR}/yaml-cpp)

# add_git_submodule(${CHARM_CONTRIB_DIR}/metis)


### TODO For using MPI in p4est and sc
# set(P4EST_ENABLE_METIS on)
# set(P4EST_ENABLE_MPI on)
# set(SC_ENABLE_MPI on)
# add_git_submodule(${CHARM_CONTRIB_DIR}/p4est)


include_directories(SYSTEM ${MPI_INCLUDE_PATH} ${METIS_INCLUDE_DIRS} ${CHARM_P4EST_INCLUDE} ${CMAKE_SOURCE_DIR}/include)
link_directories(${CHARM_P4EST_LIB})

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
        src/models/ns/charm_model_ns_geom.c 
        src/models/ns/turb/charm_model_ns_turb_sa.c 
        src/models/ns/turb/charm_model_ns_turb_sa_grad.c 
        src/models/ns/charm_model_ns_init_cond.c

        src/models/ns_jfnk/charm_model_ns_jfnk_conv.c 
        src/models/ns_jfnk/charm_model_ns_jfnk_diff_grad.c 
        src/models/ns_jfnk/charm_model_ns_jfnk_diff_integrals.c 
        src/models/ns_jfnk/charm_model_ns_jfnk_diff.c 
        src/models/ns_jfnk/charm_model_ns_jfnk_geom.c 
        src/models/ns_jfnk/charm_model_ns_jfnk_init_cond.c 
        src/models/ns_jfnk/charm_model_ns_jfnk_init.cpp 
        src/models/ns_jfnk/charm_model.c
)

target_link_libraries(charm_3d 
        p4est sc m yaml-cpp ${CHARM_P4EST_LIBRARIES} ${METIS_LIBRARIES} ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES}
)

set_property(TARGET charm_3d PROPERTY C_STANDARD 99)
set_property(TARGET charm_3d PROPERTY CXX_STANDARD 14)


add_executable(	tet2hex
        src/utils/tet2hex/tet2hex.cpp
)