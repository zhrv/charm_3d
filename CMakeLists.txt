cmake_minimum_required(VERSION 3.0)
project(charm_3d)

set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_CXX_COMPILER "mpicxx")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O2") # for production
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g -DCHARM_DEBUG") # for developing
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -std=c++14") # for production
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -DCHARM_DEBUG -std=c++14") # for developing

if (DEFINED POGGI)
    message("\nWARNING: Building for Poggi's task for RMI\n")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DPOGGI") # for RMI modeling
endif()

set(CONTRIB_DIR "${CMAKE_SOURCE_DIR}/contrib")

set(P4EST_INSTALL "${CONTRIB_DIR}/p4est/local")
set(P4EST_INCLUDE "${P4EST_INSTALL}/include" "${P4EST_INSTALL}/include/sc_builtin")
set(P4EST_LIB "${P4EST_INSTALL}/lib")

set(MPI_INCLUDE "/usr/include/openmpi")

set(MXML_INCLUDE "${CONTRIB_DIR}/mxml/local/include")

set(MXML_LIB "${CONTRIB_DIR}/mxml/local/lib")

set(YAML_INCLUDE "${CONTRIB_DIR}/yaml-cpp/local/include")

set(YAML_LIB "${CONTRIB_DIR}/yaml-cpp/local/lib")

set(CHARM_INCLUDE "${CMAKE_SOURCE_DIR}/include")

include_directories(${CHARM_INCLUDE} ${P4EST_INCLUDE} ${MPI_INCLUDE} ${MXML_INCLUDE} ${YAML_INCLUDE})
link_directories(${P4EST_LIB} ${MXML_LIB} ${YAML_LIB})

add_executable(	charm_3d

        include/charm_amr.h
        include/charm_base_func.h
        include/charm_bnd_cond.h
        include/charm_eos.h
        include/charm_fluxes.h
        include/charm_globals.h
        include/charm_limiter.h
        include/charm_models.h
        include/charm_vtk.h
        include/charm_xml.h

        src/charm_3d.c
        src/charm_connectivity.c
        src/charm_timesteps.c
        src/charm_fluxes.c
        src/charm_amr.c
        src/charm_init.c
        src/charm_init_yaml.cpp
        src/charm_output.c
        src/charm_globals.c
        src/charm_bnd_cond.c
        src/charm_geom.c
        src/charm_xml.c
        src/charm_base_func.c
        src/limiters/charm_limiter.c
        src/charm_eos.c
        src/charm_vtk.c
        src/limiters/charm_limiter_bj.c

        src/models/euler/charm_model.c
        src/models/euler/charm_model_euler_init.cpp

        src/models/ns/charm_model.c
        src/models/ns/charm_model_ns_init.cpp
        src/models/ns/charm_model_ns_conv.c
        src/models/ns/charm_model_ns_diff.c
        src/models/ns/charm_model_ns_diff_grad.c
        src/models/ns/charm_model_ns_diff_integrals.c

        src/models/ns/charm_model_ns_chem.c
)

target_link_libraries(charm_3d p4est sc m mxml yaml-cpp)
set_property(TARGET charm_3d PROPERTY C_STANDARD 99)
