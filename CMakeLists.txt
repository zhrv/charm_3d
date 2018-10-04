cmake_minimum_required(VERSION 2.6)
project(charm_dg)

set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -DCHARM_DEBUG") # for developing
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3") # for production

set(CONTRIB_DIR "${CMAKE_SOURCE_DIR}/contrib")

set(P4EST_INSTALL "${CONTRIB_DIR}/p4est/local")
set(P4EST_INCLUDE "${P4EST_INSTALL}/include" "${P4EST_INSTALL}/include/sc_builtin")
set(P4EST_LIB "${P4EST_INSTALL}/lib")

set(MPI_INCLUDE "/usr/include/openmpi")

set(MXML_INCLUDE "${CONTRIB_DIR}/mxml/local/include")

set(MXML_LIB "${CONTRIB_DIR}/mxml/local/lib")

include_directories(${P4EST_INCLUDE} ${MPI_INCLUDE} ${MXML_INCLUDE})
link_directories(${P4EST_LIB} ${MXML_LIB})

add_executable(	charm_dg
        src/charm_3d.c
		src/charm_connectivity.c
        src/charm_connectivity.h
        src/charm_globals.h
        src/charm_timestep.c
        src/charm_timestep.h
        src/charm_fluxes.c
        src/charm_fluxes.h
        src/charm_amr.c
        src/charm_amr.h
        src/charm_init.c
        src/charm_init.h
        src/charm_output.c
        src/charm_output.h
        src/charm_globals.c
        src/charm_bnd_cond.c
        src/charm_bnd_cond.h
        src/charm_geom.c
        src/charm_geom.h
		src/charm_xml.c
        src/charm_xml.h
		src/charm_base_func.c src/charm_base_func.h src/charm_timestep_conv.c src/charm_timestep_conv.h src/charm_timestep_diff.c src/charm_timestep_diff.h src/charm_limiter.c)

target_link_libraries(charm_dg p4est sc m mxml)
