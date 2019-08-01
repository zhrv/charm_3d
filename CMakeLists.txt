cmake_minimum_required(VERSION 2.6)
project(charm_3d)

set(CMAKE_C_COMPILER "mpicc")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -DCHARM_DEBUG") # for developing
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -DPOGGI") # for production

set(CONTRIB_DIR "${CMAKE_SOURCE_DIR}/contrib")

set(P4EST_INSTALL "${CONTRIB_DIR}/p4est/local")
set(P4EST_INCLUDE "${P4EST_INSTALL}/include" "${P4EST_INSTALL}/include/sc_builtin")
set(P4EST_LIB "${P4EST_INSTALL}/lib")

set(MPI_INCLUDE "/usr/include/openmpi")

set(MXML_INCLUDE "${CONTRIB_DIR}/mxml/local/include")

set(MXML_LIB "${CONTRIB_DIR}/mxml/local/lib")

set(CHARM_INCLUDE "${CMAKE_SOURCE_DIR}/include")

include_directories(${CHARM_INCLUDE} ${P4EST_INCLUDE} ${MPI_INCLUDE} ${MXML_INCLUDE})
link_directories(${P4EST_LIB} ${MXML_LIB})

add_executable(	charm_3d
        src/charm_3d.c
		src/charm_connectivity.c
        src/charm_timesteps.c
		src/charm_fluxes.c
        src/charm_amr.c
        src/charm_init.c
		src/charm_output.c
		src/charm_globals.c
        src/charm_bnd_cond.c
        src/charm_geom.c
		src/charm_xml.c
		src/charm_base_func.c
		src/limiters/charm_limiter.c
		src/charm_eos.c
		src/charm_vtk.c
		src/models/euler/charm_model.c
		src/limiters/charm_limiter_bj.c)

target_link_libraries(charm_3d p4est sc m mxml)
set_property(TARGET charm_3d PROPERTY C_STANDARD 90)
