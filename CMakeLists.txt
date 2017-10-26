cmake_minimum_required(VERSION 2.6)
project(chamr_3d)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")

set(CONTRIB_DIR "${CMAKE_SOURCE_DIR}/contrib")

set(P4EST_INSTALL "${CONTRIB_DIR}/p4est/local")
set(P4EST_INCLUDE "${P4EST_INSTALL}/include" "${P4EST_INSTALL}/include/sc_builtin")
set(P4EST_LIB "${P4EST_INSTALL}/lib")

set(MPI_INCLUDE "/usr/include/openmpi")

include_directories(${P4EST_INCLUDE} ${MPI_INCLUDE})
link_directories(${P4EST_LIB})

add_executable(	chamr_3d
	src/chamr_3d.c
		src/charm_connectivity.c src/charm_connectivity.h src/charm_globals.h src/charm_timestep.c src/charm_timestep.h src/charm_fluxes.c src/charm_fluxes.h src/charm_amr.c src/charm_amr.h src/charm_init.c src/charm_init.h src/charm_output.c src/charm_output.h src/charm_globals.c src/charm_bnd_cond.c src/charm_bnd_cond.h)

target_link_libraries(chamr_3d p4est sc mpi m)
