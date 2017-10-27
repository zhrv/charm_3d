cmake_minimum_required(VERSION 2.6)
project(chamr_3d)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")

set(CONTRIB_DIR "${CMAKE_SOURCE_DIR}/contrib")

set(P4EST_INSTALL "${CONTRIB_DIR}/p4est/local")
set(P4EST_INCLUDE "${P4EST_INSTALL}/include" "${P4EST_INSTALL}/include/sc_builtin")
set(P4EST_LIB "${P4EST_INSTALL}/lib")

set(MPI_INCLUDE "/usr/include/openmpi")

set(YAML_INCLUDE "${CONTRIB_DIR}/yaml-0.1.7/include")

set(YAML_SRC_DIR "${CONTRIB_DIR}/yaml-0.1.7/src")

include_directories(${P4EST_INCLUDE} ${MPI_INCLUDE} ${YAML_INCLUDE})
link_directories(${P4EST_LIB})

add_executable(	chamr_3d
        ${YAML_SRC_DIR}/api.c
        ${YAML_SRC_DIR}/dumper.c
        ${YAML_SRC_DIR}/emitter.c
        ${YAML_SRC_DIR}/loader.c
        ${YAML_SRC_DIR}/parser.c
        ${YAML_SRC_DIR}/reader.c
        ${YAML_SRC_DIR}/scanner.c
        ${YAML_SRC_DIR}/writer.c
        src/chamr_3d.c
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
        src/charm_grad.c
        src/charm_grad.h
)

target_link_libraries(chamr_3d p4est sc mpi m)
