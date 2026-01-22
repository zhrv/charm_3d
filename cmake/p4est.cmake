find_package(Git REQUIRED)
find_package(MPI REQUIRED)

message(STATUS "=====   Building P4EST   =====")

set(CHARM_P4EST_DIR ${CHARM_CONTRIB_DIR}/p4est)


if(NOT EXISTS ${CHARM_P4EST_DIR}/CMakeLists.txt)
  execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive -- ${CHARM_P4EST_DIR}
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMAND_ERROR_IS_FATAL ANY)
endif()

if(NOT EXISTS ${CHARM_P4EST_DIR}/local/lib/libp4est.a)
  if(NOT EXISTS ${CHARM_P4EST_DIR}/configure)
    execute_process(COMMAND ./bootstrap
      WORKING_DIRECTORY ${CHARM_P4EST_DIR}
      RESULT_VARIABLE _err
      OUTPUT_VARIABLE _out
      OUTPUT_STRIP_TRAILING_WHITESPACE)
  endif()
  if(APPLE) 
    set(ENV{LDFLAGS} "-L/usr/local/lib")
    set(ENV{CFLAGS} "-I/usr/local/include")
  endif()
  execute_process(COMMAND ./configure CC=mpicc --enable-mpi --enable-openmp --with-metis
    WORKING_DIRECTORY ${CHARM_P4EST_DIR}
    RESULT_VARIABLE _err
    COMMAND_ERROR_IS_FATAL ANY)
  
  if(NOT _err)
    execute_process(COMMAND make
      WORKING_DIRECTORY ${CHARM_P4EST_DIR}
      RESULT_VARIABLE _err
      COMMAND_ERROR_IS_FATAL ANY)
    if(NOT _err)
      execute_process(COMMAND make install
        WORKING_DIRECTORY ${CHARM_P4EST_DIR}
        RESULT_VARIABLE _err
        COMMAND_ERROR_IS_FATAL ANY)
    endif()
  endif()

  

else()

  message(STATUS "Building is nod required")

endif()

set(CHARM_P4EST_INCLUDE ${CHARM_P4EST_DIR}/local/include)
set(CHARM_P4EST_LIB ${CHARM_P4EST_DIR}/local/lib)
set(CHARM_P4EST_LIBRARIES p4est)

