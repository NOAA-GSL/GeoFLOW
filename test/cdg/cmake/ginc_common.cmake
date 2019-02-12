##############################################################
# CDG include files
##############################################################

 set(GHOME          ${CMAKE_CURRENT_SOURCE_DIR}/../../src/cdg)

#set(CDG_COMM_SRC   ${GHOME}/comm/xxx.cpp )

#set(CDG_UTILS_SRC  ${GHOME}/utils/xxx.cpp )


# Handle include files:

# Add the to the search path for include files --by target:
set(CDG_INC
          ${GHOME}/include
          ${GHOME}/blas
          ${GHOME}/comm
          ${GHOME}/sem
          ${GHOME}/utils
          ${GHOME}/grid
          ${GHOME}/pdes
          ${GHOME}/..
    )

# Add preproc defs:
set(DEF_ACC "") 
if (DO_OPENACC)
  set(DEF_ACC "-D_G_USE_OPENACC -D_G_AUTO_CREATE_DEV -D_G_AUTO_UPDATE_DEV")
endif()

set(DEF_DEBUG "") 
if (DO_DEBUG AND NOT DO_OPENACC)
  set(DEF_DEBUG "-D_G_BOUNDS_CHK")
endif()

set(DEF_GBLAS "") 
if (USE_GBLAS)
  set(DEF_GBLAS "-D_G_USE_GBLAS")
endif()

set(DEF_PAPI "") 
if (HAVE_PAPI)
  set(DEF_PAPI "-DHAVE_PAPI")
endif()

set(DEF_MPI "") 
if (USE_MPI)
  set(DEF_MPI "-D_G_USE_MPI")
endif()

add_definitions(
                ${DEF_ACC} ${DEF_DEBUG} ${DEF_GBLAS} ${DEF_PAPI} ${DEF_MPI}
               )
