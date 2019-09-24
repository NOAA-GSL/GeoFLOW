##############################################################
# CDG include files
##############################################################

set(GHOME          ${CMAKE_CURRENT_SOURCE_DIR})


# Handle include files:

# Add the to the search path for include files --by target:
set(CDG_INC
          ${GHOME}/include
          ${GHOME}/blas
          ${GHOME}/comm
          ${GHOME}/exec
          ${GHOME}/grid
          ${GHOME}/init
          ${GHOME}/io
          ${GHOME}/pdes
          ${GHOME}/sem
          ${GHOME}/user/bdy
          ${GHOME}/user/force
          ${GHOME}/user/state
          ${GHOME}/utils
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

set(DEF_GPTL "") 
if (USE_GPTL)
  set(DEF_GPTL "-D_G_USE_GPTL")
endif()

set(DEF_PAPI "") 
if (HAVE_PAPI)
  set(DEF_PAPI "-DHAVE_PAPI")
endif()

set(DEF_MPI "") 
if (USE_MPI)
  set(DEF_MPI "-D_G_USE_MPI")
endif()

set(DEF_DIM "") 
if (GDIM MATCHES "2")
  set(DEF_DIM "-D_G_IS2D")
endif()
if (GDIM MATCHES "3")
  set(DEF_DIM "-D_G_IS3D")
endif()

add_definitions(
                ${DEF_ACC} ${DEF_DEBUG} ${DEF_GBLAS} ${DEF_GPTL} ${DEF_PAPI} ${DEF_MPI} ${DEF_DIM}
               )
