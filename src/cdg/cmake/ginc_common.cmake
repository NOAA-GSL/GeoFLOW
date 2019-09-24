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

add_definitions(
                ${DEF_ACC} ${DEF_DEBUG} ${DEF_GBLAS} ${DEF_GPTL} ${DEF_PAPI} ${DEF_MPI} ${DEF_DIM}
               )
