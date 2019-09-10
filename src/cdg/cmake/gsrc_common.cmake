##############################################################
# CDG sources
##############################################################

 set(GHOME          ${CMAKE_CURRENT_SOURCE_DIR})


 set(CDG_BLAS_SRC   ${GHOME}/blas/cff_blas.F
                    ${GHOME}/blas/qmxmp.F
                    ${GHOME}/blas/dmxmp.F
                    ${GHOME}/blas/fmxmp.F
                    ${GHOME}/blas/gmtk.cpp
    )

 set(CDG_COMM_SRC   ${GHOME}/comm/gcomm.cpp 
    )

 set(CDG_GRID_SRC   ${GHOME}/grid/gdd_base.cpp
                    ${GHOME}/grid/ggrid.cpp
                    ${GHOME}/grid/ggrid_icos.cpp
                    ${GHOME}/grid/ggrid_box.cpp
                    ${GHOME}/grid/ggrid_factory.cpp
    )

 set(CDG_INIT_SRC   ${GHOME}/init/gspecbdy_factory.cpp 
    )

 set(CDG_IO_SRC     ${GHOME}/io/gio.cpp
    )

 set(CDG_SEM_SRC    ${GHOME}/sem/gelem_base.cpp
                    ${GHOME}/sem/gshapefcn_linear.cpp
                    ${GHOME}/sem/gshapefcn_embed.cpp
                    ${GHOME}/sem/gshapefcn_hostd.cpp
                    ${GHOME}/sem/gmass.cpp
                    ${GHOME}/sem/ghelmholtz.cpp
                    ${GHOME}/sem/gpdv.cpp
                    ${GHOME}/sem/gadvect.cpp
    )

 set(CDG_USER_SRC   ${GHOME}/user/bdy/ginitbdy_user.cpp 
                    ${GHOME}/user/bdy/gspecbdy_user.cpp 
                    ${GHOME}/user/bdy/gupdatebdy_user.cpp 
                    ${GHOME}/user/force/ginitfb.cpp 
                    ${GHOME}/user/force/ginitforce_direct_user.cpp 
                    ${GHOME}/user/force/ginitfps.cpp 
                    ${GHOME}/user/force/ginitfs.cpp 
                    ${GHOME}/user/force/ginitfv.cpp 
                    ${GHOME}/user/state/ginitb.cpp 
                    ${GHOME}/user/state/ginitc.cpp 
                    ${GHOME}/user/state/ginitps.cpp 
                    ${GHOME}/user/state/ginits.cpp 
                    ${GHOME}/user/state/ginitstate_direct_user.cpp 
                    ${GHOME}/user/state/ginitv.cpp 
    )

 set(CDG_UTILS_SRC  ${GHOME}/utils/gbitblock.cpp 
                    ${GHOME}/utils/geoflow.cpp
    )

# Aggregate source into CDG_SRC used in CMakeLists.txt:
 list(APPEND CDG_SRC 
                     ${CDG_BLAS_SRC} 
                     ${CDG_COMM_SRC}
                     ${CDG_GRID_SRC}
                     ${CDG_INIT_SRC}
                     ${CDG_IO_SRC}
                     ${CDG_SEM_SRC}
                     ${CDG_USER_SRC}
                     ${CDG_UTILS_SRC}
     )


