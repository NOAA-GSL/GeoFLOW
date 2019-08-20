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
#                   ${GHOME}/comm/ggfx.cpp
    )

 set(CDG_GRID_SRC   ${GHOME}/grid/gdd_base.cpp
                    ${GHOME}/grid/ggrid.cpp
                    ${GHOME}/grid/ggrid_icos.cpp
                    ${GHOME}/grid/ggrid_box.cpp
                    ${GHOME}/grid/ggrid_factory.cpp
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

 set(CDG_UTILS_SRC  ${GHOME}/utils/gbitblock.cpp 
                    ${GHOME}/utils/geoflow.cpp
    )

 set(GEOFLOW_TBOX_SRC ${GHOME}/../tbox/clock.cpp 
                      ${GHOME}/../tbox/command_line.cpp 
                      ${GHOME}/../tbox/error_handler.cpp 
                      ${GHOME}/../tbox/global_manager.cpp
                      ${GHOME}/../tbox/input_manager.cpp
                      ${GHOME}/../tbox/io_buffer.cpp
                      ${GHOME}/../tbox/mpixx.cpp
                      ${GHOME}/../tbox/pio.cpp
                      ${GHOME}/../tbox/property_tree.cpp
                      ${GHOME}/../tbox/tracer.cpp
    )


# Aggregate source into CDG_SRC used in CMakeLists.txt:
 list(APPEND CDG_SRC 
                     ${CDG_BLAS_SRC} 
                     ${CDG_COMM_SRC}
                     ${CDG_GRID_SRC}
                     ${CDG_IO_SRC}
                     ${CDG_SEM_SRC}
                     ${GEOFLOW_TBOX_SRC}
                     ${CDG_UTILS_SRC}
     )


