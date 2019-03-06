##############################################################
# CDG Test sources 
##############################################################

 set(GTHOME ${CMAKE_CURRENT_SOURCE_DIR}/src)


 # Include only sources that contain main here:
 set(CDG_TEST_MAIN  
                    ${GTHOME}/gtest_burgers.cpp
#                   ${GTHOME}/gtest_derivs.cpp
#                   ${GTHOME}/gtest_helm.cpp
#                   ${GTHOME}/gtest_mass.cpp
#                   ${GTHOME}/gtest_ggrid.cpp
#                   ${GTHOME}/gtest_gmtk.cpp
#                   ${GTHOME}/gtest_gll.cpp
#                   ${GTHOME}/gtest_tmp.cpp
#                   ${GTHOME}/gtest_ggfx.cpp
#                   ${GTHOME}/gtest_vec.cpp
#                   ${GTHOME}/gtest_blas.cpp
#                   ${GTHOME}/gtest_contig.cpp
#                   ${GTHOME}/gtest_pingpong.cpp
#                   ${GTHOME}/gtest_gmm.cpp
#                   ${GTHOME}/gtest_class_acc.cpp
#                   ${GTHOME}/gtest_binding.cpp
#                   ${GTHOME}/gtest_vec_access.cpp
    )

 # Include misc other sources:
 set(CDG_TEST_SRC_MISC ${GTHOME}/ftest_binding.f90
    )

 # Fill main source list (don't include CDG_TEST_MAIN):
 list(APPEND CDG_SRC 
                     ${CDG_TEST_SRC_MISC}
     )



