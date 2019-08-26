##############################################################
# CDG Test sources 
##############################################################

 set(GTHOME ${CMAKE_CURRENT_SOURCE_DIR}/../../test/cdg/src)


 # Include only sources that contain main here:
 set(CDG_TEST_SRC
                    ${GTHOME}/gtest_ggfx.cpp
                    ${GTHOME}/gtest_mass.cpp
                    ${GTHOME}/gtest_gll.cpp
                    ${GTHOME}/gtest_vec.cpp
                    ${GTHOME}/gtest_blas.cpp
#                   ${GTHOME}/gtest_tmp.cpp
#                   ${GTHOME}/gtest_burgers.cpp
#                   ${GTHOME}/gtest_derivs.cpp
#                   ${GTHOME}/gtest_rk.cpp
#                   ${GTHOME}/gtest_helm.cpp
#                   ${GTHOME}/gtest_ggrid.cpp
#                   ${GTHOME}/gtest_gmtk.cpp
    )

 # Include misc other sources:
 set(CDG_TEST_SRC_MISC 
                    ${GTHOME}/ftest_binding.f90
                    ${GTHOME}/gtools.cpp
    )

