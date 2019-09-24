# Print compilers being used:
message("========================>>>>>>>>>>>>>>>>> CMAKE_C_COMPILER=${CMAKE_C_COMPILER}")
message("========================>>>>>>>>>>>>>>>>> CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
message("========================>>>>>>>>>>>>>>>>> CMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}")
message("==============================================================================")
message("========================>>>>>>>>>>>>>>>>> CMAKE_C_COMPILER_ID=${CMAKE_C_COMPILER_ID}")
message("========================>>>>>>>>>>>>>>>>> CMAKE_CXX_COMPILER_ID=${CMAKE_CXX_COMPILER_ID}")
message("========================>>>>>>>>>>>>>>>>> CMAKE_Fortran_COMPILER_ID=${CMAKE_Fortran_COMPILE_ID}")

##############################################################
# GNU compilers
##############################################################
if   ( ${CMAKE_C_COMPILER_ID} STREQUAL "GNU")
  set(OMP_FLAGS "" )
  if   (DO_OPENMP)
    set( OMP_FLAGS "-fopenmp" )
  endif()

  set(INSTR_FLAGS "" )
  if   (DO_AUTO_PROF)
    set( INSTR_FLAGS "-finstrument-functions" )
  endif()

  set(EXTRA_FLAGS "${OMP_FLAGS} ${INSTR_FLAGS}")

  set(CMAKE_C_FLAGS_RELEASE "-O0 -std=c++11 -finline-limit=2000 --param max-inline-insns-single=100000 -Winline ${EXTRA_FLAGS}")
  set(CMAKE_C_FLAGS_DEBUG   "-g  -std=c++11 -finline-limit=2000 --param max-inline-insns-single=100000 -Winline ${EXTRA_FLAGS}")
  set(CMAKE_C_FLAGS         ${CMAKE_C_FLAGS_RELEASE})

  set(CMAKE_CXX_FLAGS_RELEASE "-O0 -finline-limit=2000 --param max-inline-insns-single=100000  -Winline ${EXTRA_FLAGS}")
  set(CMAKE_CXX_FLAGS_DEBUG   "-g  -finline-limit=2000 --param max-inline-insns-single=100000 -Winline ${EXTRA_FLAGS}")
  set(CMAKE_CXX_FLAGS         ${CMAKE_CXX_FLAGS_RELEASE})
  use_cxx11()

  set(CMAKE_Fortran_FLAGS_RELEASE "-O0 -finline-limit=2000 --param max-inline-insns-single=100000 -Winline ${EXTRA_FLAGS}")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-g  -finline-limit=2000 --param max-inline-insns-single=100000 -Winline ${EXTRA_FLAGS}")
  set(CMAKE_Fortran_FLAGS         ${CMAKE_Fortran_FLAGS_RELEASE})
endif()

##############################################################
# Intel compiler
##############################################################
if (${CMAKE_C_COMPILER_ID} STREQUAL "Intel")
  set( OMP_FLAGS "" )
  if   (DO_OPENMP)
    set( OMP_FLAGS "-fopenmp" )
  endif()

  set(INSTR_FLAGS "" )
  if   (DO_AUTO_PROF)
    set( INSTR_FLAGS "-finstrument-functions" )
  endif()

  set(EXTRA_FLAGS "${OMP_FLAGS} ${INSTR_FLAGS}")

  set(CMAKE_C_FLAGS_RELEASE "-O2    -std=c++11 -diag-disable 161 ${EXTRA_FLAGS}")
  set(CMAKE_C_FLAGS_DEBUG   "-O0 -g -std=c++11 -diag-disable 161 ${EXTRA_FLAGS}" )
  set(CMAKE_C_FLAGS         ${CMAKE_C_FLAGS_RELEASE})

  set(CMAKE_CXX_FLAGS_RELEASE "-O2    -std=c++11 -diag-disable 161 ${EXTRA_FLAGS}")
  set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g -std=c++11 -diag-disable 161 ${EXTRA_FLAGS}" )
  set(CMAKE_CXX_FLAGS         ${CMAKE_CXX_FLAGS_RELEASE})

  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 ${EXTRA_FLAGS}")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g ${EXTRA_FLAGS}" )
  set(CMAKE_Fortran_FLAGS         ${CMAKE_Fortran_FLAGS_RELEASE})
endif ()

##############################################################
# PGI compiler
##############################################################
if (${CMAKE_C_COMPILER_ID} STREQUAL "PGI")
  set( OMP_FLAGS "" )
  if   (DO_OPENMP)
    set( OMP_FLAGS "-mp" )
  endif()

  set(ACC_FLAGS "" )
  set(ACC_LINK_FLAGS "" )
  if   (DO_OPENACC)
    set( ACC_FLAGS "-acc -ta=tesla:cuda8.0" )
    set( ACC_LINK_FLAGS "-ta=tesla:cuda8.0" )
  endif()

  set(INSTR_FLAGS "" )
  if   (DO_AUTO_PROF)
    set( INSTR_FLAGS "-Minstrument:functions" )
  endif()

  set(EXTRA_FLAGS "${OMP_FLAGS} ${ACC_FLAGS} ${INSTR_FLAGS}")

  set(CMAKE_C_FLAGS_RELEASE "-O2    -Minfo ${EXTRA_FLAGS}")
  set(CMAKE_C_FLAGS_DEBUG   "-O0 -g -Minfo ${EXTRA_FLAGS}")
  set(CMAKE_C_FLAGS         ${CMAKE_C_FLAGS_RELEASE})

  set(CMAKE_CXX_FLAGS_RELEASE "-O2    -std=c++11 -Minfo ${OMP_FLAGS} ${ACC_FLAGS ${INSTR_FLAGS}}")
  set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g -std=c++11 -Minfo ${OMP_FLAGS} ${ACC_FLAGS ${INSTR_FLAGS}}" )
  set(CMAKE_CXX_FLAGS         ${CMAKE_CXX_FLAGS_RELEASE})

  set(CMAKE_Fortran_FLAGS_RELEASE "-O2  -Minfo ${EXTRA_FLAGS}")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -Minfo ${EXTRA_FLAGS}" )
  set(CMAKE_Fortran_FLAGS         ${CMAKE_Fortran_FLAGS_RELEASE})
  set(CMAKE_EXE_LINKER_FLAGS     "${CMAKE_EXE_LINKER_FLAGS}  ${ACC_LINK_FLAGS}")

endif ()

message("==============================================================================")
message("========================>>>>>>>>>>>>>>>>> CMAKE_C_FLAGS=${CMAKE_C_FLAGS}")
message("========================>>>>>>>>>>>>>>>>> CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}")
message("========================>>>>>>>>>>>>>>>>> CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}")
message("========================>>>>>>>>>>>>>>>>> CMAKE_EXE_LINKER_FLAGS=${CMAKE_EXE_LINKER_FLAGS}")
