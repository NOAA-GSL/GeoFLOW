#!/bin/bash
#
# Script to build Boost library on Travis platform
#
# ?? Why not just build it through the .travis.yml file ??
# Because the YAML markup language will not allow :
# to be used since they are special character within
# the language.  So we call this script instead.
#

wget https://dl.bintray.com/boostorg/release/1.72.0/source/boost_1_72_0.tar.gz
tar -xzf boost_1_72_0.tar.gz #> /dev/null 2>&1
cd boost_1_72_0
export TOOLSET=gcc
export MPI_WRAPPER=mpicxx
export NPROC=4
export INSTALL_LOC=`pwd`/gcc
export BUILD_LIBS=mpi,serialization
./bootstrap.sh --with-toolset=${TOOLSET} --with-libraries=${BUILD_LIBS} --prefix=${INSTALL_LOC}  #> /dev/null 2>&1
printf "using mpi : %s ;\n" "${MPI_WRAPPER}" >> project-config.jam
./b2 -j ${NPROC} install #> /dev/null 2>&1

