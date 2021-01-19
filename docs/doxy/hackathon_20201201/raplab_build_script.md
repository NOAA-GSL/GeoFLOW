```bash
#!/bin/bash  

#                                                                                                                             
# This script should be placed inside the top level GeoFLOW directory.                                                        
# The following commands will attempt to ensure a build directory exists                                                      
# within that level directory                                                                                                 
#                                                                                                                             
if [ ! -d "build" ]; then
    echo "build/ directory is not found at current level"
    exit 1
fi

# Load our Modules 
# IMPORTANT:
# You must have configured your system to use the custom modules
# for these paths to appear and work on your system
#
# export MODULEAPPS_ROOT=/home/u00u6ftz35hMh2y7S6357/opt/apps
# export MODULEPATH_ROOT=${MODULEAPPS_ROOT}/modulefiles
# module use ${MODULEPATH_ROOT}

# Using GCC 10.2 + OpenMPI 4.0.5 Compiler
module purge
module load gcc/7.5.0/cmake/3.18.4
module load compilers/gcc-10.2.0
module load gcc/10.2.0/openmpi/4.0.5
module load gcc/10.2.0/openmpi/4.0.5/boost/1.73.0
module load cuda/11.0.2

# Using NVHPC (PGI) Compiler
# module purge
# module load gcc/7.5.0/cmake/3.18.4
# module load compilers/nvhpc-20.9-mpi
# module load nvhpc-mpi/20.9/boost/1.73.0
# module load cuda/11.0.2

# Change into the "out of source" build directory
cd build

# Remove any cached files from previous attempts
rm -r cmake*
rm -r CMake*

# Configure using CMake                                                                                          
cmake -DGEOFLOW_USE_TRACER=ON -DGEOFLOW_TRACER_USE_NVTX=ON --log-level=VERBOSE ../.

# Make it                                                                                                                     
make install
```   