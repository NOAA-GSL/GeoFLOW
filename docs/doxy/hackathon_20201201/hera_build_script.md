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
module purge
module load cmake
module load intel/2020.2
module load impi/2020.2
module load cuda

# Set environment variables for CMake to pick up
export CC=mpiicc
export CXX=mpiicpc
export FC=mpiifort
export BOOST_ROOT=/scratch2/BMC/gsd-hpcs/Bryan.Flynt/opt/boost/intel-2020.2

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