
## Compile and Run

### Download GitHub Repository
Your first step is to get the code from the GitHub reposity
```console
User@machine: > git clone https://github.com/NOAA-GSD/GeoFLOW.git
```

### Branch feature/hackathon
For the hackathon the code has been reduced to a small kernel executing 
the portion of code we plan to concentrate on. This code (and the respective 
inputs) are only contained within the **feature/hackathon** branch of the GitHub 
repository so make sure you are inside the correct branch before compiling.

```console
User@machine: > cd GeoFLOW
User@machine: > git checkout feature/hackathon
```

### Compiling
Compiling the code is identical to the instructions provided on the main 
[README](../../README.md). A pre-compiled Boost Library will/has be/been 
provided on the machines we will work on to eliminate that hurdle and it's 
recommended we all use the matching compiler those are built with.  

### Provided Boost Library (on Hera)
```console
export BOOST_ROOT=/scratch2/BMC/gsd-hpcs/Bryan.Flynt/opt/boost/intel-2020.2
```

### Provided Boost Library (on Raplab)
To use GCC 10.2 with OpenMPI 4.0.5 the Boost Library is
```console
export BOOST_ROOT=/home/u00u6ftz35hMh2y7S6357/opt/apps/gcc-10_2_0/openmpi-4_0_5/boost/1.73.0
```

To use NVHPC 20.9 the Boost Library is
```console
export BOOST_ROOT=/home/u00u6ftz35hMh2y7S6357/opt/apps/nvhpc-20_9_mpi/boost/1.73.0
```

### Turning On/Off NVidia Profiler (NVTX)
The hackathon code has the ability to be instrumented using the 
[NVIDIA Tools Extension SDK (NVTX)](https://developer.nvidia.com/blog/cuda-pro-tip-generate-custom-application-profile-timelines-nvtx/)
profiling calls inserted automatically by the CMake build system. To accomplish 
this the Tracer must be turned on using the **GEOFLOW_USE_TRACER** variable
followed by the **GEOFLOW_TRACER_USE_NVTX** variable to specify what the tracer 
should do.  Turning on these two CMake variables the build system will detect 
the CUDA library, insert the calls to NVTX, compile and link the proper libraries.

## Complete Example (on Hera)

### Clone the repository and create a **build** directory.

```console
User@machine: > git clone https://github.com/NOAA-GSD/GeoFLOW.git
User@machine: > cd GeoFLOW
User@machine: GeoFLOW> git checkout feature/hackathon
User@machine: GeoFLOW> mkdir build
User@machine: GeoFLOW> ls
analysis  build  cmake  CMakeLists.txt  DISCLAIMER  docs  extern  LICENCE.md  README.md  scripts  src  test
```

### Build the code using the proper environment setup.  

To make it easier, I typically use the following script which I run from within 
the top level GeoFLOW directory. For your first attempt, I'd recommend to copy 
and past each line to see messages if any errors result.

[Hera Build Script](hera_build_script.md)

[Raplab Build Script](raplab_build_script.md)


### Running the test stub

```console
User@machine: GeoFLOW> cd build/bin
User@machine: GeoFLOW/build/bin> ls
hackathon_driver  hack_input.jsn
User@machine: GeoFLOW/build/bin> module load intel/2020.2 impi/2020.2
User@machine: GeoFLOW/build/bin> ./hackathon_driver
main: ---------------------------derivative OK: 2.01001e-14 : direction=1 method: old
main: ---------------------------derivative OK: 2.01001e-14 : direction=1 method: new
```

