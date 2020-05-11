# GeoFLOW
Geo FLuid Object Workbench

## Compilation from Source on UNIX
To compile the source distribution, you need at least the following to build the executable:
* [CMake](https://cmake.org/) version 3.1.0 or later to generate the Makefile for your platform 
* C++ compiler, supporting the C++11 standard or later
    * OpenMP directives (optional)
* MPI Library compatible with your C++ compiler    
* [Boost C++](https://www.boost.org/) headers with the mpi and serialization libraries compiled
    * This requirment may be removed in the future


## Compilation is done by performing the following steps

1. Download the repository:
```console
git clone https://github.com/NOAA-GSD/GeoFLOW.git
```

2. Create a build directory (for instance inside the source tree): 
```console
cd GeoFLOW
mkdir build
cd build
```

3. Configure your system to have the proper libraries visable to the CMake build system:  
    - Boost 
        - CMake will check BOOST_ROOT then system paths for Boost libraries
    ```console
	export BOOST_ROOT=/my/dir/to/boost/1.71.0                     # Boost Installation
	```
	- MPI 
	    - CMake will check path of *mpiexec* for compilers then envirnment variables
	```console 
    export MPI_C_COMPILER=/my/path/to/mpi/wrapped/mpicc           # MPI Wrapped C Compiler
    export MPI_CXX_COMPILER=/my/path/to/mpi/wrapped/mpicxx        # MPI Wrapped C++ Compiler
    export MPI_Fortran_COMPILER=/my/path/to/mpi/wrapped/mpifort   # MPI Wrapped F Compiler
    ```
**Note:** These are environment variables used by CMake and not part of the GeoFLOW application

4. Run cmake to generate the Makefile:
```console
cmake ..
```
cmake tries to determine the platform you use, and will look for the requires tools. It will report if something is missing.

5. Compile the program by running make:
```console
make
```

6. Optional: Run tests cases.
```console
make test
```

7. Optional: Generate the documentation. 
```console
make docs
```






