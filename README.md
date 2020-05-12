# GeoFLOW
Geo FLuid Object Workbench

## Compilation from Source on UNIX

To compile the source distribution, you need at least the following to build the executables:
* [CMake](https://cmake.org/) version 3.9.0 or newer to generate the Makefile for your platform 
* C++ compiler, supporting the C++11 standard or newer
    * Support for OpenMP (optional)
* MPI library compatible with your C++ compiler    
* [Boost C++](https://www.boost.org/) library with the *mpi* and *serialization* libraries compiled
    * This requirment may be removed in the future

Optional tasks can be performed using the following additional applications:
* [DOxygen](http://www.doxygen.nl/) to generate API documentation
* [ParaView](https://www.paraview.org/) to visualize grid and solution files
* [Python 3](https://www.python.org/) for post processing grid and solution files


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
**Note:** These rules and variables used by CMake and not part of the GeoFLOW application
    - Boost 
        - CMake will check BOOST_ROOT then system paths for Boost libraries
    ```console
	export BOOST_ROOT=/my/dir/to/boost/1.71.0                     # Boost Installation
	```
	- MPI 
	    - CMake will check path of *mpiexec* for compilers then environment variables
	```console 
    export MPI_C_COMPILER=/my/path/to/mpi/wrapped/mpicc           # MPI Wrapped C Compiler
    export MPI_CXX_COMPILER=/my/path/to/mpi/wrapped/mpicxx        # MPI Wrapped C++ Compiler
    export MPI_Fortran_COMPILER=/my/path/to/mpi/wrapped/mpifort   # MPI Wrapped F Compiler
    ```

4. Run cmake to generate the Makefile:
```console
cmake ..
```
CMake tries to determine the platform you use, and will look for the required tools. It will report if something is missing.

5. Compile and install the program by running make:
```console
make install
```

6. Optional: Run test cases.
```console
make test
```

7. Optional: Generate the documentation. 
```console
make docs
```
To view the detailed API documentation point a web browser at
"docs/html/index.html" within the *build* directory.  
**Note:** The API documentation is intended for developers who 
want to see details of the code inside the GeoFLOW framework.
It is not intended for GeoFLOW users who simply want to run an 
application.  


## Examples and Program Usage

For examples and instructions on how to run GeoFLOW see our 
examples page [here](docs/doxy/examples.md)





