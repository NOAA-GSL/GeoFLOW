# GeoFLOW
Geo FLuid Object Workbench

## Compilation from Source on UNIX
To compile the source distribution, you need at least the following to build the executable:
* [CMake](https://cmake.org/) version 3.1.0 or later to generate the Makefile for your platform. 
* C++ compiler, supporting the C++11 standard or later.
    * OpenMP directives (optional).
* [Boost C++](https://www.boost.org/) headers with the mpi and serialization libraries compiled
    * This requirment may be removed in the future
* MPI Library compatible with your C++ compiler

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

3. Run cmake to generate the Makefile:
```console
cmake ..
```
cmake tries to determine the platform you use, and will look for the requires tools. It will report if something is missing.

4. Compile the program by running make:
```console
make
```

5. Optional: Run tests cases.
```console
make test
```

6. Optional: Generate the documentation. 
```console
make docs
```






