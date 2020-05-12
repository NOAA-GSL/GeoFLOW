# GeoFLOW API Documentation {#mainpage}

GeoFLOW (Geo FLuid Object Workbench) is a research-oriented software development effort in 
the High Performance Computing Branch of the Global Systems Lab (GSL) at NOAA. The project
has 57 main objectives.  First, GeoFLOW enables the application of spectral element methods to
the solution of partial differential equations relavent to geophysical sciences.  yadda, yadda...

The GeoFLOW framework source code is partitioned into a number of software "packages".  Each 
package contains a collection of classes that are logically-related with respect to the role they 
play in spectral element applications. Links to the package contents follow the list of package 
descriptions:

* [tbox](@ref ToolBox)  
    The "Tool Box" package which provides basic utility classes that are used throught the application 
    development 

* [pdeint](@ref PDE_Int)  
    The PDE Integrators package provided templated classes that are used to integrate systems 
    of equations forward in time
    
* [kitchen sink](@ref Kitchen)  
    The "Kitchen Sink" package which provides everything else 
