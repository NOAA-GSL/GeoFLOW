
## Input Configuration

### Test Problem
The test case consists of using the element basis functions to approximate
the derivative to a known function (or field) over the grid.  This excersises
the derivative calculation routines which are the primary focus of the hackathon.
All element basis functions and test case parameters are controled through a 
simplified input file that the executable expects to find within the same 
execution directory. 

### Input File: [hack_input.jsn](../../../src/cdg/apps/hack_input.hpp)
The input file **hack_input.jsn** for the hackathon executable is placed by 
the build+install process into the same directory as the compiled binary. 
As the file name extension indicates it is a standard JSON formatted file 
which describes the parameters of the run. The majority of this abbreviated
file can be ignored but a couple parameters will be of interest to hackathon 
participants.

```json
"exp_order"            : [6, 4, 5],
```
The **exp_order** controls the polynomial expansions within each element 
of the grid. In the above example, the first coordinate direction of the 
reference element contains a 6th order polynomial while the second and 
third coordinate directions contain a 4th and 5th order polynomial 
respectively. To properly approximate the derivative and recieve a "Passing"
message at the conclusion of the test run these values must be equal to 
or larger than the **poly** values described in the next section. In 
general, it makes since to have the expansion polynomial equal in all 
directions and greater than the **poly** value for testing (i.e. [4,4,4]).

```json
  "poly_test": {
    "poly"     : [2, 2, 1],
    "ncycles"  : 1000,
    "idir"     : 1
  },
```
The **poly_test** block controls the test case functionality. The known 
function takes the form of f(x,y,z) = x<sup>p</sup> y<sup>q</sup> z<sup>r</sup> 
with the polynomial powers being represented by the 3 values within **poly**. The
**ncylces** variable controls how many times the same derivative calculations
are repeated to provide a level of work that is sufficient to properly profile. 

```json
"num_elems" : [100, 100, 1],
```
The **num_elems** parameter controls the number of elements in each direction
for a cubed region of the domain. This computation isn't part of the kernel
development but can be usefull to control the total number of elements to
perform the test calculations on.

### Quick Reference
| Parameter |    Description   |    Example   |
|-----------|:----------------:|:------------:|
|exp_order  | Expansion polynomial order in reference space for each direction  |   [6, 4, 5]  |
|poly       | Polynomial order of known function we are approximating           |   [2, 2, 2]  |
|ncycles    | Number of repeated derivative calculations to perform             |      1000    |
|idir       | Direction to calculate the derivative in                          | 1 or 2 or 3  |
|num_elems  | Number of elements within a cubed region of the grid              | [50, 50, 50] |

### Example File
```json
{ 
  "exp_order"            : [4, 4, 4],
  "grid_type"            : "grid_box",
  "IO_implementation"    : "gio",
  "poly_test": {
    "poly"     : [2, 2, 1],
    "ncycles"  : 1000,
    "idir"     : 1
  },
  "fully_periodic": {
    "bdy_class"        : "uniform",
    "base_type"        : ["GBDY_PERIODIC"],
    "istate"           : [[0, 1, 2]]
  },
  "grid_box": {
    "grid_name" : "grid_box",
    "xyz0"      : [0.0, 0.0, 0.0],
    "delxyz"    : [1.0, 1.0, 1.0],
    "num_elems" : [100, 100, 1],
    "bdy_x_0"   : "fully_periodic",
    "bdy_x_1"   : "fully_periodic",
    "bdy_y_0"   : "fully_periodic",
    "bdy_y_1"   : "fully_periodic",
    "bdy_z_0"   : "fully_periodic",
    "bdy_z_1"   : "fully_periodic"
  },
  "gio": {
    "ivers"   : 0,
    "multivar" : false,
    "io_type"  : "collective",
    "wtime"    : 6,
    "wtask"    : 5,
    "wfile"    : 2048
  }

}
```
