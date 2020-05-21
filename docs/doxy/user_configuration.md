# GeoFLOW User Configuration: How to set up a problem and interact with the system

# I. Introduction.

 The primary point of interaction with GeoFLOW is the user-configurable JSON file.
 The default file name is 'input.jsn', and it is assumed to lie in the run 
 (working) directory. On the command line, this default config file name may be changed 
 by using the "-i" command line argument:

> ./geoflow_cg -i my_new_config_file.jsn


 Inside this config text file, there are several items that must be set in order
 to run a problem. User must at a minumum: <br>
A. Specify a solver name that tells us which initial value-boundary value 
    (IVBV) partial differential equation (PDE) set to solve; <br>
B. Specify a time-stepping scheme <br>
C. Specify which grid to use, and probably also the basis expansion order <br>
D. Specify a set of boundary conditions <br>
E. Specify a set of initial conditions <br>
F. Specify the evolution (integration) time period <br>
G. Configure output: specify which variables to output, and their cadences. <br>

The user may also wish to: 

H. Specify external forcing <br>
I. Specify grid terrain. <br>

Each of these specifications is described here in turn below. One additional
quantity that is not user-configurble at run-time is the problem dimension, GDIM,
which is specified by defining the preprocessor variables _G_IS2D, or _G_IS3D 
at build time. This is described in XXX. We simply refer to this as necessary
below.

# II. Config file example.

Before getting started, here is a sample configuration file. Let's assume it's 
given the name "input.jsn" if we need that later. We will refer to this as 
we proceed, and the user will have a concrete example to serve as a template.
For the time being, we assume that we are configuring only GeoFLOW Spectral Element 
discretizations.

```json
{
  "pde_name"             : "pde_burgers",
  "exp_order"            : [4, 4, 4],
  "exp_order_type"       : "constant",
  "grid_type"            : "grid_box",
  "initstate_block"      : "initstate_icosnwave",
  "initforce_block"      : "",
  "use_forcing"          : false,
  "initstate_type"       : "direct",
  "initforce_type"       : "",
  "restart_index"        : 0,
  "benchmark"            : false,
  "do_comparison"        : true,
  "observer_list"        : ["binary_observer", "diag_observer"],
  "IO_implementation"    : "gio",
  "terrain_type"         : "myboxterrain",
  "myboxterrain": {
    "name"    : "boxgauss_range",
    "x0"      : [0.5],
    "y0"      : [0.0],
    "xsigma"  : [0.1],
    "ysigma"  : [0.1],
    "h0"      : [0.15]
  },
  "initvabc": {
    "name"  : "abc",
    "kdn"   : 2,
    "kup"   : 3,
    "kpower": 0.5,
    "A"     : 0.9,
    "B"     : 1.0,
    "C"     : 1.1,
    "E0"    : 1.0
  },
  "compinit": {
    "initv"     : "initvabc",
    "initdt"    : "random",
    "inittemp"  : "zero",
  },
  "initstate_icosgaussburgers": {
    "latitude0"  : [0.0, -40.0],
    "longitude0" : [-90.0,-90.0],
    "sigma"      : [0.25, 0.25],
    "c0"         : [1.0,1.0],
    "alpha"      : 0.0,
    "u0"         : 1.0,
    "nlumps"     : 1
  },
  "initstate_boxdirgauss": {
    "x0"     : 0.5,
    "y0"     : 0.5,
    "z0"     : 0.5,
    "sigma"  : 0.05,
    "u0"     : 1.0,
    "adv_vel": [1.0, 0.0, 0.0]
  },
  "initstate_boxnwave": {
    "x0"        : [0.35 , 0.65],
    "y0"        : [0.35 , 0.65],
    "z0"        : [0.35 , 0.65],
    "ULparam"   : [1.0  , 1.0],
    "t0"        : [0.002, 0.008],
    "planar"    : [false, false],
    "prop_dir_x": [0    , 0],
    "prop_dir_y": [0    , 0],
    "prop_dir_z": [0    , 0]
  },
  "initstate_icosnwave": {
    "latitude0"  : [10,  40 ],
    "longitude0" : [280, 330],
    "Uparam"     : [1.0, 1.0],
    "t0"         : [0.08, 0.16]
  },
  "pde_burgers": {
    "doheat"      : false,
    "bpureadv"    : false,
    "bconserved"  : false,
    "forcing_comp": [0, 1, 2]
  },
  "dirichlet_fixed": {
    "bdy_class"        : "uniform",
    "base_type"        : "GBDY_DIRICHLET",
    "value"            : [0.0, 0.0, 0.0],
    "bdy_config_method": ""
  },
  "grid_icos": {
    "grid_name"   : "grid_icos",
    "ilevel"      : 2,
    "refine_type" : "GICOS_BISECTION",
    "radius"      : 1.0,
    "maxit"       : 128,
    "tol"         : 1.0e-8,
    "norm_type"   : "GCG_NORM_INF"
  },
  "grid_sphere" : {
    "grid_name"         : "grid_sphere",
    "radiusi"           : 1.0,
    "radiuso"           : 2.0,
    "num_radial_elems"  : 10,
    "bdy_inner"         : "dirichlet_fixed",
    "bdy_outer"         : "sponge_layer",
    "use_state_init"    : false
    "bdy_init_method  " : "",
    "bdy_update_method" : "",
  },
  "fully_periodic": {
    "base_type"        : "GBDY_PERIODIC"
    "bdy_class"        : "uniform"
    "bdy_config_method": ""
  },
  "compound_inflow_example": {
    "bdy_class"        : "user",
    "base_type"        : "GBDY_DIRICHLET",
    "bdy_config_method": "circular_jet"
  },
  "burgers_inflow": {
    "bdy_class"        : "uniform",
    "base_type"        : "GBDY_INFLOWT",
    "bdy_config_method": ""
  },
  "grid_box": {
    "grid_name" : "grid_box",
    "xyz0"      : [0.0, 0.0, 0.0],
    "delxyz"    : [1.0, 1.0, 1.0],
    "num_elems" : [16, 16, 1],
    "bdy_x_0"   : "fully_periodic",
    "bdy_x_1"   : "fully_periodic",
    "bdy_y_0"   : "fully_periodic",
    "bdy_y_1"   : "fully_periodic",
    "bdy_z_0"   : "fully_periodic",
    "bdy_z_1"   : "fully_periodic",
    "use_state_init"    : true,
    "bdy_init_method  " : "",
    "bdy_update_method" : "",
    "maxit"     : 128,
    "tol"       : 1.0e-8,
    "norm_type" : "GCG_NORM_INF"
  },
  "stepper_props": {
    "stepping_method"  : "GSTEPPER_EXRK",
    "time_deriv_order" : 4,
    "extrap_order"     : 2,
    "variable_dt"      : false,
    "courant"          : 0.04
  },
  "gio": {
    "ivers"   : 0,
    "multivar" : false,
    "io_type"  : "collective",
    "wtime"    : 6,
    "wtask"    : 5,
    "wfile"    : 2048
  },
  "curlv": {
    "state_index" : [0, 1, 2],
    "names"       : ["omegax", "omegay", "omegaz"],
    "mathop"      : "curl"
  },
  "curlvmag": {
    "state_index" : [0, 1, 2],
    "names"       : ["omega"],
    "mathop"      : "curlmag"
  },
  "vmag": {
    "state_index" : [0, 1, 2],
    "names"       : ["vmag"],
    "mathop"      : "vmag"
  },
  "gradv": {
    "state_index" : [0],
    "names"       : ["gradx", "grady", "gradz"],
    "mathop"      : "grad"
  },
  "binary_observer": {
    "observer_name"      : "gio_observer",
    "indirectory"        : "outs",
    "outdirectory"       : "outs",
    "cadence_type"       : "cycle",
    "cycle_interval"     : 800,
    "time_interval"      : 0.01,
    "state_index"        : [0, 1],
    "state_names"        : ["u1", "u2", "u3"],
    "grid_names"         : ["xgrid", "ygrid", "zgrid"],
    "agg_state_name"     : "u",
    "agg_grid_name "     : "grid",
    "derived_quantities" : ["vmag", "curlvmag"],
    "filename_size"      : 2048
  },
  "diag_observer": {
    "observer_name"      : "global_diag_basic",
    "indirectory"        : ".",
    "outdirectory"       : ".",
    "interval_freq_fact" : 10.0,
    "treat_as_1d"        : false
  },
  "dissipation_traits": {
    "nu"      : 1.0e-1,
    "nu_type" : "constant"
  },
  "time_integration": {
    "integ_type" : "cycle",
    "time_end"   : 0.150,
    "cycle_end"  : 1,
    "dt"         : 2.0e-6
  }
} 
```

# III. Configuration steps.

## A. Specify a solver name.

  To specify a PDE solver, simply provide it in the input.jsn file by 
setting:

```json
    "pde_name" : "pde_burgers"
```

  Currently, only "pde_burgers" is enabled, which refers to a linear or nonlinear 
  advection-diffusion solver. In the current configuration, there must be a JSON 
  block in the file with this same name that will configure the solver. If you 
  look below in input.jsn, you'll see this block in the file. This block 
  tells us that the "pde_burgers" can also solve the heat equation (with a 
  single state variable), or the linear advection/diffusion. If these are both 
  set to False, then the default nonlinear multi-dimensional Burgers equation is
  solved. The flag "bconserved", not currently available, tells us to use a 
  conserved form of the equations, if applicable. The final configuration item,
  "forcing_comp" tells us which components to apply forcing to, if we use forcing. 
  This will be discussed later. Other solvers may require additional configuration.

## B. Specify a time-stepping scheme

  In the current configuration, the time stepping configuration is not "attached"
  to the PDE solver; it a future version it likely will be, and it will take an
  arbitrary name. For now, this block is called "stepper_props", and contains the 
  following configurable items
    	
```json
    "stepper_props": {
      "stepping_method"  : "GSTEPPER_EXRK",
      "time_deriv_order" : 4,
      "extrap_order"     : 2,
      "variable_dt"      : false,
      "courant"          : 0.04
    },
```

  where the items are described here:

  |expression         |  description                                      |
  |------------------ |---------------------------------------------------|
  |"stepping_method"  |May be one of the valid entries in src/cdg/include/gtypes.h: in the sGStepperType array. At present, only explicit Runge-Kutta integration is enabled.|
  |"time_deriv_order" | refers to the time truncation order of the scheme (for all available orders. The current scheme allows up to order 5.| 
  |"extrap_order"     | refers to the extrapolation order, if this method were to be enabled.|
  |"variable_dt"      | tells the stepper that the timestep may vary from step to step. If true, a method is called that limits the size of the step based on the current state.|
  |"courant"          | is the Courant number: a "fudge" factor that multiplies the timestep from the timestep method when "variable_dt" = true, so that the Courant condition isn't violated|

## C. Select and configure grid.

  The grid type is specified at the top of the JSON file by setting the 
  "grid_type" variable. Like with PDE specification, this is a quantity that
  will also name a self-referential configation block within the file. 

```json
    "grid_type"   : "grid_icos",
```

  Valid "grid_type"'s are: "grid_icos", or "grid_box".  

  If GDIM=2, then "grid_box" will refer to a planar Cartesian grid; if GDIM=3, 
  it will refer to a 3d Cartesian box grid. 

  If GDIM=2, the "grid_icos" refers to a 2d spherical surface grid, constructed from
  a base icosohedron; if GDIM=3, "grid_icos" uses the 2d "grid_icos" grid at the 
  surface, and constructs elements radially from this "base" 2d grid.

### (i) Configuring the box grids:

  The self-referential config block for the "grid_box" grids  look like:

```json
      "grid_box": {
        "grid_name" : "grid_box",
        "xyz0"      : [0.0, 0.0, 0.0],
        "delxyz"    : [1.0, 1.0, 1.0],
        "num_elems" : [16, 16, 1],
        "bdy_x_0"   : "fully_periodic",
        "bdy_x_1"   : "fully_periodic",
        "bdy_y_0"   : "fully_periodic",
        "bdy_y_1"   : "fully_periodic",
        "bdy_z_0"   : "fully_periodic",
        "bdy_z_1"   : "fully_periodic",
        "bdy_init_method  " : "",
        "bdy_update_method" : "",
        "use_state_init"    : true,
        "maxit"     : 128,
        "tol"       : 1.0e-8,
        "norm_type" : "GCG_NORM_INF"
      },
```

  where entries and their descriptions are:

  |expression           |  description                                      |
  |---------------------|---------------------------------------------------|
  |"grid_name"          | name of grid from src/cdg/grid/ggrid_factory.ipp ; e.g., "grid_box"|       
  |"xyz0"               | tuple of startig lower left of box. Must be of size >= GDIM.|
  |"delxyz"             | tuple of widths in each direction. Must be of size >= GDIM.|
  |"num_elems"          | number of elements in each direction. Must be of size >= GDIM.|
  |"bdy_x_0"            | name of config block defining boundary conditions on west face|
  |"bdy_x_1"            | name of config block defining boundary conditions on east face|
  |"bdy_y_0"            | name of config block defining boundary conditions on south face|
  |"bdy_y_1"            | name of config block defining boundary conditions on north face|
  |"bdy_z_0"            | name of config block defining boundary conditions on bottom face.  Used if GDIM=3|
  |"bdy_z_1"            | name of config block defining boundary conditions on top face Used if GDIM=3|
  |"bdy_init_method"    | if specified, determines how boundaries are initialized|
  |"bdy_update_method"  | if specified, determines how boundaries are updated|
  |"use_state_init"     | If specified, the boundary initialization for DIRICHLET-type conditions will be taken from the initial conditions|
  |"maxit"              | max no. Krylov iterations if using terrain|
  |"tol"                | Krylov loop tolerance is using terrain|
  |"norm_type"          | norm to use in establishing Krylov loop residual. Valid values are provided in src/pdeint.in_solver_base.hpp: "GCG_NORM_INF", "GCG_NORM_EUC", "GCG_NORM_L2", "GCG_NORM_L1". Used if doing terrain.|
  

  For the GeoFLOW Spectrarl Element-like discretizations, there is one "grid-like"
  quantity that is important for the user to be able set. This is the polynomial
  expansion order within each direction within each element. This is currently
  specified outside of the "grid_type" config block. It is set in the config file example
  with the JSON parameters:

  ```json
  "exp_order"            : [4, 4, 4],
  "exp_order_type"       : "constant",
  ```
  At this time, the only "exp_order_type" available is "constant", and this means that
  the expansion order for each direction is the same for all elements. In the future, 
  it's possible we will loosen this restriction. The "exp_order" vector must be of 
  size at least GDIM in the JSON file, and the expansion order may be specified 
  independently in each direction.


## D. Boundary condition specification.

  In the "grid_box" configuration, each enumerated natural boundary is configured 
  to use a "fully_peridoc" boundary condition in each direction via a statement like:
      
```json
  "bdy_x_0"           : "fully_periodic"
  "bdy_x_1"           : "fully_periodic"
```

  The name on the right specifies a self-referenced configuration block in the 
  JSON file that actually configures that boundary. It's important to remember 
  that these blocks govern _specification_ of the boundary, not, e.g., boundary 
  _initialization_. Boundary initialization and evolution is considered below. 

  These boundary spec blocks, in general, have the form:

```json
       "fully_periodic": {
         "base_type"        : "GBDY_PERIODIC"
         "bdy_class"        : "uniform"
         "bdy_config_method": ""
       },
```

  Valid "base_type" values are provided in 

     src/cdg/include/gtypes.h

  in the the sGBdyType array. Currently, BDY_DIRICHLET (or GBDY_INFLOWT, 
  these are the same thing), and GBDY_PERIODIC work currently. The types GBDY_0FLUX 
  and GBDY_SPONGE will be implemented soon. 

  Valid "bdy_class" expressions may be: "uniform" or "user". 
     
  | expression | description                                     |
  |------------|-------------------------------------------------|
  |"uniform"   | the condition is applied automatically to the entire natural boundary|
  |"user"      | if specified the user MUST specify the "bdy_config_method".|
       
  The named "bdy_config_method" is used to set the boundary condition type. 
  Methods that may be provided are specified in the factory associated with 
  boundary specification: 
    src/cdg/init/gspecbdy_factory.cpp. 
  The steps by which to add a new boundary specification method are discussed later.

  "GBDY_PERIODIC" boundary conditions are somewhat special. They may be specified
  for only "grid_box" boundaries, and, if applied, must be applied to the entire 
  boundary; i.e., it will require that the "bdy_class" is "uniform". Also, they must 
  be applied pairwise between associated boundaries. So the following would be 
  incorrect:

```json
  "bdy_x_0"   : "fully_periodic",
  "bdy_x_1"   : "burgers_inflow".
```

  Note that if the user wants "bdy_class" to be "uniform", since this is the default,
  neither the "bdy_class" nor the "bdy_config_method" need to be provided in the
  boundary spec block; these two parameters need be specified only if "bdy_class" 
  is set to "user". 

  Note that a "user" specified spec method may require other parameters
  to be specified in the spec block.
       

### Boundary condition initialization and updating. 

  If a boundary type is GBDY_DIRICHLET/GBDY_INFLOWT or GBDY_SPONGE (when 
  ready), the boundary values must be initialized and evolved. Here we describe
  how this is done.
      
  Two variables in the "grid_box" config block determine how the boundaries are
  initialized:

```json
       "use_state_init"    : true,
       "bdy_init_method  " : "",
       "bdy_update_method" : "",
```

  If "use_state_init" = true, then the next two parameters are not required.
  This options suggests that the boundaries are initialized (and evolved)
  by using the initial (and evolving) state vector.

  If "use_state_init" = false, then the boundary initialization method must be
  provided by setting the "bdy_init_method".  The valid methods can be found in 

    src/cdg/init/ginitbdy_factory.ipp

  If "use_state_init" = false, then the boundary update method should also be
  provided by setting the "bdy_update_method".  The valid methods can be found in 

    src/cdg/init/gupdatebdy_factory.ipp

  We discuss later how to write new boundary initial and update methods and 
  make them available in the system.

### (ii) Configuring spherical grids.

  When GDIM=2, we use the following configuration block for a spherical
  surface:

```json
      "grid_icos": {
        "grid_name"         : "grid_icos",
        "radius"            : 1.0,
        "ilevel"            : 2,
        "refine_type"       : "GICOS_BISECTION",
      },
```

 where

  | expression         | description                                     |
  |--------------------|-------------------------------------------------|
  |"grid_name"         |"grid_icos"|
  |"radius"            |scalar radius|
  |"refine_type"       |type of refinement: either "GICOS_BISECTION"" or "GICOS_LAGRANGIAN" |
  |"ilevel"            |if "refine_type" = "GICOS_BISECTION", then this is the number of bisection levels, yielding a grid with 60*(2^ilevel) elements.  If "refine_type" = "GICOS_LAGRANGIAN", ilevel gives the number of divisions along a triangle edge, yielding 60*(ilevel+1)^2 elements on the grid.|

  Note: there are no boundaries on the 2d sphere, so no boundary conditions need to
  be specified.


  When GDIM=3, we use the following configuration block for a spherical
  grid:

```json
      "grid_sphere" : {
        "grid_name"         : "grid_sphere",
        "radiusi"           : 1.0,
        "radiuso"           : 2.0,
        "num_radial_elems"  : 10,
        "refine_type"       : "GICOS_BISECTION",
        "ilevel"            : 2,
        "num_radial_elems"  : 10,
        "bdy_inner"         : "dirichlet_fixed",
        "bdy_outer"         : "sponge_layer",
        "use_state_init"    : false
        "bdy_init_method  " : "",
        "bdy_update_method" : "",
      },
```

  where entries and their descriptions are:
       
  | expression          | description                                     |
  |---------------------|-------------------------------------------------|
  |"grid_name"          | "grid_icos"|
  |"radiusi"            | scalar inner radius|
  |"radiuso"            | scalar outer radius (not used if GDIM=2)|
  |"refine_type"        | type of refinement: either "GICOS_BISECTION"" or "GICOS_LAGRANGIAN" 
  |"ilevel"             | if "refine_type" = "GICOS_BISECTION", then this is the number of bisection levels, yielding a grid with 60*(2^ilevel) elements.  If "refine_type" = "GICOS_LAGRANGIAN", ilevel gives the number of divisions along a triangle edge, yielding 60*(ilevel+1)^2 elements on the grid.|
  |"num_radial_elems"   | number of elements in radial direction|
  |"bdy_inner"          | name of config block defining boundary conditions on inner surface|
  |"bdy_outer"          | name of config block defining boundary conditions on outersurface|
  |"bdy_init_method"    | if specified, determines how boundaries are initialized|
  |"bdy_update_method"  | if specified, determines how boundaries are updated|
  |"use_state_init"     | If specified, the boundary initialization for DIRICHLET type conditions will be taken from the initial conditions|
      
  Note: The 3d spherical grid is built upon the 2d surface grid by extending elements
  radially. The boundary condition spec, initilization, and update methods work
  the same was as for box grids.

E. Specify state initial conditions.

  In the sample config file in Sec. II, there are two primary statements that control
  state initialization:

```json
  "initstate_block"      : "initstate_icosnwave",
  "initstate_type"       : "direct",
```
  where

  | expression          | description                                     |
  |---------------------|-------------------------------------------------|
  | "initstate_block"   | specifes JSON block that configures initial conditions|
  | "initstate_type"    | if "direct", then named block will specify a method available in src/cdg/init/ginitstate_factory.ipp that initializes entire state vector in one call; if "component", then named block will specify named config blocks that initializes state vector component- by component|

  In the sample JSON file, we assume we are calling the method "initstate_icosnwave"
  in the ginitstate_factory.ipp, and this method is of type "direct", and so
  sets the entire state vector, which may consist of velocity, pressure, density
  temperature, ... components. The named direct state initialization block
  has the form

  ```json 
  "initstate_icosnwave": {
    "latitude0"  : [10,  40 ],
    "longitude0" : [280, 330],
    "Uparam"     : [1.0, 1.0],
    "t0"         : [0.08, 0.16]
  },
  ```
  This block not only names the method in the factory, but it also provides
  tunable parameters for the method. The "direct" methods are a good way to
  set an entire state with ``canned'' but tunable idealized initial data.

  On the other hand, assume we have a PDE solver that accommodates a complicated state 
  consisting of velocity, temperature and density components that we want to 
  initialize independently. Then we might set 

```json
  "initstate_block"      : "compinit",
  "initstate_type"       : "component",
  "initvabc": {
    "name"  : "abc",
    "kdn"   : 2,
    "kup"   : 3,
    "kpower": 0.5,
    "A"     : 0.9,
    "B"     : 1.0,
    "C"     : 1.1,
    "E0"    : 1.0
  },
  "compinit": {
    "initv"     : "initvabc",
    "initdt"    : "random",
    "inittemp"  : "zero",
  },

``
  which tells us to use the named block, "initvabc" to set the velocity components, 
  and to initialize the total density, "initdt", using the built-in scheme to set 
  random values, and, finally,  to set the temperature, "inittemp" using the built-in zeroing
  method. The "initvabc" block tells us to use the method named "abc" to set the velocity
  components for an Arnoldi, Beltrami, Childress ABC) flow,  and this method is defined 
  in src/cdg/init/ginitstate_fatory.ipp. If we had, say a
  magnetic field vector as well, we could also initialize it to use the the same "abc" method,
  if we wanted to or we could use a different vector initialization scheme. We reserve for 
  another chapter a discussion about how a user may construct and ``wire in'' a new initialization
  scheme.

  The currently available JSON names for corresponding component types, 
  GStateCompType in the src/cdg/gtypes.h file are provided here:

  | JSON componet name  | GStateCompType   | description
  |---------------------|------------------|--------------------------|
  | "initv"             | GSC_KINETIC      | velocity components      |
  | "initb"             | GSC_MAGNETIC     | magnetic field components|
  | "initdt"            | GSC_DENSITYT     | total density            |
  | "initd1"            | GSC_DENSITY1     | 1-density component      |
  | "initd2"            | GSC_DENSITY2     | 2-density component      |
  | "inittemp"          | GSC_TEMPERATURE  | temperature  component   |
  | "initc"             | GSC_PRESCRIBED   | prescribed component (e.g.,
                                             advection velocity components)   |


F. Specify evolution time.

  The user must be able to specify a number of time cycles for which to 
  integrate the PDE. This is done in the JSON file with the following block:

  ```json
  "time_integration": {
    "integ_type" : "cycle",
    "time_end"   : 0.150,
    "cycle_end"  : 1,
    "dt"         : 2.0e-6
  }
  ```
  Valid values for "integ_type" are "cycle" or "time", which means that the
  iterator will be based on the time _cycle_ or the evolutionary _time_. 
  So, if "integ_tupe" is "cycle", then the itegration will proceed from the 
  current cycle to the "cycle_end" value. If "integ_type" is "time", then
  the integration will proceed from the current value to "time_end". The time
  step size, "dt", is always dimensioned.

  A related issue is the restart. Each binary output file  contains a complete
  snapshot of the PDE state that also encodes the timestamp and cycle number
  at which the dump was made. The time or cycle cadence of these output files 
  is described in Sec. G, but we point out that upon restart, the initial time
  cycle and time stamp are taken from the restart file. Each restart file is
  created with a time index, so to restart from one of these indexed files, the
  user must simply specify the restart index:

  ```json
  "restart_index"        : 50,
  ```
  in the JSON input file. If the "restart_index" is 0, this always starts
  a new run.

G. Configure output and ``Observers''.

  Output is necessary if we want to analyze PDE solutions, to monitor a run, or
  if we simply want to be able to restart a solution to continue an integration.
  In GeoFLOW, each type of output is considered to be an ``observer'', which
  refers to something that peeks at the solution while it's being integrated.

  Currently, there are two classes of observers: binary output or restart observers,
  and diagnostic observers, and there is one implementation of each. Configuration
  of each of these is done in the sample config file by way of the following statements:

  ```json
  "observer_list"        : ["binary_observer", "diag_observer"],
  "binary_observer": {
    "observer_name"      : "gio_observer",
    "indirectory"        : "outs",
    "outdirectory"       : "outs",
    "cadence_type"       : "cycle",
    "cycle_interval"     : 800,
    "time_interval"      : 0.01,
    "state_index"        : [0, 1],
    "state_names"        : ["u1", "u2", "u3"],
    "grid_names"         : ["xgrid", "ygrid", "zgrid"],
    "derived_quantities" : ["vmag", "curlvmag"],
    "agg_state_name"     : "u",
    "agg_grid_name "     : "grid",
    "filename_size"      : 2048
  },
  "diag_observer": {
    "observer_name"      : "global_diag_basic",
    "indirectory"        : ".",
    "outdirectory"       : ".",
    "interval_freq_fact" : 10.0,
    "treat_as_1d"        : false
  },
  ```
  The "observer_list" simply contains a vector of the configuration blocks, which 
  carry out the configuration of the observer, src/pdeint/observer_factory.ipp.

  "binary_observer" is an observer that handles binary (restart) output.
  (parallel) IO and Posix IO. The "diag_observer" observer is a basic 
  diagnostic capability that will be generalized in the futurre to more
  cleanly accommodate state vectors with many different component types, 
  and it will likely also become configurable.

  Following are the descriptions of the "binary_observer" parameters:

  | expression           | description                                     |
  |----------------------|-------------------------------------------------|
  | "observer_name"      | name of observer in src/pdeint/observer_factory.cpp|
  | "indirectory"        | input directory name (for reads)|
  | "outdirectory"       | output directory name (for writes)|
  | "cadence_type"       | either "cycle" or time, as in Sec. F.|
  | "cycle_interval"     | cycle cadence if outputting in terms of "cycle"'s|
  | "time_interval"      | time cadence if outputting in terms of evol. time|
  | "state_index"        | array of state component indices to output|
  | "state_names"        | file name prefixes for component indices in "state_index" |
  | "grid_names"         | names for the grid files; there should be >= GDIMS entries|
  | "agg_state_name"     | ``aggregate'' state file name if all state components are placed in a single file. Whether this can be done is determined by the capabilities of the specified "IO_implementation"|
  | "agg_grid_name"      | ``aggregate'' grid file name if all grid components are placed in a single file. Whether this can be done is determined by the capabilities of the specified "IO_implementation".|
  | "derived_quantities" | quantities to compute from basic state for output only. These are not evolved. These names refer to self-referential config blocks within the JSON file|

  | "filename_size"      | maximum filename size|

  Under the hood, the "binary_observer" can take any type of IO implementation, in principle. 
  The implementation used is specified in the "IO_implementation" variable in the JSON 
  file, which 
  names an implementation type in src/pdeint/io_factory.ipp. The implementation determines, for
  instance, how IO is done, whether collectively or task-by-task. Currently, there is only one such
  implementation available, the "gio" implementation, and it provides IO according
  to the GIO class in src/cdg/io/gio.hpp. In the future, it may be desirable to have, say, 
  a NetCDF IO implementation, which will have its own set of configuration parameters.
  
  The following are the allowed parameters for configuring the "gio" IO implementation:

  ```json
  "gio": {
    "ivers"    : 0,
    "multivar" : false,
    "io_type"  : "collective",
    "wtime"    : 6,
    "wtask"    : 5,
    "wfile"    : 2048
  },
  ```
  with the configuration parameters described here:

  | expression           | description                                     |
  |----------------------|-------------------------------------------------|
  | "ivers"    | version number. Default is 0, and that's the only one available currently|
  | "multivar" | flag telling GIO to place all state variables in the same file|
  | "io_type"  | may be "collective" or "POSIX" for collective (MPIIO) or by-task IO, respetively|
  | "wtime"    | width of time index field in filename|
  | "wtask"    | width of MPI task field when doing "POSIX" IO|
  | "wfile"    | maximum width of filename|
  
H. Configure external forcing.

  External forcing configuration works almost exactly like the state initialization described 
  in Sec. E. But because the procedure by which the forcing may be updated is not
  yet generalized, we postpone further discussion of this topic. 

I. Configure terrain.

  GeoFLOW provides the ability to generate grid terrain. This is done by starting with
  a basic grid (described in Sec, C.), then providing the grid terrain as a boundary 
  condition to the 
  Laplace equation for each grid terrain coordinate. This yields a Poisson equation for
  each grid coordinate. Once we solve for the new grid coordinates, we then modify the
  original grid with the new grid solutions. 

  Terrain configuration consists of setting the parameter "terrain_type" to
  a config block name. The config block then contains the name of a valid
  terrain specification method in src/cdg/gspecterrain_factory.ipp, and the
  associated parameters.

  For example, assume GDIM=2, and we have provided a "grid_box" grid as in the
  JSON file sample. Then

  ```json
  "terrain_type"         : "myboxterrain",
  "myboxterrain": {
    "name"    : "boxgauss_range",
    "x0"      : [0.25, 0.75],
    "y0"      : [0.0, 0.0],
    "xsigma"  : [0.1, 0.2]
    "ysigma"  : [0.1, 0.2],
    "h0"      : [0.15, 0.4]
  },
  ```
  specifies that we will use the "boxgauss_range" method to generate terrain. This method
  allows one to configure a superposition of as many Gaussian ``hills'' as desired, with 
  different locations, 
  FWHM values and heights, in order to model a mountain range. When GDIM=2, "y0" should be
  0.0, and the y-FWHM will not be used. The same procedure is used to specify terrain 
  for 2d and 3d spherical grids.

  Valid terrain "name" parameters for spherical and box grids are set in 
  src/cdg/init/gspecterrain_factory.ipp.

  Note that the three parameters in the grid config block 
     ```json
     "maxit"     : 128,
     "tol"       : 1.0e-8,
     "norm_type" : "GCG_NORM_INF"
     ```
  control the iterative Krylov solution for the Poisson grid solution, as described
  in the "grid_box" table in Sec. C..


