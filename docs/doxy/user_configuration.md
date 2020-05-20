  GeoFLOW: How to set up a problem and interact with the system

#I. Introduction.

    The entry point to run all problems is to use a user-configuralbe JSON file.
The default file name is 'input.jsn' it is assumed to lie in the the run 
(working) directory. On the command line, this default config file name may be changed 
by using the "-i" command line argument:

> ./geoflow_cg -i my_new_config_file.jsn


  Inside this config text file, there are several items that must be set in order
to run a problem. User must at a minumum:
(A) specify a solver name that tells us which initial value-boundary value 
    (IVBV) partial differential equation set to solve; 
(B) specify a time-stepping scheme
(C) specify which grid to use, and probably also the basis expansion order
(D) specify a set of initial conditions
(E) specify a set of boundary conditions
(F) specify the evolution (integration) time period
(G) configure output: specify which variables to output, and their cadences.

Occasionally, the user may also wish to: 

(H) specify external forcing
(I) specify grid terrain.

Each of these specifications is described here in turn below. One additional
quantity that is not user-configurble at run-time is the problem dimension, GDIM,
which is specified by defining the preprocessor variables _G_IS2D, or _G_IS3D 
at build time. This is described in XXX. We simply refer to this as necessary
below.

#II. Config file example.

Before getting started, here is a sample configuration file. Let's assume it's 
given the name "input.jsn" if we need that later. We will refer to this as 
we proceed, and the user will have a concrete example to serve as a template.
For the time being, we assume that we are onfiguring GeoFLOW Spectral Element 
discretizations only.

{
  "pde_name"             : "pde_burgers",
  "exp_order"            : [4, 4, 4],
  "exp_order_type"       : "constant",
  "grid_type"            : "grid_icos",
  "initstate_block"      : "initstate_icosnwave",
  "initforce_block"      : "",
  "use_forcing"          : false,
  "initstate_type"       : "direct",
  "initforce_type"       : "",
  "restart_index"        : 0,
  "benchmark"            : false,
  "do_comparison"        : true,
  "observer_list"        : ["gio_observer", "global_diag_basic"],
  "terrain_type"         : ""
  "initv": {
    "name"  : "abc",
    "kdn"   : 2,
    "kup"   : 3,
    "kpower": 0.5,
    "A"     : 0.9,
    "B"     : 1.0,
    "C"     : 1.1,
    "E0"    : 1.0
  },
  "mystateinit": {
    "initv"   : "initvss",
    "initps"  : "zero",
    "initc"   : "zero"
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
  "gio_observer": {
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
    "time_field_width"   : 6,
    "task_field_width"   : 5,
    "filename_size"      : 2048
  },
  "global_diag_basic": {
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


#III. Describe each configuration step.

    ##(A) Specify a solver name.

        To specify a PDE solver, simply provide it in the input.jsn file by 
setting:

    "pde_name" : "pde_burgers"

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
    this will be discussed later.

    ##(B) Specify a solver name.

        In the current configuration, the time stepping configuration is not "attached"
    to the PDE solver; it a future version it likely will be, and it will take an
    arbitrary name. For now, this block is called "stepper_props", and contains the 
    following configurable items
    	
    "stepper_props": {
      "stepping_method"  : "GSTEPPER_EXRK",
      "time_deriv_order" : 4,
      "extrap_order"     : 2,
      "variable_dt"      : false,
      "courant"          : 0.04
    },

    where the items are described here:

    "stepping_method": May be one of the valid entries in src/cdg/include/gtypes.h:
    in the sGStepperType array. At present, only explicit Runge-Kutta integration is enabled.

    "time_deriv_order": refers to the time truncation order of the scheme (for all 
    available orders. The current scheme allows up to order 5. 

    "extrap_order": refers to the extrapolation order, if this method were to be enabled.
    
    "variable_dt": tells the stepper that the timestep may vary from step to step. If true,
    a method is called that limits the size of the step based on the current state.

    "courant": is the Courant number: a "fudge" factor that multiplies the timestep
    from the timestep method when "variable_dt" = true, so that the Courant condition
    isn't violated

    ##(C) Select and configure grid.

        The grid type is specified at the top of the JSON file by setting the 
    "grid_type" variable. Like with PDE specification, this is a quantity that
    will also name a self-referential configation block within the file. 

    "grid_type"   : "grid_icos",

    Valid "grid_type"'s are: "grid_icos", or "grid_box".  

    If GDIM=2, then "grid_box" will refer to a planar Cartesian grid; if GDIM=3, 
    it will refer to a 3D Cartesian box grid. 

    If GDIM=2, the "grid_icos" refers to a 2D spherical surface grid, constructed from
    a base icosohedron; if GDIM=3, "grid_icos" uses the 2D "grid_icos" grid at the 
    surface, and constructs elements radially from this "base" 2D grid.

      ###(i) Configuring the box grids:

      The self-referential config block for the "grid_box" grids  look like:

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

      where entries and their descriptions are:

      "grid_name"         : "grid_box", not required
      "xyz0"              : tuple of startig lower left of box. Must be of size >= GDIM.
      "delxyz"            : tuple of widths in each direction. Must be of size >= GDIM.
      "num_elems"         : number of elements in each direction. Must be of size >= GDIM.
      "bdy_x_0"           : name of config block defining boundary conditions on west face
      "bdy_x_1"           : name of config block defining boundary conditions on east face
      "bdy_y_0"           : name of config block defining boundary conditions on south face
      "bdy_y_1"           : name of config block defining boundary conditions on north face
      "bdy_z_0"           : name of config block defining boundary conditions on bottom face. 
                            Used if GDIM=3
      "bdy_z_1"           : name of config block defining boundary conditions on top face
                            Used if GDIM=3
      "bdy_init_method"   : if specified, determines how boundaries are initialized
      "bdy_update_method" : if specified, determines how boundaries are updated
      "use_state_init"    : If specified, the boundary initialization for DIRICHLET type 
                            conditions will be taken from the initial conditions
      "maxit"             : max no. Krylov iterations if using terrain
      "tol"               : Krylov loop tolerance is using terrain
      "norm_type"         : norm to use in establishing Krylov loop residual. Valid values
                            are provided in src/pdeint.in_solver_base.hpp: "GCG_NORM_INF",
 			    "GCG_NORM_EUC", "GCG_NORM_L2", "GCG_NORM_L1". Used if doing terrain.
      

       #### Boundary condition specification block. 

           In the "grid_box" configuration, each enumerated natural boundary is configured 
       to use a "fully_peridoc" boundary condition in each direction via a statement like:
      
       "bdy_x_0"           : "fully_periodic"
       "bdy_x_1"           : "fully_periodic"

       The name on the right specifies a self-referenced configuration block in the 
       JSON file that actually configures that boundary. It's important to remember 
       that these blocks govern _specification_ of the boundary, not, e.g., boundary 
       _initialization_. Boundary initialization and evolution is considered below. 

       These boundary spec blocks, in general, have the form:

       "fully_periodic": {
         "base_type"        : "GBDY_PERIODIC"
         "bdy_class"        : "uniform"
         "bdy_config_method": ""
       },

       Valid "base_type" values are provided in 

         src/cdg/include/gtypes.h

       in the the sGBdyType array. Currently, BDY_DIRICHLET (or GBDY_INFLOWT, 
       these are the same thing), and GBDY_PERIODIC work currently. The types GBDY_0FLUX 
       and GBDY_SPONGE will be implemented soon. 

       Valid "bdy_class" expressions may be: "uniform" or "user". 
     
       "uniform" : the condition is applied automatically to the entire natural boundary
       "user"    : if specified the user MUST specify the "bdy_config_method".
       
       The named "bdy_config_method" is used to set the boundary condition type. 
       Methods that may be provided are specified in the factory associated with 
       boundary specification: 
         src/cdg/init/gspecbdy_factory.cpp. 
       The steps by which to add a new boundary specification method is discussed later.

       "GBDY_PERIODIC" boundary conditions are somewhat special. They may be specified
       for only "grid_box" boundaries, and, if applied, must be applied to the entire 
       boundary; i.e., it will require that the "bdy_class" is "uniform". Also, they must 
       be applied pairwise between associated boundaries. So the following would be 
       incorrect:

       "bdy_x_0"   : "fully_periodic",
       "bdy_x_1"   : "burgers_inflow".

       Note that if the user wants "bdy_class" to be "uniform", since this is the default,
       neither the "bdy_class" nor the "bdy_config_method" need to be provided in the
       boundary spec block; these two parameters need be specified only if "bdy_class" 
       is set to "user". 

       Note that a "user" specified spec method may require other parameters
       to be specified in the spec block.
       

       #### Boundary condition initialization and updating. 

           If a boundary type is GBDY_DIRICHLET/GBDY_INFLOWT or GBDY_SPONGE (when 
       ready), the boundary values must be initialized and evolved. Here we describe
       how this is done.
      
       Two variables in the "grid_box" config block determine how the boundaries are
       initialized:

       "use_state_init"    : true,
       "bdy_init_method  " : "",
       "bdy_update_method" : "",

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

      ###(ii) Configuring spherical grids:

      When GDIM=2, we use the following configuration block for a spherical
      surface:
      "grid_icos": {
        "grid_name"         : "grid_icos",
        "radius"            : 1.0,
        "ilevel"            : 2,
        "refine_type"       : "GICOS_BISECTION",
      },

      where

      "grid_name"         : "grid_box", not required
      "radius"            : scalar radius
      "refine_type"       : type of refinement: either "GICOS_BISECTION"" or "GICOS_LAGRANGIAN" 
      "ilevel"            : if "refine_type" = "GICOS_BISECTION", then this is the number
                            of bisection levels, yielding a grid with 60*(2^ilevel) elements.
                            If "refine_type" = "GICOS_LAGRANGIAN", ilevel gives the number of
                            divisions along a triangle edge, yielding 60*(ilevel+1)^2 
                            elements on the grid.

      Note: there are no boundaries on the 2d sphere, so no boundary conditions need to
      be specified.


      When GDIM=3, we use the following configuration block for a spherical
      surface:

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

      where entries and their descriptions are:
       
      "grid_name"         : "grid_box", not required
      "radiusi"           : scalar inner radius
      "radiuso"           : scalar outer radius (not used if GDIM=2)
      "refine_type"       : type of refinement: either "GICOS_BISECTION"" or "GICOS_LAGRANGIAN" 
      "ilevel"            : if "refine_type" = "GICOS_BISECTION", then this is the number
                            of bisection levels, yielding a grid with 60*(2^ilevel) elements.
                            If "refine_type" = "GICOS_LAGRANGIAN", ilevel gives the number of
                            divisions along a triangle edge, yielding 60*(ilevel+1)^2 
                            elements on the grid.
      "num_radial_elems"  : number of elements in radial direction
      "bdy_inner"         : name of config block defining boundary conditions on inner surface
      "bdy_outer"         : name of config block defining boundary conditions on outersurface
      "bdy_init_method"   : if specified, determines how boundaries are initialized
      "bdy_update_method" : if specified, determines how boundaries are updated
      "use_state_init"    : If specified, the boundary initialization for DIRICHLET type 
                            conditions will be taken from the initial conditions
      
       Note: The 3d spherical grid is built upon the 2d surface grid by extending elements
       radially. The boundary condition spec, initilization, and update methods work
       the same was as for box grids.


    ##(D) Select and configure grid.



