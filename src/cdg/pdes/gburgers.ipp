//==================================================================================
// Module       : gburgers.ipp
// Date         : 10/18/18 (DLR)
// Description  : Object defining a multidimensional Burgers (advection-diffusion) 
//                PDE. 
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include "gburgers.hpp"



//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Instantiate with 2 fields, so a 2d problem
// ARGS   : grid  : grid object
//          u, v, : fields
//          tmp   : tmp vectors of same size as u, v
// RETURNS: none
//**********************************************************************************
GBurgers::GBurgers(GGrid &ggrid, GTVector<GFTYPE> &u, GTVector<GFTYPE> &v, GTVector<GFTYPE> &tmp) :
ggrid_      (&ggrid),
u_          (&u),
v_          (&v),
tmp_        (&tmp),
bdycallback_(NULLPTR)
{
  init2d();
} // end of constructor method (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (2)
// DESC   : Instantiate with 3 fields, so a 3d problem. NOTE: do we allow 2.5d?
// ARGS   : grid    : grid object
//          u, v, w : fields
//          tmp     : tmp vectors of same size as u, v
// RETURNS: none
//**********************************************************************************
GBurgers::GBurgers(GGrid &ggrid, GTVector<GFTYPE> &u, GTVector<GFTYPE> &v, GTVector<GFTYPE> &w, GTVector<GFTYPE> &tmp) :
ggrid_      (&ggrid),
u_          (&u),
v_          (&v),
v_          (&w),
tmp_        (&tmp),
bdycallback_(NULLPTR)
{
  init3d();
} // end of constructor method (2)



//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
GBurgers::~GBurgers()
{
} // end, destructor


//**********************************************************************************
//**********************************************************************************
// METHOD : dt_impl
// DESC   : Compute time step
// ARGS   : t : time
//          u : state
//          dt: timestep, returned
// RETURNS: none.
//**********************************************************************************
void GBurgers::dt_impl(const Time &t, State &u, Time &dt)
{
    dt = 1.0;

} // end of method dt_impl


//**********************************************************************************
//**********************************************************************************
// METHOD : dudt_impl
// DESC   : Compute RHS for explicit schemes
// ARGS   : t  : time
//          u  : state
//          dt : time step
// RETURNS: none.
//**********************************************************************************
void GBurgers::dudt_impl(const Time &t, State &u, Time &dt, Derivative &dudt)
{


} // end of method dudt_impl


//**********************************************************************************
//**********************************************************************************
// METHOD : dfdu_impl
// DESC   : Compute Jacobian dF/du for F == RHS
// ARGS   : t  : time
//          u  : state
//          dt : time step
// RETURNS: none.
//**********************************************************************************
void GBurgers::dfdu_impl(const Time &t, State &u, Time &dt, Jacobian &dfdu)
{

} // end of method dfdu_impl


//**********************************************************************************
//**********************************************************************************
// METHOD : set_bdy_callback
// DESC   : Set the callback object and method pointers for setting bdy conditions
// ARGS   : ptr2obj : pointer to object that hosts method callback
//          callback: method name
// RETURNS: none.
//**********************************************************************************
void GBurgers::set_bdy_callback(std::function<void(GGrid &)> &callback)
{
  bdycallback_  = &callback;
} // end of method set_bdy_callback


//**********************************************************************************
//**********************************************************************************
// METHOD : init2d
// DESC   : Initialize base state/base icosahedron
// ARGS   : none.
// RETURNS: none.
//**********************************************************************************
void GBurgers::init2d()
{
  GString serr = "GBurgers::init2d: ";


} // end of method init2d


//**********************************************************************************
//**********************************************************************************
// METHOD : init3d
// DESC   : Initialize for 3d elements
// ARGS   : none.
// RETURNS: none.
//**********************************************************************************
void GBurgers::init3d()
{

} // end, method init3d


