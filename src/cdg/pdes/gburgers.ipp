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
// METHOD : Constructor method 
// DESC   : Instantiate with grid + state + tmp 
// ARGS   : grid  : grid object
//          u     : state 
//          tmp   : tmp vectors of same size as u, v
// RETURNS: none
//**********************************************************************************
GBurgers::GBurgers(GGrid &grid, State &u, GTVector<GTVectorGFTYPE>*> &tmp) :
ggrid_      (&ggrid)
{
  static_assert(std::is_same<State,GTVector<GTVectorGFTYPE>>>::value,
               "State is of incorrect type"); 

  init();
  
} // end of constructor method 



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
// DESC   : Compute Jacobian dF/du for F == RHS. Used to linearize
//          equations for purely implicit schemes.
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


