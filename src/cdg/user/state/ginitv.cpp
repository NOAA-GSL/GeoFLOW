//==================================================================================
// Module       : ginitv.cpp
// Date         : 7/16/19 (DLR)
// Description  : Velocity inititial condition implementations 
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#include "ginitv.hpp"


namespace ginitv {


//**********************************************************************************
//**********************************************************************************
// METHOD : impl_abc_box
// DESC   : Inititialize velocity with Arnoldi-Beltrami-Childress (ABC)
//          initial conditions for box grids, 2d and 3d.
// ARGS   : ptree  : main property tree
//          sconfig: ptree block name containing variable config
//          grid   : grid object
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          u      : velocity-state to be initialized.
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_abd_box(const PropertyTree &ptree, GString &sconfig, GGrid &grid, Time &time, State &utmp, State &ub, State &u)
{

  assert(grid.gtype() == GE_REGULAR && "Box grids required");

  PropertyTree vtree = ptree.getPropertyTree(sconfig);

  GINT    kdn    = vtree.getValue  <GINT>("kdn); 
  GINT    kup    = vtree.getValue  <GINT>("kup"); 
  GINT    p      = vtree.getValue  <GINT>("kpower",3);
  GSIZET  nn ;
  GFTYPE  A      = vtree.getArray<GFTYPE>("A", 1.1); 
  GFTYPE  B      = vtree.getArray<GFTYPE>("B"  2.3); 
#if defined(_G_IS3D)
  GFTYPE  C      = vtree.getArray<GFTYPE>("C", 2.6); 
#endif
  GFTYPE  u0     = vtree.getArray<GFTYPE>("u0", 1.0); 
  GFTYPE  pi2, x, y, z; 
  GTVector<GTVector<GFTYPE>> 
         *xnodes = &grid.xNodes();

  nn = xnodes[0]->size();

  // Stream fcn is 
  //   psi = Sum_i { -A cos(2pi*ki*x) + B sin(2pi*ki*y) / ki^p }
  // Compute vel components s.t. ux = d psi / dy, uy = -d psi / dx
#if defined(_G_IS2D)

  *u[0] = 0.0;
  *u[1] = 0.0;
  for ( GSIZET j=0; j<nn; j++ ) {
    x = (*xnodes[0])[j]; y = (*xnodes[1])[j]; 
    for ( GINT k=kdn; k<=kup; k++ ) {
      pi2         = 2.0*PI*k;
      (*u[0])[j] +=  B*pi2*cos(pi2*y) / pow(k,p);
      (*u[1])[j] += -A*pi2*sin(pi2*x) / pow(k,p);
    }
  }
  
#elif defined(_G_IS3D)

  *u[0] = 0.0;
  *u[1] = 0.0;
  *u[2] = 0.0;
  for ( GSIZET j=0; j<nn; j++ ) {
    x = (*xnodes[0])[j]; y = (*xnodes[1])[j]; 
    for ( GINT k=kdn; k<kup; k++ ) {
      pi2         = 2.0*PI*k;
      (*u[0])[j] +=  ( B*cos(pi2*y) + C*sin(pi2*z) ) / pow(k,p);
      (*u[1])[j] +=  ( A*sin(pi2*x) + C*cos(pi2*z) ) / pow(k,p);
      (*u[2])[j] +=  ( A*cos(pi2*x) + B*sin(pi2*y) ) / pow(k,p);
    }
  }

#endif

  GMTK::normalize_euclidean(u, NULLPTR, 0, u0);

  return TRUE;

} // end, method impl_rand

//**********************************************************************************
//**********************************************************************************
// METHOD : impl_rand
// DESC   : Inititialize velocity with Gaussian-randomized values
// ARGS   : ptree  : main property tree
//          sconfig: ptree block name containing variable config
//          grid   : grid object
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          u      : state to be initialized.
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_rand(const PropertyTree &ptree, GString &sconfig, GGrid &grid, Time &time, State &utmp, State &ub, State &u)
{

  return FALSE;

} // end, method impl_rand


} // end, ginitv  namespace
