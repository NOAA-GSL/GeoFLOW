//==================================================================================
// Module       : ginitfv.cpp
// Date         : 7/16/19 (DLR)
// Description  : Velocity force inititial condition implementations 
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#include "ginitfv.hpp"

namespace ginitfv {


//**********************************************************************************
//**********************************************************************************
// METHOD : impl_abc_box
// DESC   : Inititialize velocity with Arnold-Beltrami-Childress (ABC)
//          initial conditions for box grids, 2d and 3d.
// ARGS   : ptree  : main property tree
//          sconfig: ptree block name containing variable config
//          grid   : grid object
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          uf     : velocity-state to be initialized.
// RETURNS: TRUE on success; else FALSE
//**********************************************************************************
GBOOL impl_abc_box(const PropertyTree &ptree, GString &sconfig, GGrid &grid, Time &time, State &utmp, State &ub, State &uf)
{

  assert(grid.gtype() == GE_REGULAR && "Box grids required");


  GINT         kdn, kup, p, pdef; 
  GSIZET       nn ;
  GFTYPE       A, B, C, f0, pi2, x, y, z;
  PropertyTree vtree;
  GTVector<GTVector<GFTYPE>>
              *xnodes = &grid.xNodes();

#if defined(_G_IS3D)
  pdef = 2;
#else
  pdef = 3;
#endif

  vtree = ptree.getPropertyTree(sconfig);
  kdn   = vtree.getValue<GINT>("kdn");
  kup   = vtree.getValue<GINT>("kup");
  p     = vtree.getValue<GINT>("kpower",pdef);
  A     = vtree.getValue<GFTYPE>("A", 0.9);
  B     = vtree.getValue<GFTYPE>("B", 1.0);
#if defined(_G_IS3D)
  C     = vtree.getValue<GFTYPE>("C", 1.1);
#endif
  f0    = vtree.getValue<GFTYPE>("f0", 1.0);
  nn    = (*xnodes)[0].size();

  // Stream fcn is
  //   psi = Sum_i { -A cos(2pi*ki*x) + B sin(2pi*ki*y) / ki^p }
  // Compute vel components s.t. ux = d psi / dy, uy = -d psi / dx
#if defined(_G_IS2D)

  *uf[0] = 0.0;
  *uf[1] = 0.0;
  for ( GSIZET j=0; j<nn; j++ ) {
    x = (*xnodes)[0][j]; y = (*xnodes)[1][j];
    for ( GINT k=kdn; k<=kup; k++ ) {
      pi2         = 2.0*PI*k;
      (*uf[0])[j] +=  B*pi2*cos(pi2*y) / pow(k,p);
      (*uf[1])[j] += -A*pi2*sin(pi2*x) / pow(k,p);
    }
  }
 
#elif defined(_G_IS3D)

  *uf[0] = 0.0;
  *uf[1] = 0.0;
  *uf[2] = 0.0;
  for ( GSIZET j=0; j<nn; j++ ) {
    x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j]; 
    for ( GINT k=kdn; k<kup; k++ ) {
      pi2         = 2.0*PI*k;
      (*uf[0])[j] +=  ( B*cos(pi2*y) + C*sin(pi2*z) ) / pow(k,p);
      (*uf[1])[j] +=  ( A*sin(pi2*x) + C*cos(pi2*z) ) / pow(k,p);
      (*uf[2])[j] +=  ( A*cos(pi2*x) + B*sin(pi2*y) ) / pow(k,p);
    }
  }

#endif

  GMTK::normalize_euclidean(uf, NULLPTR, 0, f0);

  return TRUE;

} // end, method impl_abc_box




//**********************************************************************************
//**********************************************************************************
// METHOD : impl_rand
// DESC   : Inititialize velocity with Gaussian-randomized values
// ARGS   : ptree  : initial condition property tree
//          sconfig: ptree block name containing variable config
//          grid   : grid object
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          uf     : velocity-state to be initialized.
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_rand(const PropertyTree &ptree, GString &sconfig, GGrid &grid, Time &time, State &utmp, State &ub, State &uf)
{

  return FALSE;

} // end, method impl_rand


} // end, ginitfv  namespace
