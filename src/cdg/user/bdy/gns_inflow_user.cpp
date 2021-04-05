//==================================================================================
// Module       : ginitbdy_user.cpp
// Date         : 7/11/19 (DLR)
// Description  : Boundary initialization function implementations provided
//                by user
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#include "gns_inflow_user.hpp"


namespace GInflowBdyMethods {




//**********************************************************************************
//**********************************************************************************
// METHOD : myinflow
// DESC   : Place holder template for inflow update methods
// ARGS   : 
//          grid  : grid
//          t     : time
//          id    : canonical bdy id
//          utmp  : tmp arrays
//          u     : current state, overwritten here
//          ub    : bdy vectors (one for each state element)
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL myinflow(GGrid &grid, Time &time, const GINT id, State &utmp, State &u, State &ub)
{
#if 0
  Time             tt = t;
  GString          serr = "impl_mybdyinit: ";


  GTVector<GTVector<GSIZET>> *igbdy = &grid.igbdy_binned()[id];


  // Set from State vector, u and others that we _can_ set:
  for ( auto k=0; k<u.size(); k++ ) { 
    for ( auto j=0; j<(*igbdy)[GBDY_DIRICHLET].size() 
       && ub[k] != NULLPTR; j++ ) {
      (*ub[k])[j] = (*u[k])[(*igbdy)[GBDY_DIRICHLET][j]];
    }
    for ( auto j=0; j<(*igbdy)[GBDY_NOSLIP].size() 
       && ub[k] != NULLPTR; j++ ) {
      (*ub[k])[j] = 0.0;
    }
  }

#endif
  
   return TRUE;

} // end, myinflow




} // end, GInflowBdyMethods
