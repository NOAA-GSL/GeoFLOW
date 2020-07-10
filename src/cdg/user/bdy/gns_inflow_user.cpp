//==================================================================================
// Module       : ginitbdy_user.cpp
// Date         : 7/11/19 (DLR)
// Description  : Boundary initialization function implementations provided
//                by user
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#include "ginitbdy_user.hpp"


namespace ginitbdy {




//**********************************************************************************
//**********************************************************************************
// METHOD : impl_mybdyinit
// DESC   : Initialize bdy using state initialization method 
//          specified by ptree
// ARGS   : ptree : main prop tree
//          grid  : grid
//          stinfo: StateInfo
//          t     : time
//          utmp  : tmp arrays
//          u     : current state, overwritten here
//          ub    : bdy vectors (one for each state element)
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_mybdyinit(const PropertyTree &ptree, GGrid &grid, StateInfo &stinfo, Time &time, State &utmp, State &u, State &ub)
{
#if 0
  Time             tt = t;
  GString          serr = "impl_mybdyinit: ";


  GTVector<GTVector<GSIZET>> *igbdy = &grid.igbdy();


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

  return TRUE;

#else

  /* 
     Fill in here; change function name 
     if desired. Add (unique) function name to 
     ginitb_factory.ipp.
  */


  return FALSE;
#endif

} // end, impl_mybdyinit




} // end, ginitbdy namespace
