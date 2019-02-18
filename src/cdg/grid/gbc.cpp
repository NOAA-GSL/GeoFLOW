//==================================================================================
// Module       : gbc
// Date         : 2/6/19 (DLR)
// Description  : GeoFLOW CDG boundary condition update object
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#include "gbc.hpp"


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Default constructor
// ARGS   : grid: GGrid operator
// RETURNS: none
//**********************************************************************************
GBC::GBC(GGrid &grid)
:
grid_                    (&grid), 
update_dirichlet_callback_ (NULLPTR)
{
} // end of constructor method (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
GBC::~GBC()
{
} // end, destructor


//**********************************************************************************
//**********************************************************************************
// METHOD : update_bdy_state
// DESC   : Update (Dirichlet) bdy vectors, ub
// ARGS   : t    : time
//          u    : current state
//          ub   : bdy vectors (one for each state element)
// RETURNS: none.
//**********************************************************************************
void GBC::update_bdy_state(GFTYPE &t, GTVector<GTVector<GFTYPE>*> &u,
                           GTVector<GTVector<GFTYPE>*> &ub)
{
  GTVector<GTVector<GSIZET>> *igbdy = &grid_->igbdy();

  assert((*igbdy)[GBDY_NEUMANN].size() == 0
      && (*igbdy)[GBDY_OUTFLOW].size() == 0
      && (*igbdy) [GBDY_SPONGE].size() == 0
      && "Invalid GBdyType specified" );

  // Set Dirichlet conditions:

  // ...GBDY_NOSLIP:
  for ( GSIZET k=0; k<u.size(); k++ ) { 
    for ( GSIZET j=0; j<(*igbdy)[GBDY_NOSLIP].size(); j++ ) {
      (*ub[k])[(*igbdy)[GBDY_NOSLIP][j]] = 0.0;
    }
  }

  // ...GBDY_0FLUX:
  for ( GSIZET k=0; k<u.size(); k++ ) { 
    for ( GSIZET j=0; j<(*igbdy)[GBDY_0FLUX].size(); j++ ) {
      (*ub[k])[(*igbdy)[GBDY_0FLUX][j]] = 0.0;
    }
  }

  // ...GBDY_REFLECT:
  if ( (*igbdy)[GBDY_REFLECT].size() > 0 ) {
    do_reflective_bcs(*grid_, t, u, ub);
  }

  // ...GBDY_DIRICHLET:
  // (These are possibly time dependent, and are often
  //  set with data outside the scope of the solvers)
  if ( update_dirichlet_callback_ != NULLPTR) {
    (*update_dirichlet_callback_)(t, u, ub);
  }


} // end of update_bdy_state


//**********************************************************************************
//**********************************************************************************
// METHOD : do_reflective_bcs
// DESC   : Update reflective bdy vectors, ub
// ARGS   : grid : GGrid object
//          t    : time
//          u    : current state
//          ub   : bdy vectors (one for each state element)
// RETURNS: none.
//**********************************************************************************
void GBC::do_reflective_bcs(GGrid &grid, GFTYPE t, GTVector<GTVector<GFTYPE>*> &u,
                            GTVector<GTVector<GFTYPE>*> &ub)
{

  assert(FALSE && "Reflective conditions not available");

} // end of do_reflective_bcs


