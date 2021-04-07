//==================================================================================
// Module       : gsponge_bdy.hpp
// Date         : 7/5/20 (DLR)
// Description  : Object that handles the the time updates for
//                'sponge' boundaries.
//
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : UpdateBdyBase.
//==================================================================================




//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method 
// DESC   : 
// RETURNS: none
//**********************************************************************************
template<typename Types>
GSpongeBdy<Types>::GSpongeBdy(typename GSpongeBdy<Types>::Traits &traits) :
UpdateBdyBase<Types>(),
traits_                 (traits)
{
  // Do some checks:
  assert(traits_.istate  .size() == traits_.farfield.size());
  assert(traits_.exponent.size() == 1 
      || traits_.exponent.size() == traits_.farfield.size());
  assert(traits_.sigma   .size() == 1 
      || traits_.sigma   .size() == traits_.farfield.size());
  assert( abs(traits_.idir) > 0 && traits_.idir <= GDIM );

} // end of constructor method 


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename Types>
GSpongeBdy<Types>::~GSpongeBdy()
{
} // end, destructor


//**********************************************************************************
//**********************************************************************************
// METHOD : update_impl
// DESC   : Entry method for doing a sponge-layer update
// ARGS   : 
//          eqn   : equation implementation
//          grid  : grid object (necessary?)
//          time  : timestep
//          utmp  : tmp vectors
//          u     : state vector
// RETURNS: none.
//**********************************************************************************
template<typename Types>
GBOOL GSpongeBdy<Types>::update_impl(
                              EqnBasePtr &eqn,
                              Grid       &grid,
                              Time       &time,
                              State      &utmp,
                              State      &u)
{
   GString    serr = "GSpongeBdy<Types>::update_impl: ";
   GBOOL      bret = FALSE;
   GGridBox<Types>   *box    = dynamic_cast<GGridBox<Types>*>(&grid);
   GGridIcos<Types>  *sphere = dynamic_cast<GGridIcos<Types>*>(&grid);


   if ( box != NULLPTR ) {
     bret = update_cart(eqn, grid, time, utmp, u);
   }
   else if ( sphere != NULLPTR ) {
     bret = update_sphere(eqn, grid, time, utmp, u);
   }
   else {
     assert(FALSE && "Invalid grid");
   }

   return bret;

} // end of method update_impl


//**********************************************************************************
//**********************************************************************************
// METHOD : update_cart
// DESC   : Method for doing a sponge-layer update on Cartesian grids
// ARGS   : 
//          eqn   : equation implementation
//          grid  : grid object (necessary?)
//          time  : timestep
//          utmp  : tmp vectors
//          u     : state vector
// RETURNS: none.
//**********************************************************************************
template<typename Types>
GBOOL GSpongeBdy<Types>::update_cart(
                              EqnBasePtr &eqn,
                              Grid       &grid,
                              Time       &time,
                              State      &utmp,
                              State      &u)
{
  GString          serr = "GSpongeBdy<Types>::update_cart: ";
  Time             tt = time;
  GINT             idstate;
  GSIZET           ind;
  GFTYPE           beta, ifact, rtst, sig0, sgn;
  GTVector<GSIZET> *igbdy = &traits_.ibdyvol;

  GTVector<GTVector<GFTYPE>> 
                  *xnodes = &grid.xNodes();

  assert(GDIM == 3);

  // Get parameters from ptree:

  ifact    = 1.0/(traits_.ro - traits_.rs[0]);

  // This method applies a sponge layer to only the outer
  // part of a artesian grid. The first values in
  // traits.rs, is used to set the coord value in the 
  // direction traits.idir that defines the sponge layer.
  // The value traits.idir is signed , so that 
  //   traits.idir X ( r - rs ) > 0 defines the r values
  // that sit in the layer.

  sgn = traits_.idir / abs(traits_.idir);

  assert(FALSE);

  // Update state due to sponge layer:
  // Note: This is equiavalent to adding a dissipation 
  //       term to the RH of the operator-split equation, s.t.:
  //        du/dt = -sig(r) (u - u_infinity)
  //       where
  //        sig(r) = sig_0 [(r - rs)/(ro - rs)]^exponent
  //       and u_infinity is the far-field solution
  // Note: We may have to re-form this scheme if we use semi-implicit
  //       or implicit time stepping methods!
  for ( auto k=0; k<traits_.istate.size(); k++ ) { // for each state component
    idstate = traits_.istate[k];
    for ( auto j=0; j<u[idstate]->size(); j++ ) { // for all grid points
      rtst = sgn * ( (*xnodes)[abs(traits_.idir-1)][k] - traits_.rs[0] );
      beta = rtst > 0 ? pow(ifact*fabs(rtst),traits_.exponent[k]) : 0.0; // check if in sponge layer
//    (*u[idstate])[j]  = (1.0-beta)*(*u[idstate])[j] + beta*traits_.farfield[k];
      sig0 = traits_.sigma.size() > 1 ? traits_.sigma[k] : traits_.sigma[0];
      (*u[idstate])[j] -= sig0*beta*( (*u[idstate])[j] - traits_.farfield[k] );
    }
  }


  return TRUE;

} // end of method update_cart


//**********************************************************************************
//**********************************************************************************
// METHOD : update_sphere
// DESC   : Method for doing a sponge-layer update on spherical 
//          (3d only) grids
// ARGS   : 
//          eqn   : equation implementation
//          grid  : grid object (necessary?)
//          time  : timestep
//          utmp  : tmp vectors
//          u     : state vector
// RETURNS: none.
//**********************************************************************************
template<typename Types>
GBOOL GSpongeBdy<Types>::update_sphere (
                              EqnBasePtr &eqn,
                              Grid       &grid,
                              Time       &time,
                              State      &utmp,
                              State      &u)
{
  GString          serr = "GSpongeBdy<Types>::update_sphere: ";
  Time             tt = time;
  GINT             idstate;
  GSIZET           ind;
  GFTYPE           beta, ifact, sig0;
  GFTYPE           r, x, y, z;
  GTVector<GSIZET> *igbdy = &traits_.ibdyvol;

  GTVector<GTVector<GFTYPE>> 
                  *xnodes = &grid.xNodes();



  assert(GDIM == 3);

  // Get parameters from ptree:

  ifact    = 1.0/(traits_.ro - traits_.rs[0]);

  // This method applies a sponge layer to only the outer
  // part of a spherical grid. Thus, only first values in
  // traits.rs, and traits.ro are used to define inner and
  // outer radii (which is grid radius)

  // Update state due to sponge layer:
  // Note: This is equiavalent to adding a dissipation 
  //       term to the RH of the operator-split equation, s.t.:
  //        du/dt = -sig(r) (u - u_infinity)
  //       where
  //        sig(r) = sig_0 [(r - rs)/(ro - rs)]^traits_.exponent
  //       and u_infinity is the far-field solution
  // Note: We may have to re-form this scheme if we use semi-implicit
  //       or implicit time stepping methods!
  // Note: traits.idir is ignored here, since, for the sphere,
  //       sponge layers are only defined in the radial 
  //       direction in idirection of outer boundary
  for ( auto k=0; k<traits_.istate.size(); k++ ) { // for each state component
    idstate = traits_.istate[k];
    for ( auto j=0; j<u[idstate]->size(); j++ ) { // for all grid points
      x    = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
      r    = sqrt(x*x + y*y + z*z); 
      beta = r >= traits_.rs[0] ? pow(ifact*(r-traits_.rs[0]),traits_.exponent[k]) : 0.0; // check if in sponge layer
//    (*u[idstate])[j]  = (1.0-beta)*(*u[idstate])[j] + beta*traits_.farfield[k];
      if ( traits_.sigma.size() > 1 )
        (*u[idstate])[j] -= traits_.sigma[k]*beta*( (*u[idstate])[j] - traits_.farfield[k] );
      else
        (*u[idstate])[j] -= traits_.sigma[0]*beta*( (*u[idstate])[j] - traits_.farfield[k] );
    }
  }
 
  // Set bdy vectors:
  for ( auto k=0; k<traits_.istate.size(); k++ ) {
    idstate = traits_.istate[k];

    // Set from initialized State vector,
    for ( auto j=0; j<igbdy->size(); j++ ) {
      ind = (*igbdy)[j];
      (*u[idstate])[ind] = traits_.farfield[k];
    }
  }


  return TRUE;
} // end of method update_sphere


