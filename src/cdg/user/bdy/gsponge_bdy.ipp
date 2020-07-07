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
template<typename TypePack>
GSpongeBdy<TypePack>::GSpongeBdy(GSpongeBdy<TypePack>::Traits &traits) :
UpdateBdyBase<TypePack>(),
traits_                 (traits)
{
  // Do some checks:
  assert(traits_.istate  .size() == traits_.farfield.size());
  assert(traits_.exponent.size() == 1 
      || traits_.exponent.size() == traits_.farfield.size());
  assert(traits_.sigma   .size() == 1 
      || traits_.signa   .size() == traits_.farfield.size());
  assert( abs(traits_.idir) > 0 && traits_.idir <= GDIM );

} // end of constructor method 


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename TypePack>
GSpongeBdy<TypePack>::~GSpongeBdy()
{
} // end, destructor


//**********************************************************************************
//**********************************************************************************
// METHOD : update_impl
// DESC   : Entry method for doing a sponge-layer update
// ARGS   : 
//          grid  : grid object (necessary?)
//          stinfo: state info structure
//          time  : timestep
//          utmp  : tmp vectors
//          u     : state vector
//          ub    : bdy vector
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
GBOOL GSpongeBdy<TypePack>::update_impl(
                              Grid      &grid,
                              StateInfo &stinfo,
                              Time      &time,
                              State     &utmp,
                              State     &u,
                              State     &ub)
{
   GString    serr = "GSpongeBdy<TypePack>::update_impl: ";
   GBOOL      bret = FALSE;
   GGridBox   *box = std::dynamic_cast<GGridBox*>(grid_);
   GGridIcos  *sphere = std::dynamic_cast<GGridIcos*>(grid_);

   if ( box != NULLPTR ) {
     bret = update_cart(grid, stinfo, time, utmp, u, ub);
   }
   else if ( sphere != NULLPTR ) {
     bret = update_sphere(grid, stinfo, time, utmp, u, ub);
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
//          grid  : grid object (necessary?)
//          stinfo: state info structure
//          time  : timestep
//          utmp  : tmp vectors
//          u     : state vector
//          ub    : bdy vector
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
GBOOL GSpongeBdy<TypePack>::update_cart(
                              Grid      &grid,
                              StateInfo &stinfo,
                              Time      &time,
                              State     &utmp,
                              State     &u,
                              State     &ub)
{
  GString          serr = "GSpongeBdy<TypePack>::update_cart: ";
  Time             tt = time;
  GINT             idstate;
  GFTYPE           beta, ifact, rtst, sig0;
//GTVector<GTVector<GINT>> 
//                *igbdycf = &grid.igbdycf_binned(); 
//GTVector<GTVector<GSIZET>> 
//                *igbdy = &grid.igbdy_binned();

  GTVector<GTVector<GFTYPE>> 
                  *xnodes = &grid.xNodes();

  assert(GDIM == 3);

  // Get parameters from ptree:

  ifact    = 1.0/(traits_.ro[0] - traits_.rs[0]);

  // This method applies a sponge layer to only the outer
  // part of a artesian grid. The first values in
  // traits.rs, is used to set the coord value in the 
  // direction traits.idir that defines the sponge layer.
  // The value traits.idir is signed , so that 
  //   traits.idir X ( r - rs ) > 0 defines the r values
  // that sit in the layer.

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
    if ( stinfo.compdesc[idstate] == GSC_PRESCRIBED
      || stinfo.compdesc[idstate] == GSC_NONE ) continue;
    for ( auto j=0; j<u[idstate]->size(); j++ ) { // for all grid points
//    igb  = (*igbdy)  [GBDY_SPONGE][j];
//    igf  = (*igbdycf)[GBDY_SPONGE][j];
      rtst = traits_.idir*( (*xnodes)[abs(idir-1)][k] - rs[0] );
      beta = rtst > 0 ? pow(ifact*fabs(rtst),exponent) : 0.0; // check if in sponge layer
//    (*u[idstate])[j]  = (1.0-beta)*(*u[idstate])[j] + beta*farfield[k];
      sig0 = traits_.sigma.size() > 1 ? traits_.sigma[k] : traits_.sigma[0];
      (*u[idstate])[j] -= sig0*beta*( (*u[idstate])[j] - farfield[k] );
    }
  }

   return bret;

} // end of method update_cart


//**********************************************************************************
//**********************************************************************************
// METHOD : update_sphere
// DESC   : Method for doing a sponge-layer update on spherical 
//          (3d only) grids
// ARGS   : 
//          grid  : grid object (necessary?)
//          stinfo: state info structure
//          time  : timestep
//          utmp  : tmp vectors
//          u     : state vector
//          ub    : bdy vector
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
GBOOL GSpongeBdy<TypePack>::update_sphere (
                              Grid      &grid,
                              StateInfo &stinfo,
                              Time      &time,
                              State     &utmp,
                              State     &u,
                              State     &ub)
{
  GString          serr = "GSpongeBdy<TypePack>::update_sphere: ";
  Time             tt = time;
  GINT             idstate;
  GFTYPE           beta, ifact, sig0;
  GFTYPE           r, x, y, z;
//GTVector<GTVector<GINT>> 
//                *igbdycf = &grid.igbdycf_binned(); 
//GTVector<GTVector<GSIZET>> 
//                *igbdy = &grid.igbdy_binned();

  GTVector<GTVector<GFTYPE>> 
                  *xnodes = &grid.xNodes();

  assert(GDIM == 3);

  // Get parameters from ptree:

  ifact    = 1.0/(traits_.ro[0] - traits_.rs[0]);

  // This method applies a sponge layer to only the outer
  // part of a spherical grid. Thus, only first values in
  // traits.rs, and traits.ro are used to define inner and
  // outer radii (which is grid radius)

  // Update state due to sponge layer:
  // Note: This is equiavalent to adding a dissipation 
  //       term to the RH of the operator-split equation, s.t.:
  //        du/dt = -sig(r) (u - u_infinity)
  //       where
  //        sig(r) = sig_0 [(r - rs)/(ro - rs)]^exponent
  //       and u_infinity is the far-field solution
  // Note: We may have to re-form this scheme if we use semi-implicit
  //       or implicit time stepping methods!
  // Note: traits.idir is ignored here, since, for the sphere,
  //       sponge layers are only defined in the radial 
  //       direction in idirection of outer boundary
  for ( auto k=0; k<traitstraits__istate.size(); k++ ) { // for each state component
    idstate = traits_.istate[k];
    for ( auto j=0; j<u[idstate]->size(); j++ ) { // for all grid points
      x    = (*xnodes)[0][k]; y = (*xnodes)[1][k]; z = (*xnodes)[2][k];
      r    = sqrt(x*x + y*y + z*z); 
//    igb  = (*igbdy)  [GBDY_SPONGE][j];
//    igf  = (*igbdycf)[GBDY_SPONGE][j];
      beta = r >= rs[0] ? pow(ifact*(r-rs[0]),exponent) : 0.0; // check if in sponge layer
//    (*u[idstate])[j]  = (1.0-beta)*(*u[idstate])[j] + beta*farfield[k];
      if ( traits_.sigma.size() > 1 )
        (*u[idstate])[j] -= traits_.sigma[k]*beta*( (*u[idstate])[j] - farfield[k] );
      else
        (*u[idstate])[j] -= traits_.sigma[0]*beta*( (*u[idstate])[j] - farfield[k] );
    }
  }


} // end of method update_sphere


