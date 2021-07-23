//==================================================================================
// Module       : gsponge_bdy.hpp
// Date         : 7/5/20 (DLR)
// Description  : Object that handles the the time updates for
//                'sponge' boundaries (Rayleigh filtering).
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
binit_                   (FALSE),
xmax_                      (0.0),
traits_                 (traits)
{
  // Do some checks:
  assert(traits_.farfield.size() == traits_.istate.size());
  assert(traits_.falloff .size() == traits_.istate.size());
  assert(traits_.exponent.size() == traits_.istate.size());
  assert(abs(traits_.idir) > 0 && abs(traits_.idir) <= GDIM );


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
   GridBox   *box    = dynamic_cast<GridBox*>(&grid);
   GridIcos  *sphere = dynamic_cast<GridIcos*>(&grid);

   if ( !binit_ ) init(grid);

   if ( box != NULLPTR ) {
     bret = update_box(eqn, grid, time, utmp, u);
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
// METHOD : update_box
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
GBOOL GSpongeBdy<Types>::update_box(
                              EqnBasePtr &eqn,
                              Grid       &grid,
                              Time       &time,
                              State      &utmp,
                              State      &u)
{
  GString          serr = "GSpongeBdy<Types>::update_box: ";
  Time             dt, tt = time;
  GINT             adir, idstate;
  GSIZET           j, ind;
  Ftype            eps = 1e2*std::numeric_limits<Ftype>::epsilon();
  Ftype            ifact, rate, rtst, sgn;
  vector<GSIZET> *igbdy = &traits_.ibdyvol;

  GTVector<GTVector<Ftype>> 
                  *xnodes = &grid.xNodes();


  ifact    = 1.0/fabs(xmax_ - traits_.xstart);
  dt       = tt - told_;

  // This method applies a sponge layer to only the outer
  // part of a Cartesian grid. 
  // The value traits.idir is signed , so that 
  //   traits.idir X ( r - rs ) > 0 defines the r values
  // that reside in the layer.

  sgn = traits_.idir / abs(traits_.idir);
  adir = abs(traits_.idir);
  assert(traits_.falloff.size() >= traits_.istate.size());
  assert(traits_.exponent.size() >= traits_.istate.size());

  // Update state due to sponge layer:
  // Note: This is equiavalent to adding a dissipation 
  //       term to the RH of the operator-split equation, s.t.:
  //        du/dt = -sig(r) (u - u_infinity)
  //       where
  //        sig(r) = sig_0 [(r - rs)/(rmax - rs)]^exponent
  //       and u_infinity is the far-field solution
  // Note: We have to re-examine this scheme if we use semi-implicit
  //       or implicit time stepping methods!
Ftype del;
  for ( auto k=0; k<traits_.istate.size(); k++ ) { // for each state component
    idstate = traits_.istate[k];
    for ( auto jj=0; jj<isponge_.size(); jj++ ) { // for all grid points
//  for ( auto j=0; j<(*xnodes)[0].size(); j++ ) { // for all grid points
      j = isponge_[jj];
//    rtst = sgn * ( (*xnodes)[adir-1][j] - traits_.xstart );
      rtst =       ( (*xnodes)[adir-1][j] - traits_.xstart );
//    rate = rtst > 0 ? pow(ifact*fabs(rtst),traits_.exponent[k]) : 0.0; // check if in sponge layer
      rate = pow(ifact*fabs(rtst),traits_.exponent[k])* traits_.falloff[k]; 
//    (*u[idstate])[j]  = (1.0-rate)*(*u[idstate])[j] + rate*traits_.farfield[k];
// del = dt*traits_.falloff[k]*rate*( (*u[idstate])[j] - traits_.farfield[k] );
//if ( idstate == 0 ) 
//cout << "del = " << del << endl;
      (*u[idstate])[j] -= rate*( (*u[idstate])[j] - traits_.farfield[k] );
    }
  }

  told_ = tt;

  return TRUE;

} // end of method update_box


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
  Time             dt, tt = time;
  GINT             idstate;
  GSIZET           j, ind;
  Ftype           ifact, rate;
  Ftype           r, x, y, z;
  vector<GSIZET> *igbdy = &traits_.ibdyvol;

  GTVector<GTVector<Ftype>> 
                  *xnodes = &grid.xNodes();


  assert(GDIM == 3);

  // Get parameters from ptree:

  ifact    = 1.0/(xmax_ - traits_.xstart);
  dt       = tt = told_;
  assert(traits_.falloff.size() >= traits_.istate.size());
  assert(traits_.exponent.size() >= traits_.istate.size());

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
    for ( auto jj=0; jj<isponge_.size(); jj++ ) { // for all grid points
      j    = isponge_[jj];
      x    = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
      r    = sqrt(x*x + y*y + z*z); 
      rate = r >= traits_.xstart ? pow(ifact*(r-traits_.xstart),traits_.exponent[k]) : 0.0; // check if in sponge layer
//    (*u[idstate])[j]  = (1.0-rate)*(*u[idstate])[j] + rate*traits_.farfield[k];
      (*u[idstate])[j] -= dt*traits_.falloff[k]*rate*( (*u[idstate])[j] - traits_.farfield[k] );
    }
  }
  
 
#if 0
  // Set bdy vectors:
  for ( auto k=0; k<traits_.istate.size(); k++ ) {
    idstate = traits_.istate[k];

    // Set from initialized State vector,
    for ( auto j=0; j<igbdy->size(); j++ ) {
      ind = (*igbdy)[j];
      (*u[idstate])[ind] = traits_.farfield[k];
    }
  }
#endif
  told_ = tt;

  return TRUE;

} // end of method update_sphere


//**********************************************************************************
//**********************************************************************************
// METHOD : init
// DESC   : initialization method
// ARGS   : grid: Grid object
// RETURNS: none.
//**********************************************************************************
template<typename Types>
void GSpongeBdy<Types>::init(Grid &grid)
{
   GString    serr = "GSpongeBdy<Types>::init: ";
   GINT       adir, idir;
   GSIZET     j, n, nxy;
   Ftype      eps = 1.0e2*numeric_limits<Ftype>::epsilon();
   Ftype      del, r, sgn;
   GridBox   *box    = dynamic_cast<GridBox*>(&grid);
   GridIcos  *sphere = dynamic_cast<GridIcos*>(&grid);
   GTVector<GSIZET>
              itmp;
   GTVector<GTVector<Ftype>>  *xnodes = &grid.xNodes();
   GTMatrix<GTVector<Ftype>>  *dXdXi  = &grid.dXdXi();
   typename Grid::GElemList   *elems  = &grid.elems();

   assert( ( box != NULLPTR || sphere != NULLPTR )&& "Invalid grid");

   itmp.resize((*xnodes)[0].size());

   idir = traits_.idir; // is signed
   adir = abs(idir);
   sgn  = idir / adir;
 

   assert(grid.ispconst()); // Object requires const p


   // First, find indices of all grid points that
   // reside in the sponge layer:
   n = 0;
   if ( box != NULLPTR ) { // is a box grid
     assert ( adir >= 1 && adir <= GDIM );
     for ( auto j=0; j<(*xnodes)[0].size(); j++ ) {
       r = (*xnodes)[adir-1][j] ;
//     del = (*xnodes)[adir-1][j] - traits_.xstart;
       if ( (sgn > 0 && r >= traits_.xstart) 
         || (sgn < 0 && r <= traits_.xstart) )
          itmp[n++] = j;
     }
     xmax_ = idir > 0 ? (*xnodes)[adir-1].max() 
                      : (*xnodes)[adir-1].min();
     isponge_.resize(n);
     isponge_.set(itmp.data(), n);
   }
   else if ( GDIM == 3 && sphere != NULLPTR ) { // is a spherical grid
     assert ( adir == 1 ); // radial direction only
     for ( auto j=0; j<(*xnodes)[0].size(); j++ ) {
       r = 0.0;
       for ( auto i=0; i<GDIM; i++ ) r += pow((*xnodes)[i][j],2);
       r = sqrt(r);
       del = r - traits_.xstart;
       if ( (sgn > 0 && del >= 0.0) || (sgn < 0 && del <= 0.0) )
          itmp[n++] = j;
     }
     for ( auto i=0; i<GDIM; i++ ) r += pow((*xnodes)[i].amax(),2);
     xmax_ = sqrt(r);
 
     isponge_.resize(n);
     isponge_.set(itmp.data(), n);
   }

   binit_ = TRUE;

} // end of method init

