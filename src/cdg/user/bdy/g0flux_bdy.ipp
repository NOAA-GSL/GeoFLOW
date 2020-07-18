//==================================================================================
// Module       : g0flux_bdy.hpp
// Date         : 7/5/20 (DLR)
// Description  : Object that handles the the time updates for
//                0-flux boundaries. Acts on kinetic
//                vector.
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
G0FluxBdy<Types>::G0FluxBdy(typename G0FluxBdy<Types>::Traits &traits) :
UpdateBdyBase<Types>(),
bcomputed_               (FALSE),
nstate_                      (0),
traits_                 (traits)
{

  assert( traits_.istate.size() == GDIM 
       && "Kinetic vector must be specified");

  nstate_ = traits_.istate.size();

} // end of constructor method 


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename Types>
G0FluxBdy<Types>::~G0FluxBdy()
{
} // end, destructor


//**********************************************************************************
//**********************************************************************************
// METHOD : update_impl
// DESC   : Entry method for updating 0-Flux bdy conditions
// ARGS   : 
//          grid  : grid object (necessary?)
//          stinfo: state info structure
//          time  : timestep
//          utmp  : tmp vectors
//          u     : state vector
//          ub    : bdy vector
// RETURNS: none.
//**********************************************************************************
template<typename Types>
GBOOL G0FluxBdy<Types>::update_impl(
                              Grid      &grid,
                              StateInfo &stinfo,
                              Time      &time,
                              State     &utmp,
                              State     &u,
                              State     &ub)
{
   GString    serr = "G0FluxBdy<Types>::update_impl: ";
   GBdyType   itype;
   GINT       nid, k;
   GSIZET     il, ind;
   Ftype      sum, xn;
   GTVector<GTVector<Ftype>>  *bdyNormals;
   GTVector<GINT>             *idepComp;
   GTVector<GSIZET>           *igbdy;
// GTVector<GTVector<GSIZET>> *ilbdy;



  if ( traits_.compute_once && bcomputed_ ) return TRUE;


#if 0
  for ( auto j=0; j<traits_.istate.size(); j++ ) {
    assert(ub[istate[j]] != NULL && "Illegal bdy vector");
  }
#endif

  // Handle 0-Flux bdy conditions. This
  // is computed by solving
  //    vec{n} \cdot vec{u} = 0
  // for 'dependent' component set in grid.
  //
  // Note: We may want to switch the order of the
  //       following loops to have a better chance
  //       of vectorization. Unrolling likely
  //       won't occur:
  bdyNormals = &grid.bdyNormals(); // bdy normal vector
  idepComp   = &grid.idepComp();   // dependent components
  igbdy      = &traits_.ibdy;
  ilbdy      = &grid_->ilbdy_binned()[traits_.bdyid];
  itype      = GBDY_0FLUX;
  for ( auto j=0; j<igbdy->size() ) {
    ind = (*igbdy)[i][j];         // index into volume array
//  il  = (*ilbdy)[i][itype][j];  // index into bdy array (for normals, e.g.)
    nid = (*idep)[ind];           // dependent vector component
    xn  = (*bdyNormals)[nid][ind];// n_id == normal component for dependent vector comp
    sum = 0.0;
    for ( auto k=0; k<nstate_; k++ ) { // for each vector component
      sum -= (*bdyNormals)[k][il] * (*u[istate[k]])[ind];

      // Set all comps here, then reset the dependent one:
      (*ub[k])[ind] = (*u[k])[ind]; 
    }
//  (*ub[nid])[ind] = ( sum + (*n)[nid][il] * (*u[nid])[ind] ) / xn;
    (*u[nid])[ind] = ( sum + (*n)[nid][il] * (*ua[nid])[ind] ) / xn;
  }

  bcomputed_ = TRUE;

  return TRUE;

} // end of method update_impl


