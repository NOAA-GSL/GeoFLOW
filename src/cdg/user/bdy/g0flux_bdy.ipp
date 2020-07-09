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
template<typename TypePack>
G0FluxBdy<TypePack>::G0FluxBdy(G0FluxBdy<TypePack>::Traits &traits) :
UpdateBdyBase<TypePack>(),
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
template<typename TypePack>
G0FluxBdy<TypePack>::~G0FluxBdy()
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
template<typename TypePack>
GBOOL G0FluxBdy<TypePack>::update_impl(
                              Grid      &grid,
                              StateInfo &stinfo,
                              Time      &time,
                              State     &utmp,
                              State     &u,
                              State     &ub)
{
   GString    serr = "G0FluxBdy<TypePack>::update_impl: ";
   GBdyType   itype;
   GINT       id, ind, k;
   GSIZET     il;
   Ftype      sum, xn;
   GTVector<GTVector<Ftype>>  *bdyNormals;
   GTVector<GINT>             *idepComp;
   GTVector<GTVector<GSIZET>> *igbdy;
   GTVector<GTVector<GSIZET>> *ilbdy;


  GTVector<GTVector<GSIZET>> *igbdy = &grid.igbdy_binned();

  if ( traits_.compute_once && bcomputed_ ) return TRUE;

  for ( auto j=0; j<traits_.istate.size(); j++ ) {
    assert(ub[istate[j]] != NULL && "Illegal bdy vector");
  }

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
  igbdy      = &grid_->igbdy_binned();
  ilbdy      = &grid_->ilbdy_binned();
  itype      = GBDY_0FLUX;
  for ( auto j=0; j<(*igbdy)[itype].size() ) {
    ind = (*igbdy)[itype][j]; // index into long vector array
    il  = (*ilbdy)[itype][j]; // index into bdy array (for normals, e.g.)
    id  = (*idep)[ind];       // dependent vector component
    xn  = (*n)[id][ind];      // n_id == normal component for dependent vector comp
    sum = 0.0;
    for ( auto k=0; k<nstate_; k++ ) { // for each vector component
      sum -= (*n)[k][il] * (*u[istate[k]])[ind];
//    if ( k != (*idep)[ind] ) sum -= (*n)[idstate][il] * (*u[idstate])[ind];

      // Set all comps here, then reset the dependent one:
      (*ub[k])[ind] = (*u[k])[ind]; 
    }
    (*ub[id])[ind] = ( sum + (*n)[id][il] * (*u[id])[ind] ) / xn;
  }

  bcomputed = TRUE;

  return TRUE;

} // end of method update_impl


