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
   GINT       idd, k;
   GSIZET     il, iloc, ind;
   Ftype      sum, xn;
   GTVector<GTVector<Ftype>>  *bdyNormals;
// GTVector<GTVector<Ftype>>  *xnodes;
   GTVector<GINT>             *idep;
   GTVector<GSIZET>           *igbdy;
// GTVector<GTVector<GSIZET>> *ilbdy;



  if ( traits_.compute_once && bcomputed_ ) return TRUE;


  // Handle 0-Flux bdy conditions. This
  // is computed by solving
  //    vec{n} \cdot vec{u} = 0
  // for 'dependent' component set in grid.
  //
  // Note: We may want to switch the order of the
  //       following loops to have a better chance
  //       of vectorization. Unrolling likely
  //       won't occur:
//xnodes     = &grid.xNodes();
  bdyNormals = &grid.bdyNormals(); // bdy normal vector
  idep       = &grid.idepComp();   // dependent components
  igbdy      = &traits_.ibdyvol;
//ilbdy      = &grid_->ilbdy_binned()[traits_.bdyid];
  itype      = GBDY_0FLUX;
  for ( auto j=0; j<igbdy->size(); j++ ) {
    iloc = traits_.ibdyloc[j];     // index into bdy arrays
    ind  = (*igbdy)[j];            // index into volume array
    idd  = (*idep)[iloc];          // dependent vector component

    if ( idd >= 0 ) {
      xn   = (*bdyNormals)[idd][iloc];// n_idd == normal component for dependent vector comp
      sum  = 0.0;
      for ( auto k=0; k<nstate_; k++ ) { // for each indep vector component
        sum += idd == traits_.istate[k] ? 0.0 
                    : (*bdyNormals)[k][iloc] * (*u[traits_.istate[k]])[ind];
      }
    
      // Ensure v.n = 0:
//    (*u[idd])[ind] = ( sum + (*bdyNormals)[idd][iloc] * (*u[idd])[ind] ) / xn;
      (*u[idd])[ind] = -sum / xn;
    }
    else {
        // Set v == 0, say, on a vertex:
        for ( auto k=0; k<nstate_; k++ ) (*u[traits_.istate[k]])[ind] = 0.0;
    }
  }

  bcomputed_ = TRUE;


  return TRUE;

} // end of method update_impl


