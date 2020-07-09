//==================================================================================
// Module       : g0flux_bdy.hpp
// Date         : 7/5/20 (DLR)
// Description  : Object that handles the the time updates for
//                0-flux boundaries
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
bcomput_once_            (FALSE),
traits_                 (traits)
{

  assert( traits_.istate.size() == GDIM 
       && "Kinetic vector must be specified");

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
GBOOL G0FluxBdy<TypePack>::update_impl(
                              Grid      &grid,
                              StateInfo &stinfo,
                              Time      &time,
                              State     &utmp,
                              State     &u,
                              State     &ub)
{
   GString    serr = "G0FluxBdy<TypePack>::update_impl: ";
   GBOOL      bret;
   GINT       idstate, ind;

  GTVector<GTVector<GSIZET>> *igbdy = &grid.igbdy_binned();

  if ( bcompute_once_ && bcomputed_ ) return TRUE;


  // Set boundary vector with initialized state:
  for ( auto n=0; n<traits_.istate.size() && bret; n++ ) { 
    idstate = traits_.istate[n];
    if ( stinfo.compdesc[idstate] == GSC_PRESCRIBED
      || stinfo.compdesc[idstate] == GSC_NONE ) continue;
    }
    // Set from initialized State vector, 
    for ( auto j=0; j<(*igbdy)[GBDY_0FLUX].size()
       && ub[idstate] != NULLPTR; j++ ) {
      ind = (*igbdy)[GBDY_0FLUX][j];
      (*ub[idstate])[j] = 0
    }
  }
  bcomputed = bret;

  return bret;

} // end of method update_impl


