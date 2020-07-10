//==================================================================================
// Module       : goutflow_bdy.hpp
// Date         : 7/5/20 (DLR)
// Description  : Object that handles the the time updates for
//                outflow boundaries
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
GOutflowBdyBdy<TypePack>::GOutflowBdyBdy(GOutflowBdyBdy<TypePack>::Traits &traits) :
UpdateBdyBase<TypePack>(),
traits_                 (traits)
bcomputed_               (FALSE)
{

  assert(FALSE && "Not yet implemented");

} // end of constructor method 


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename TypePack>
GOutflowBdyBdy<TypePack>::~GOutflowBdyBdy()
{
} // end, destructor


//**********************************************************************************
//**********************************************************************************
// METHOD : update_impl
// DESC   : Entry method for updateing a simple outflow bdy condition
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
GBOOL GOutflowBdyBdy<TypePack>::update_impl(
                              Grid      &grid,
                              StateInfo &stinfo,
                              Time      &time,
                              State     &utmp,
                              State     &u,
                              State     &ub)
{
   GString    serr = "GOutflowBdyBdy<TypePack>::update_impl: ";
   GINT       idstate, ind;

  GTVector<GTVector<GSIZET>> *igbdy = &grid.igbdy_binned();

  if ( traits_.compute_once && bcomputed_ ) return TRUE;

  assert( u.size() == stinfo.compdesc.size() && "State info structure invalid");

  // Set from State vector, u:
  for ( auto n=0; n<traits_.istate.size(); n++ ) { // apply to specified state comps
    idstate = traits_.istate[n];
    if ( stinfo.compdesc[idstate] == GSC_PRESCRIBED
      || stinfo.compdesc[idstate] == GSC_NONE ) continue;
    for ( auto j=0; j<(*igbdy)[GBDY_OUTFLOW].size()
       && ub[idstate] != NULLPTR; j++ ) {
      ind = (*igbdy)[GBDY_OUTFLOW][j];
      (*ub[idstate])[j] = (*u[idstate])[(*igbdy)[ind];
    }
  }
  bcomputed_ = TRUE;

  return TRUE;

} // end of method update_impl


