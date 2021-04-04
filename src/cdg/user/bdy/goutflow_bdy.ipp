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
template<typename Types>
GOutflowBdy<Types>::GOutflowBdy(typename GOutflowBdy<Types>::Traits &traits) :
UpdateBdyBase<Types>(),
traits_                 (traits)
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
template<typename Types>
GOutflowBdy<Types>::~GOutflowBdy()
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
template<typename Types>
GBOOL GOutflowBdy<Types>::update_impl(
                              Grid      &grid,
                              StateInfo &stinfo,
                              Time      &time,
                              State     &utmp,
                              State     &u,
                              State     &ub)
{
   GString    serr = "GOutflowBdy<Types>::update_impl: ";
   GINT       idstate;
   GSIZET     ind;

  GTVector<GSIZET> *igbdy = &traits_.ibdyvol;


  assert( u.size() == stinfo.icomptype.size() && "State info structure invalid");

  // Set from State vector, u:
  for ( auto n=0; n<traits_.istate.size(); n++ ) { // apply to specified state comps
    idstate = traits_.istate[n];
    if ( stinfo.icomptype[idstate] == GSC_PRESCRIBED
      || stinfo.icomptype[idstate] == GSC_NONE ) continue;
    for ( auto j=0; j<igbdy->size()
       && ub[idstate] != NULLPTR; j++ ) {
      ind = (*igbdy)[j];
      (*u[idstate])[ind] = (*u[idstate])[ind];
    }
  }

  return TRUE;

} // end of method update_impl


