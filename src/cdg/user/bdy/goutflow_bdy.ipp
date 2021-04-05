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
//          eqn   : equation implenetation
//          grid  : grid object (necessary?)
//          time  : timestep
//          utmp  : tmp vectors
//          u     : state vector
// RETURNS: none.
//**********************************************************************************
template<typename Types>
GBOOL GOutflowBdy<Types>::update_impl(
                              EqnBasetPtr &eqn,
                              Grid        &grid,
                              Time        &time,
                              State       &utmp,
                              State       &u)
{
   GString    serr = "GOutflowBdy<Types>::update_impl: ";
   GINT       idstate;
   GSIZET     ind;

  GTVector<GSIZET> *igbdy = &traits_.ibdyvol;



  // Set from State vector, u:
  for ( auto n=0; n<traits_.istate.size(); n++ ) { // apply to specified state comps
    idstate = traits_.istate[n];
    for ( auto j=0; j<igbdy->size()
       && ub[idstate] != NULLPTR; j++ ) {
      ind = (*igbdy)[j];
      (*u[idstate])[ind] = (*u[idstate])[ind];
    }
  }

  return TRUE;

} // end of method update_impl


