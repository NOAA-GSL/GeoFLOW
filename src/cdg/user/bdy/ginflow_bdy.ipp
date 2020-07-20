//==================================================================================
// Module       : ginflow_frominit_bdy.ipp
// Date         : 7/5/20 (DLR)
// Description  : Object that handles the the time updates for
//                inflow boundaries that are set by initialization method
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
GInflowBdy<Types>::GInflowBdy(typename GInflowBdy<Types>::Traits &traits) :
UpdateBdyBase<Types>(),
bcomputed_               (FALSE),
traits_                 (traits)
{


} // end of constructor method 


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename Types>
GInflowBdy<Types>::~GInflowBdy()
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
template<typename Types>
GBOOL GInflowBdy<Types>::update_impl(
                              Grid      &grid,
                              StateInfo &stinfo,
                              Time      &time,
                              State     &utmp,
                              State     &u,
                              State     &ub)
{
   GString    serr = "GInflowBdy<Types>::update_impl: ";
   GBOOL      bret;

   if ( traits_.use_init ) {
     bret = update_from_init(grid, stinfo, time, utmp, u, ub);
   }
   else {
     bret = update_from_user(grid, stinfo, time, utmp, u, ub);
   }

   return bret;

} // end of method update_impl


//**********************************************************************************
//**********************************************************************************
// METHOD : update_from_user
// DESC   : Entry method for updating inflow bdys from 
//          user-specified function
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
GBOOL GInflowBdy<Types>::update_from_user(
                              Grid      &grid,
                              StateInfo &stinfo,
                              Time      &time,
                              State     &utmp,
                              State     &u,
                              State     &ub)
{
  GString    serr = "GInflowBdy<Types>::update_from_user: ";
  GBOOL      bret;

  assert(traits_.callback != NULLPTR);  

  bret = traits_.callback(grid, stinfo, time, traits_.bdyid,  utmp, u, ub);
  bcomputed_ = bret;

  return bret;

} // end of method update_from_user


//**********************************************************************************
//**********************************************************************************
// METHOD : update_from_init
// DESC   : Entry method for updating inflow bdys from 
//          initialization state. Useful for testing.
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
GBOOL GInflowBdy<Types>::update_from_init(
                              Grid      &grid,
                              StateInfo &stinfo,
                              Time      &time,
                              State     &utmp,
                              State     &u,
                              State     &ub)
{
   GString    serr = "GInflowBdy<Types>::update_from_init: ";
   GBOOL      bret;
   GINT       idstate, ind;
   State      tmp;

  GTVector<GSIZET> *igbdy = &traits_.ibdyvol;

  if ( traits_.compute_once && bcomputed_ ) return TRUE;

  assert( u.size() == stinfo.icomptype.size() && "State info structure invalid");

  assert(utmp.size() >= u.size());
  if ( unew_.size() < u.size() ) {
    unew_.resizem(u.size());
    tmpnew_.resizem(utmp.size()-u.size());
  }
  for ( auto j=0; j<u.size(); j++ ) unew_ [j] = utmp[j];
  for ( auto j=0; j<utmp.size()-u.size(); j++ ) tmpnew_[j] = utmp[u.size()+j];

  // Call initialization method with utmp:
  bret = GInitStateFactory<Types>::init(traits_.ptree, grid, stinfo, time, tmpnew_, ub, unew_);

  // Set boundary vector with initialized state:
  for ( auto n=0; n<traits_.istate.size() && bret; n++ ) { 
    idstate = traits_.istate[n];
    if ( stinfo.icomptype[idstate] == GSC_PRESCRIBED
      || stinfo.icomptype[idstate] == GSC_NONE ) continue;
    
    // Set from initialized State vector, 
    for ( auto j=0; j<igbdy->size()
       && ub[idstate] != NULLPTR; j++ ) {
      ind = (*igbdy)[j];
//    (*ub[idstate])[j] = (*unew_[idstate])[ind];
      (*u [idstate])[ind] = (*unew_[idstate])[ind];
    }
  }
  bcomputed_ = bret;

  return bret;

} // end of method update_from_init


