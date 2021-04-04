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
  // Allocate bdy arrays foreach state component:
  bdydata_.resize(traits_.istate.size());
  bdydata_ = NULLPTR;
  for ( auto j=0; j<traits_.istate.size(); j++ ) {
    bdydata_[j]  = new GTVector<Ftype>(traits_.ibdyvol.size());
  }

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

  for ( auto j=0; j<bdydata_.size(); j++ ) {
    if ( bdydata_[j] != NULLPTR ) delete bdydata_[j];
  }

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
   GINT       idstate;
   GSIZET     ind;

   // Compute bdy data, if necessary:
   bret = compute_bdy_data(grid, stinfo, time, utmp, u, bdydata_);
   assert(bret);
   
   // Apply bdy data:
   for ( auto n=0; n<traits_.istate.size() && bret; n++ ) { 
     idstate = traits_.istate[n];
    
     // Set bdy values computed above:
     for ( auto j=0; j<traits_.ibdyvol.size(); j++ ) {
       ind = traits_.ibdyvol[j];
       (*u[idstate])[ind] = (*bdydata_[n])[j];
     }
   }

   return bret;

} // end of method update_impl


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_bdy_data
// DESC   : Compute boundary data from state.
// ARGS   : 
//          grid  : grid object (necessary?)
//          stinfo: state info structure
//          time  : timestep
//          utmp  : tmp vectors
//          u     : state vector
//          ub    : bdy vector(s), returned
// RETURNS: none.
//**********************************************************************************
template<typename Types>
GBOOL GInflowBdy<Types>::compute_bdy_data(
                              Grid      &grid,
                              StateInfo &stinfo,
                              Time      &time,
                              State     &utmp,
                              State     &u,
                              State     &ub)
{
   GString    serr = "GInflowBdy<Types>::compute_bdy_data: ";
   GBOOL      bret;
   GINT       idstate;
   GSIZET     ind;
   State      tmpnew, unew;

  GTVector<GSIZET> *igbdy = &traits_.ibdyvol;

  if ( traits_.compute_once && bcomputed_ ) return TRUE;


  assert(utmp.size() >= u.size());
  if ( unew.size() < u.size() ) {
    unew.resize(u.size());
    tmpnew.resize(utmp.size()-u.size());
  }
  for ( auto j=0; j<u.size(); j++ ) unew [j] = utmp[j];
  for ( auto j=0; j<tmpnew.size(); j++ ) tmpnew[j] = utmp[unew.size()+j];

  // Call initialization method with utmp:
   if ( traits_.use_init ) {
     bret = GInitStateFactory<Types>::init(traits_.ptree, grid, stinfo, time, tmpnew, ub, unew);
   }
   else {
     bret = traits_.callback(grid, stinfo, time, traits_.bdyid, tmpnew, unew,  ub);
   }
   assert(bret);

  // Set boundary vector with initialized state:
  for ( auto n=0; n<traits_.istate.size() && bret; n++ ) { 
    idstate = traits_.istate[n];
    
    // Set from initialized State vector, 
    for ( auto j=0; j<igbdy->size(); j++ ) {
      ind = (*igbdy)[j];
      (*ub[n])[j] = (*unew[idstate])[ind];
    }
  }
  bcomputed_ = true;

  return bret;

} // end of method compute_bdy_data


