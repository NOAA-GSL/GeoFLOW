//==================================================================================
// Module       : gdirichlet_bdy.hpp
// Date         : 7/5/20 (DLR)
// Description  : Object that handles the the time updates for
//                Dirichlet boundaries
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
GDirichletBdy<Types>::GDirichletBdy(typename GDirichletBdy<Types>::Traits &traits) :
UpdateBdyBase<Types>(),
bcomputed_               (FALSE),
traits_                 (traits)
{

  assert( traits_.value.size() == traits_.istate.size() 
       && "State info structure invalid");

} // end of constructor method 


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename Types>
GDirichletBdy<Types>::~GDirichletBdy()
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
GBOOL GDirichletBdy<Types>::update_impl(
                              Grid      &grid,
                              StateInfo &stinfo,
                              Time      &time,
                              State     &utmp,
                              State     &u,
                              State     &ub)
{
   GString    serr = "GDirichletBdy<Types>::update_impl: ";
   GINT       idstate;
   GSIZET     ind;

  GTVector<GSIZET> *igbdy = &traits_.ibdyvol;

  if ( traits_.compute_once && bcomputed_ ) return TRUE;


  // Set boundary vector to corresp. value:
  for ( auto k=0; k<traits_.istate.size(); k++ ) { 
    idstate = traits_.istate[k];
    if ( stinfo.icomptype[idstate] == GSC_PRESCRIBED
      || stinfo.icomptype[idstate] == GSC_NONE ) continue;
    
    // Set from initialized State vector, 
    for ( auto j=0; j<igbdy->size()
       && ub[idstate] != NULLPTR; j++ ) {
      ind = (*igbdy)[j];
//    (*ub[idstate])[j] = traits_.value[k];
      (*ub[idstate])[ind] = traits_.value[k];
    }
  }
  bcomputed_ = TRUE;

  return TRUE;

} // end of method update_impl


