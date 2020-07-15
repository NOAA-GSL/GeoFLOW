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
template<typename TypePack>
GDirichletBdy<TypePack>::GDirichletBdy(GDirichletBdy<TypePack>::Traits &traits) :
UpdateBdyBase<TypePack>(),
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
template<typename TypePack>
GDirichletBdy<TypePack>::~GDirichletBdy()
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
GBOOL GDirichletBdy<TypePack>::update_impl(
                              Grid      &grid,
                              StateInfo &stinfo,
                              Time      &time,
                              State     &utmp,
                              State     &u,
                              State     &ub)
{
   GString    serr = "GDirichletBdy<TypePack>::update_impl: ";
   GINT       idstate, ind;

  GTVector<GTVector<GSIZET>> *igbdy = &grid.igbdy_binned()[traits_.bdyid];

  if ( traits_.compute_once && bcomputed_ ) return TRUE;


  // Set boundary vector to corresp. value:
  for ( auto k=0; k<traits_.istate.size(); k++ ) { 
    idstate = traits_.istate[n];
    if ( stinfo.compdesc[idstate] == GSC_PRESCRIBED
      || stinfo.compdesc[idstate] == GSC_NONE ) continue;
    }
    // Set from initialized State vector, 
    for ( auto j=0; j<(*igbdy)[GBDY_DIRICHLET].size()
       && ub[idstate] != NULLPTR; j++ ) {
      ind = (*igbdy)[GBDY_DIRICHLET][j];
      (*ub[idstate])[j] = traits_.value[k];
    }
  }
  bcomputed = TRUE;

  return TRUE;

} // end of method update_impl


