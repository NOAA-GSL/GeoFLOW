//==================================================================================
// Module       : gnoslip_bdy.hpp
// Date         : 7/5/20 (DLR)
// Description  : Object that handles the the time updates for
//                no-slip boundary conditions. Acts on kinetic
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
GNoSlipBdy<TypePack>::GNoSlipBdy(GNoSlipBdy<TypePack>::Traits &traits) :
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
GNoSlipBdy<TypePack>::~GNoSlipBdy()
{
} // end, destructor


//**********************************************************************************
//**********************************************************************************
// METHOD : update_impl
// DESC   : Entry method for hahdling no-slip bdyconditions
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
GBOOL GNoSlipBdy<TypePack>::update_impl(
                              Grid      &grid,
                              StateInfo &stinfo,
                              Time      &time,
                              State     &utmp,
                              State     &u,
                              State     &ub)
{
  GString    serr = "GNoSlipBdy<TypePack>::update_impl: ";

  GTVector<GTVector<GSIZET>> *igbdy = &grid.igbdy_binned()[traits_.bdyid];

  if ( traits_.compute_once && bcomputed_ ) return TRUE;

  for ( auto k=0; k<nstate_; k++ ) { // for each vector component
    for ( auto j=0; j<(*igbdy)[GBYD_0FLUX].size() ) { // all bdy points
      ind = (*igbdy)[itype][j]; // index into long vector array
      (*ub[k])[ind] = 0.0;
    }
  }

  bcomputed = TRUE;

  return TRUE;

} // end of method update_impl


