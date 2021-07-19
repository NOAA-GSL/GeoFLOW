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
//          time  : timestep
//          utmp  : tmp vectors
//          u     : state vector
// RETURNS: none.
//**********************************************************************************
template<typename Types>
GBOOL GDirichletBdy<Types>::update_impl(
                              EqnBasePtr &eqn,
                              Grid       &grid,
                              Time       &time,
                              State      &utmp,
                              State      &u)
{
   GString    serr = "GDirichletBdy<Types>::update_impl: ";
   GINT       idstate;
   GSIZET     ind;

   vector<GSIZET> *igbdy = &traits_.ibdyvol;



  // Set boundary vector to corresp. value:
  for ( auto k=0; k<traits_.istate.size(); k++ ) { 
    idstate = traits_.istate[k];
    
    // Set from initialized State vector, 
    for ( auto j=0; j<igbdy->size(); j++ ) {
      ind = (*igbdy)[j];
      (*u[idstate])[ind] = traits_.value[k];
    }
  }

  return TRUE;

} // end of method update_impl


