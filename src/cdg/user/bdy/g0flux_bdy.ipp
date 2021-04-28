//==================================================================================
// Module       : g0flux_bdy.hpp
// Date         : 7/5/20 (DLR)
// Description  : Object that handles the the time updates for
//                0-flux boundaries. Acts on kinetic
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
template<typename Types>
G0FluxBdy<Types>::G0FluxBdy(typename G0FluxBdy<Types>::Traits &traits) :
UpdateBdyBase<Types>(),
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
template<typename Types>
G0FluxBdy<Types>::~G0FluxBdy()
{
} // end, destructor


//**********************************************************************************
//**********************************************************************************
// METHOD : update_impl
// DESC   : Entry method for updating 0-Flux bdy conditions
// ARGS   : 
//          eqn   : equation implementation
//          grid  : grid object (necessary?)
//          time  : timestep
//          utmp  : tmp vectors
//          u     : state vector
// RETURNS: none.
//**********************************************************************************
template<typename Types>
GBOOL G0FluxBdy<Types>::update_impl(
                              EqnBasePtr &eqn,
                              Grid       &grid,
                              Time       &time,
                              State      &utmp,
                              State      &u)
{
   GString    serr = "G0FluxBdy<Types>::update_impl: ";
   GBdyType   itype;
   GINT       idd, k;
   GSIZET     il, iloc, ind;
   Ftype      sum, xn;
   GTVector<GTVector<Ftype>>  *bdyNormals;
// GTVector<GTVector<Ftype>>  *xnodes;
   GTVector<GINT>             *idep;




  // Handle 0-Flux bdy conditions. This
  // is computed by solving
  //    vec{n} \cdot vec{u} = 0
  // for 'dependent' component set in grid.
  //
  // Note: We may want to switch the order of the
  //       following loops to have a better chance
  //       of vectorization. Unrolling likely
  //       won't occur:
  bdyNormals = &grid.bdyNormals(); // bdy normal vector
  idep       = &grid.idepComp();   // dependent components
  itype      = GBDY_0FLUX;
  for ( auto j=0; j<traits_.ibdyloc.size(); j++ ) {
    iloc = traits_.ibdyloc[j];     // index into bdy arrays
    ind  = traits_.ibdyvol[j];     // index into volume array
    idd  = (*idep)[iloc];          // dependent vector component


    xn   = (*bdyNormals)[idd][iloc];// n_idd == normal component for dependent vector comp
    sum  = 0.0;
    for ( auto k=0; k<traits_.istate.size(); k++ ) { // for each indep vector component
        sum += traits_.istate[k] != idd 
             ? (*bdyNormals)[k][iloc] * (*u[traits_.istate[k]])[ind]
             : 0.0;
    }
    (*u[idd])[ind] = -sum / xn; // Ensure v.n = 0:
   
    if ( GET_NDTYPE(traits_.ibdydsc[iloc]) == GElem_base::VERTEX ) {
      for ( auto k=0; k<traits_.istate.size(); k++ ) { 
        (*u[traits_.istate[k]])[ind] = 0.0;
      }
    }

  }


  return TRUE;

} // end of method update_impl


