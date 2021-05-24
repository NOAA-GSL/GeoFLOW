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
// GBdyType   itype;
   GINT       idd, k, nv;
   GSIZET     il, iloc, ind;
   Ftype      sum, tiny, us, ut;
   Ftype      xs, ys, zs;
   Ftype      xt, yt, zt;
   GTPoint<Ftype>             n, s, pu, t;
   VVecFtype                  *bdyNormals;
   GTVector<VVecFtype>        *bdyTangents;
// GTVector<GTVector<Ftype>>  *xnodes;

   tiny  = 100.0*std::numeric_limits<Ftype>::epsilon();

  // Handle 0-Flux bdy conditions. This
  // is enforced by solving
  //    vec{n} \cdot vec{u} = 0.
  //
  // Note: The following isn't good for vectorization;
  //       will revisit later.
  
  bdyNormals = &grid.bdyNormals();  // bdy normal vector
  bdyTangents= &grid.bdyTangents(); // bdy tangent plane vectors
  
  nv = grid.gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;

  for ( auto j=0; j<traits_.ibdyloc.size(); j++ ) {
    iloc = traits_.ibdyloc[j];      // index into bdy arrays
    ind  = traits_.ibdyvol[j];      // index into volume array

    n .assign(*bdyNormals,j);
    t .assign((*bdyTangents)[0],j);
    if ( GDIM > 2 ) {
      s .assign((*bdyTangents)[1],j);
    }
    pu.assign(u, nv, ind);

#if defined(_G_IS2D)

    ut           = pu.dot(t);  // u.tangent: tangential vel.
    zt           = 1.0/(t.x1*n.x2 - t.x2*n.x1);
    (*u[0])[ind] =  n.x2*ut * zt;
    (*u[1])[ind] = -n.x1*ut * zt;
  
#elif defined(_G_IS3D)

    ut           = pu.dot(t);  // u.tangent1: tangential vel.
    us           = pu.dot(s);  // u.tangent2: tangential vel.
    (*u[2])[ind] =  t.x3*ut + s.x3*us;
    (*u[1])[ind] =  s.x3 > tiny 
                 ?    n.x1*ut + s.x2*(*u[2])[ind]
                 :    (*u[1])[ind];
    (*u[0])[ind] =  n.x1 > tiny 
                 ?    -n.x2*(*u[0])[ind] - n.x3*(*u[2])[ind]
                 :    (*u[0])[ind];

#endif

  }


  return TRUE;

} // end of method update_impl


