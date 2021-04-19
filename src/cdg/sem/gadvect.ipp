//==================================================================================
// Module       : gadvect.hpp
// Date         : 11/11/18 (DLR)
// Description  : Represents the SEM discretization of the advection operator:
//                u.Grad p  This is a nonlinear operator, so should not derive 
//                from GLinOp. This operator requires that grid consist of
//                elements of only one type.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved.
// Derived From : none
//==================================================================================

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Constructor accepting grid operator only
// ARGS   : grid  : Grid object
// RETURNS: none
//**********************************************************************************
template<typename TypePack>
GAdvect<TypePack>::GAdvect(Grid &grid)
:
grid_           (&grid)
{
} // end of constructor method (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename TypePack>
GAdvect<TypePack>::~GAdvect()
{
} // end, destructor


//**********************************************************************************
//**********************************************************************************
// METHOD : apply
// DESC   : Compute application of this operator to input vector:
//            po = u.Grad p
//          NOTE: Require GDIM velocity components for GE_DEFORMED or
//                GE_REGULAR elements, and require GDIM+1 components
//                for GE_2DEMBEDDED elements.
// ARGS   : p   : input p field
//          u   : input vector field
//          utmp: tmp arrays, sie >= 2
//          po  : output (result) vector
//             
// RETURNS:  none
//**********************************************************************************
template<typename TypePack>
void GAdvect<TypePack>::apply(StateComp &p, const State &u, State &utmp, StateComp  &po) 
{
    
  GSIZET nxy = grid_->gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;

  assert(u.size() >= nxy && "Insufficient number of velocity components");

  typename Grid::BinnedBdyIndex *igb = &grid_->igbdy_binned();

  if ( p.size() <= 1 ) { // p is a constant 
    po = 0.0;
    return;
  }

  if ( u[0] != NULLPTR ) {
    grid_->deriv(p, 1, *utmp[0], po);
#if defined(GEOFLOW_USE_NEUMANN_HACK)
GMTK::zero<Ftype>(po,(*igb)[1][GBDY_0FLUX]);
GMTK::zero<Ftype>(po,(*igb)[3][GBDY_0FLUX]);
#endif
    po.pointProd(*u[0]); // do u_1 * dp/dx_1)
  }
  else {
    po = 0.0;
  }
  for ( auto j=1; j<u.size(); j++ ) { 
    if ( u[j] == NULLPTR ) continue;
    grid_->deriv(p, j+1, *utmp[1], *utmp[0]);
#if defined(GEOFLOW_USE_NEUMANN_HACK)
GMTK::zero<Ftype>(*utmp[0],(*igb)[0][GBDY_0FLUX]);
GMTK::zero<Ftype>(*utmp[0],(*igb)[2][GBDY_0FLUX]);
#endif
    utmp[0]->pointProd(*u[j]);
    po += *utmp[0];
  }
  po.pointProd(*grid_->massop().data()); // multiply by mass

} // end of method apply


