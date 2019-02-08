//==================================================================================
// Module       : init_heat.cpp
// Date         : 2/4/19 (DLR)
// Description  : Initializes heat equation with Gaussian 'bump'
// Copyright    : Copyright 2019. Colorado State University. All rights reserved.
// Derived From :
//==================================================================================

//**********************************************************************************
//**********************************************************************************
// ARGS   : grid  : GGrid object
//          t     : time at which solution is computed
//          ptree : master property tree
//          sblock: name of property tree object that contains 
//                  parameters required for this function
// RETURNS: none
//**********************************************************************************

void init_pde(GGrid &grid, Time &t, const tbox::PropertyTree& ptree, GString &sblock, State &u);
{

  PropertyTree initptree;
  
  initptree   = ptree.getPropertyTree(sblock);

  // Get initialization parameters from sblock:
  GFTYPE r2;
  GFTYPE x0    = initptree.getValue("x0");
  GFTYPE y0    = initptree.getValue("y0");
  GFTYPE z0    = initptree.getValue("z0");
  GFTYPE u0    = initptree.getValue("u0");
  GFTYPE sigma = initptree.getValue("sigma");

  GTVector<GFTYPE> *x;
  GTVector<GFTYPE> *y;
  GTVector<GFTYPE> *z;

  GSISET i, j, k, n;
#if defined(_G_IS2D)
  if ( grid.gtype() == GE_2DEMBEDDED ) { 
    x = &grid.xNodes(0);
    y = &grid.xNodes(1);
    z = &grid.xNodes(2);
    for ( j=0, n=0; j<x->size(); j++ ) {
      for ( i=0; i<x->size(); i++, n++ ) {
         r2 = pow((*x)[n] - x0,2) 
            + pow((*y)[n] - y0,2)
            + pow((*z)[n] - z0,2) ;
         u[0][n] = u0*exp(-r2/(sigma*sigma)); 
      }
    }
  }
  else { // GE_REGULAR
    x = &grid.xNodes(0);
    y = &grid.xNodes(1);
    for ( j=0, n=0; j<x->size(); j++ ) {
      for ( i=0; i<x->size(); i++, n++ ) {
         r2 = pow((*x)[n] - x0,2) 
            + pow((*y)[n] - y0,2);
         u[0][n] = u0*exp(-r2/(sigma*sigma)); 
      }
    }
  }
#elif defined(_G_IS3D)
    x = &grid.xNodes(0);
    y = &grid.xNodes(1);
    z = &grid.xNodes(2);
    for ( k=0, n=0; k<x->size(); k++ ) {
      for ( j=0, n=0; j<x->size(); j++ ) {
        for ( i=0; i<x->size(); i++, n++ ) {
          r2 = pow((*x)[n] - x0,2) 
             + pow((*y)[n] - y0,2)
             + pow((*z)[n] - z0,2) ;
          u[0][n] = u0*exp(-r2/(sigma*sigma)); 
        }
      }
    }
#endif

} // end, init_heat

