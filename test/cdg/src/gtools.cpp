//==================================================================================
// Module       : gtools.cpp
// Date         : 5/5/19 (DLR)
// Description  : GeoFLOW initialization tools
// Copyright    : Copyright 2019. Colorado State University. All rights reserved.
// Derived From : 
//==================================================================================
#include "gtools.h"

//**********************************************************************************
//**********************************************************************************
// METHOD: init_ggfx
// DESC  : Initialize gather/scatter operator
// ARGS  : ptree   : main property tree
//         grid    : GGrid object
//         ggfx    : gather/scatter op, GGFX
//**********************************************************************************
void init_ggfx(PropertyTree &ptree, GGrid &grid, GGFX<GFTYPE> &ggfx)
{
  GFTYPE                         delta;
  GFTYPE                         rad;
  GMorton_KeyGen<GNODEID,GFTYPE> gmorton;
  GTPoint<GFTYPE>                dX, porigin, P0;
  GTVector<GNODEID>              glob_indices;
  GTVector<GTVector<GFTYPE>>    *xnodes;
  GString                        sgrid;
  std::vector<GFTYPE>            pstd;
  PropertyTree                   gtree;

  sgrid = ptree.getValue<GString>("grid_type");
  gtree = ptree.getPropertyTree(sgrid);

  if ( sgrid == "grid_box" ) {
    pstd = gtree.getArray<GFTYPE>("xyz0");
    P0   = std;
  }
  if ( sgrid == "grid_icos" ) {
    rad   = gtree.getValue<GFTYPE>("radius");
    P0.x1 = -rad ;
    P0.x2 = -rad ;
    P0.x3 = -rad ;
  }
  if ( sgrid == "grid_sphere" ) {
    rad   = gtree.getValue<GFTYPE>("radiuso");
    P0.x1 = -rad ;
    P0.x2 = -rad ;
    P0.x3 = -rad ;
  }

  // First, periodize coords if required to, 
  // before labeling nodes:
  if ( typeid(*grid_) == typeid(GGridBox) ) { 
    static_cast<GGridBox*>(grid_)->periodize();
  }

  delta  = grid.minnodedist();
  dX     = 0.1*delta;
  xnodes = &grid.xNodes();
  glob_indices.resize(grid.ndof());

  // Integralize *all* internal nodes
  // using Morton indices:
  gmorton.setIntegralLen(P0,dX);

  gmorton.key(glob_indices, *xnodes);
  GComm::Synch(comm_);

#if 0
  GPP(comm_,"init_ggfx: glob_indices=" << glob_indices);
#endif

  // Initialize gather/scatter operator:
  GBOOL bret;
  bret = ggfx.init(glob_indices);
  assert(bret && "Initialization of GGFX operator failed");

  // Unperiodize nodes now that connectivity map is
  // generated, so that coordinates mean what they should:
  if ( typeid(*grid_) == typeid(GGridBox) ) { // periodize coords
    static_cast<GGridBox*>(grid_)->unperiodize();
  }

} // end method init_ggfx

