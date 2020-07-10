//==================================================================================
// Module       : gutils.ipp
// Date         : 1/31/19 (DLR)
// Description  : GeoFLOW utilities namespace
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

namespace geoflow
{


//**********************************************************************************
//**********************************************************************************
// METHOD : smooth
// DESC   :
//          
// DESC   : Computes a weighted average.
//              Calculates:
//                u <- DSS(M_L u) / DSS(M_L),
//          where u is the field expressed in local format;
//          M_L is the local mass matrix (unassembled);
//          DSS is the application of Q Q^T, or the direct stiffness operator.
// ARGS   : 
//          grid : GGrid object
//          tmp  : tmp space 
//          op   : GGFX_OP 
//          u    : Locally expressed field to smooth
// RETURNS: none.
//************************************************************************************
template<typename T>
void smooth(GGrid &grid, GGFX_OP op, GTVector<T> &tmp, GTVector<T> &u)
{
  GGFX<T> *ggfx = &grid.get_ggfx();
  tmp = u;
 
  u.pointProd(*(grid.massop().data()));
  tmp = *(grid.imassop().data());
  ggfx->doOp(tmp, op);  // DSS(mass_local)

  u.pointProd(tmp);      // M_assembled u M_L

} // end, smooth method

//**********************************************************************************
//**********************************************************************************
// METHOD : smooth
// DESC   :
//          
// DESC   : Computes a weighted average.
//              Calculates:
//                u <- DSS(M_L u) / DSS(M_L),
//          where u is the field expressed in local format;
//          M_L is the local mass matrix (unassembled);
//          DSS is the application of Q Q^T, or the direct stiffness operator.
// ARGS   : 
//          ptree : main prop tree
//          xmin  : vector with coord minima, allocated here if necessary
//          xmax  : vector with coord maxima, allocated here if necessary
// RETURNS: none.
//************************************************************************************
template<typename T>
void coord_dims(const ropertyTree &ptree, GTVector<T> &xmin, GTVector<T> &xmax)
{
  GTPoint<T>   P0, P1, dP;
  std::vector<GFTYPE> fvec;
  GString      sgrid;
  PropertyTree gtree;

  sgrid = ptree.getValue<GString>("grid_type");
  if      ( "grid_icos"    == sgrid  ) {
    P0.resize(3); P1.resize(3); dP.resize(3);
    xmin.resize(3); xmax.resize(3);
    assert(GDIM == 2 && "GDIM must be 2");
    xmin[0] = xmax[0] = gridptree.getValue<GFTYPE>("radius");
    xmin[1] = -PI/2.0; xmax[1] = PI/2.0;
    xmin[2] = 0.0    ; xmax[2] = 2.0*PI;
  }
  if      ( "grid_sphere"  == sgrid ) {
    P0.resize(3); P1.resize(3); dP.resize(3);
    xmin.resize(3); xmax.resize(3);
    assert(GDIM == 3 && "GDIM must be 3");
    std::vector<GINT> ne(3);
    xmin[0] = gridptree.getValue<GFTYPE>("radiusi");
    xmax[0] = gridptree.getValue<GFTYPE>("radiuso");
    xmin[1] = -PI/2.0; xmax[1] = PI/2.0; // lat
    xmin[2] = 0.0    ; xmax[2] = 2.0*PI; // long
  }
  else if ( "grid_box"     == sgrid ) {
    P0.resize(GDIM); P1.resize(GDIM); dP.resize(GDIM);
    xmin.resize(GDIM); xmax.resize(GDIM);

    fvec = gridptree.getArray<GFTYPE>("xyz0");
    for ( auto j=0; j<GDIM; j++ ) P0[j] = fvec[j];
    fvec = gridptree.getArray<GFTYPE>("delxyz");
    for ( auto j=0; j<GDIM; j++ ) dP[j] = fvec[j];
    P1   = dP + P0;
    for ( autoj=0; j<GDIM; j++ ) {
      xmin[j] = P0[j];
      xmax[j] = P1[j];
    }
  }
  else {
    assert(FALSE);
  }

  
 
  u.pointProd(*(grid.massop().data()));
  tmp = *(grid.imassop().data());
  ggfx->doOp(tmp, op);  // DSS(mass_local)

  u.pointProd(tmp);      // M_assembled u M_L

} // end, coord_dims method


} // end, namespace

