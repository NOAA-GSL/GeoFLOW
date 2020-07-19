//==================================================================================
// Module       : gutils.ipp
// Date         : 1/31/19 (DLR)
// Description  : GeoFLOW utilities namespace
// Copyright    : Copyright 2019. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================

namespace geoflow
{


//**********************************************************************************
//**********************************************************************************
// METHOD : append
// DESC   : Appends add vector to base vector, modifying base
// ARGS   : 
//          base   : input vector to be modified
//          add    : vector to append to base
// RETURNS: none.
//************************************************************************************
template<typename T>
void append(GTVector<T> &base, GTVector<T> &add)
{
  GSIZET sznew;
  GTVector<T> tmp(base.size());

  sznew = base.size() + add.size();
  for ( auto i=0; i<base.size(); i++ ) {
    tmp[i] = base[i];
  }
  for ( auto i=base.size(); i<sznew; i++ ) {
    tmp[i] = add[i-base.size()];
  }

  base.resize(sznew);
  base = tmp;


} // end, unique method

//**********************************************************************************
//**********************************************************************************
// METHOD : unique
// DESC   : Finds indices of unique variables in input vector
// ARGS   : 
//          vec    : input vector
//          ibeg,
//          iend   : begining, ending indices to search in vec
//          iunique: beginning indices for unique elements. Is resized to
//                   exact amount to account for the number of unique
//                   elements found.
// RETURNS: none.
//************************************************************************************
template<typename T>
void unique(GTVector<T> &vec, GSIZET ibeg, GSIZET iend, GTVector<GSIZET> &iunique)
{
  GSIZET n=0;
  T      val;
  GTVector<GSIZET> tmp(vec.size());

  for ( auto i=ibeg; i<iend+1; i+=vec.gindex_.stride() ) {
    val = vec[i];
    if ( !iunique.containsn(val,n) ) {
      tmp[n] = val;
      n++;
    }
  }

  iunique.resize(n);
  for ( auto i=0; i<n; i++ ) iunique[i] = tmp[i];

} // end, unique method


//**********************************************************************************
//**********************************************************************************
// METHOD : coord_dims
// DESC   : Gets and or computes coord dimensions from 
//          prop tree
// ARGS   : 
//          ptree : main prop tree
//          xmin  : vector with coord minima, allocated here if necessary
//          xmax  : vector with coord maxima, allocated here if necessary
// RETURNS: none.
//************************************************************************************
template<typename T>
void coord_dims(const geoflow::tbox::PropertyTree &ptree, GTVector<T> &xmin, GTVector<T> &xmax)
{
  GTPoint<T>   P0, P1, dP;
  std::vector<GFTYPE> fvec;
  GString      sgrid;
  geoflow::tbox::PropertyTree gtree;

  sgrid = ptree.getValue<GString>("grid_type");
  if      ( "grid_icos"    == sgrid  ) {
    P0.resize(3); P1.resize(3); dP.resize(3);
    xmin.resize(3); xmax.resize(3);
    assert(GDIM == 2 && "GDIM must be 2");
    xmin[0] = xmax[0] = gtree.getValue<GFTYPE>("radius");
    xmin[1] = -PI/2.0; xmax[1] = PI/2.0;
    xmin[2] = 0.0    ; xmax[2] = 2.0*PI;
  }
  if      ( "grid_sphere"  == sgrid ) {
    P0.resize(3); P1.resize(3); dP.resize(3);
    xmin.resize(3); xmax.resize(3);
    assert(GDIM == 3 && "GDIM must be 3");
    std::vector<GINT> ne(3);
    xmin[0] = gtree.getValue<GFTYPE>("radiusi");
    xmax[0] = gtree.getValue<GFTYPE>("radiuso");
    xmin[1] = -PI/2.0; xmax[1] = PI/2.0; // lat
    xmin[2] = 0.0    ; xmax[2] = 2.0*PI; // long
  }
  else if ( "grid_box"     == sgrid ) {
    P0.resize(GDIM); P1.resize(GDIM); dP.resize(GDIM);
    xmin.resize(GDIM); xmax.resize(GDIM);

    fvec = gtree.getArray<GFTYPE>("xyz0");
    for ( auto j=0; j<GDIM; j++ ) P0[j] = fvec[j];
    fvec = gtree.getArray<GFTYPE>("delxyz");
    for ( auto j=0; j<GDIM; j++ ) dP[j] = fvec[j];
    P1   = dP + P0;
    for ( auto j=0; j<GDIM; j++ ) {
      xmin[j] = P0[j];
      xmax[j] = P1[j];
    }
  }
  else {
    assert(FALSE);
  }

} // end, coord_dims method


} // end, namespace

