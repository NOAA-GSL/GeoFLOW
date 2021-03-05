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


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_temp
// DESC   : Compute temperature from state
//             T =  e /( d * Cv ),
//          with e the sensible internal energy density, 
//          d = total density, and Cv is, e.g.
//             Cv = Cvd qd + Cvv qv + Sum_i(Cl_i ql_i) + Sum_j(Ci_j qi_j).
// ARGS   : e    : energy density
//          d    : density
//          cv   : specific heat
//          temp : temperature
// RETURNS: none.
//**********************************************************************************
template<typename T>
void compute_temp(const GTVector<T> &e, const GTVector<T> &d, const GTVector<T> &cv,  GTVector<T> &temp)
{
   GString    serr = "compute_temp: ";

   // Compute temperature:
   for ( auto j=0; j<e.size(); j++ ) {
     temp[j] = e[j] / ( d[j] * cv[j] );
   }

} // end of method compute_temp


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_p 
// DESC   : Compute total pressure from state
//              p = d ( q R) T,
//          with total density, d, q the
//          dry mass fraction, R the gas constants, 
//          and T the temperature.
// ARGS   : Temp: temperature
//          d: density
//          q: mass fraction
//          R: gas constant
//          p: pressure fluctuation field returned
// RETURNS: none.
//**********************************************************************************
template<typename T>
void compute_p(const GTVector<T> &Temp, const GTVector<T> &d, const GTVector<T> &q, GFTYPE R, GTVector<T> &p)
{
   GString    serr = "compute_p: ";

   // p' = d q R T:
   for ( auto j=0; j<p.size(); j++ ) {
     p[j] = d[j] * q[j] * R * Temp[j];
   }

} // end of method compute_p 


//**********************************************************************************
//**********************************************************************************
// METHOD : in_seg (1)
// DESC   : Determine if points, x, lie within segment 
//          defined by verts. Return indices of contained x
//          in ind vector, and return the number of points found
// ARGS   : verts: use first 2 points to define segment
//          x    : test points
//          ind  : vector of indices in x that lie on segment. Reallocated
//                 if necessary to be of insufficient size
// RETURNS: number of points found on segment (valid number of points in ind)
//**********************************************************************************
template<typename T>
GSIZET in_seg(const GTVector<GTPoint<T>> &verts, const GTVector<GTVector<T>> &x, GTVector<GSIZET> &ind)
{
  GINT       ndim;
  GSIZET     nfound;
  T          mr1, mr2;
  GTPoint<T> c(3), d(3), p(3), r1(3), r2(3);

  p    = 0.0;
  d    = verts[1] ; d -= verts[0];
  ndim = x.size();
  ind.resizem(x[0].size());
  nfound = 0;
  for ( auto i=0; i<x[0].size(); i++ ) {
    p.assign(x, i);
    r1 = p ; r1 -= verts[0]; mr1 = r1.mag();
    r2 = p ; r2 -= verts[1]; mr2 = r2.mag();
    d.cross(r1, c); // c = d X r1
    if ( (mr1+mr2) <= d.mag() ) {
      ind[nfound++] = i;
    }
  } // end, test point loop
  
  return nfound;

} // end of method in_seg (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : in_seg (2)
// DESC   : Determine if points, x, lie within segment 
//          defined by verts. Return indices of contained x
//          in ind vector, and return the number of points found
// ARGS   : verts: use first 2 points to define segment
//          x    : test points
//          ix   : indirection indices in x to consider
//          ind  : vector of indices in x that lie on segment. Reallocated
//                 if necessary to be of insufficient size
// RETURNS: number of points found on segment (valid number of points in ind)
//**********************************************************************************
template<typename T>
GSIZET in_seg(const GTVector<GTPoint<T>> &verts, const GTVector<GTVector<T>> &x, const GTVector<GSIZET> &ix, GTVector<GSIZET> &ind)
{
  GINT       ndim;
  GSIZET     nfound;
  T          mr1, mr2;
  GTPoint<T> c(3), d(3), p(3), r1(3), r2(3);

  p    = 0.0;
  d    = verts[1] ; d -= verts[0];
  ndim = x.size();
  ind.resizem(x[0].size());
  nfound = 0;
  for ( auto i=0; i<ix.size(); i++ ) {
    p.assign(x, ix[i]);
    r1 = p ; r1 -= verts[0]; mr1 = r1.mag();
    r2 = p ; r2 -= verts[1]; mr2 = r2.mag();
    d.cross(r1, c); // c = d X r1
    if ( (mr1+mr2) <= d.mag() ) {
      ind[nfound++] = ix[i];
    }
  } // end, test point loop
  
  return nfound;

} // end of method in_seg (2)


//**********************************************************************************
//**********************************************************************************

// METHOD : in_poly (1)
// DESC   : Determine if points, x, lie within polygon 
//          defined by verts. Return indices of enclosed x
//          in ind vector, and return the number of points found
// ARGS   : verts: vertices defining polygon, in 'counter-clockwise' order
//          x    : test points (x, y, ... z arrays)
//          ind  : vector of indices in x that lie on segment. Reallocated
//                 if necessary to be of insufficient size
//          eps  : 'zero' value
// RETURNS: number of points found in polygon (valid number of points in ind)
//**********************************************************************************
template<typename T>
GSIZET in_poly(GTVector<GTPoint<T>> &verts, const GTVector<GTVector<T>> &x, T eps, GTVector<GSIZET> &ind)
{
  GINT         nverts;
  GSIZET       nfound;
  T            dotdrN, dotdrc;
  GTPoint<T>   c, N, dr, r1, r2, xp;

  ind.resizem(x[0].size());
  nverts = verts.size();
  nfound = 0;
  for ( auto i=0; i<x[0].size(); i++ ) {
    for ( auto n=0; n<verts.size(); n++ ) {

      // Assume vertices in counterclockwise order, and
      // compute normal to edge segment n, np1. This vector
      // points 'inward' in a r.h. sense:
      r1 = verts[(n+1)%nverts]; r1 -= verts[n];
      r2 = verts[(n+2)%nverts]; r2  -= verts[(n+1)%nverts];
     
      r1.cross(r2, c); // c = r1 X r2, normal to plane

      c.cross(r1, N);  // N = c X r1, in plane 

      xp.assign(x, i); // test point

      // Compute vector dr =  r_test - polyVertex:
      dr = xp; dr -= verts[n];
      dotdrN = dr.dot(N);
      dotdrc = dr.dot(c);

      // If dr dot N < 0, or there is a component perp to plane of
      // polygon (dotPar !=0 ) then point is outside of polygon:
      if ( (fabs(dotdrN) > eps && dotdrN < 0) 
        ||  fabs(dotdrc) > eps ) {
        ind[nfound++] = i;
      }
    } // end, vert loop
  }   // end, test point loop

  return nfound;   

} // end of method in_poly (1)


//**********************************************************************************
//**********************************************************************************

// METHOD : in_poly (2)
// DESC   : Determine if points, x, lie within polygon 
//          defined by verts. Return indices of enclosed x
//          in ind vector, and return the number of points found
// ARGS   : verts: vertices defining polygon, in 'counter-clockwise' order
//          x    : test points (x, y, ... z arrays)
//          ix   : indirection indices in x to consider
//          ind  : vector of indices in x that lie on segment. Reallocated
//                 if necessary to be of insufficient size
//          eps  : 'zero' value
// RETURNS: number of points found in polygon (valid number of points in ind)
//**********************************************************************************
template<typename T>
GSIZET in_poly(GTVector<GTPoint<T>> &verts, const GTVector<GTVector<T>> &x, const GTVector<GSIZET> &ix, T eps, GTVector<GSIZET> &ind)
{
  GINT         nverts;
  GSIZET       nfound;
  T            dotdrN, dotdrc;
  GTPoint<T>   c, N, dr, r1, r2, xp;

  ind.resizem(x[0].size());
  nverts = verts.size();
  nfound = 0;
  for ( auto i=0; i<ix.size(); i++ ) {
    for ( auto n=0; n<verts.size(); n++ ) {

      // Assume vertices in counterclockwise order, and
      // compute normal to edge segment n, np1. This vector
      // points 'inward' in a r.h. sense:
      r1 = verts[(n+1)%nverts]; r1 -= verts[n];
      r2 = verts[(n+2)%nverts]; r2  -= verts[(n+1)%nverts];
     
      r1.cross(r2, c); // c = r1 X r2, normal to plane

      c.cross(r1, N);  // N = c X r1, in plane 

      xp.assign(x, ix[i]); // test point

      // Compute vector dr =  r_test - polyVertex:
      dr = xp; dr -= verts[n];
      dotdrN = dr.dot(N);
      dotdrc = dr.dot(c);

      // If dr dot N < 0, or there is a component perp to plane of
      // polygon (dotPar !=0 ) then point is outside of polygon:
      if ( (fabs(dotdrN) > eps && dotdrN < 0) 
        ||  fabs(dotdrc) > eps ) {
        ind[nfound++] = ix[i];
      }
    } // end, vert loop
  }   // end, test point loop

  return nfound;   

} // end of method in_poly (2)


} // end, namespace

