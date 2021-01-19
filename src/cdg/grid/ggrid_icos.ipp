//==================================================================================
// Module       : ggrid_icos.ipp
// Date         : 8/31/18 (DLR)
// Description  : Object defining a (global) icosahedral grid, that in 2d
//                uses (extrinsic) gnomonic projections to locate element vertices.
//                Vertices always reside on sphere, so centroids will
//                not (will reside within). In 3d, the base is computed from
//                the same procedure as in 2d, but we use isoparameteric
//                representation on the sphere.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved.
// Derived From : GGrid.
//==================================================================================

#include "tbox/tracer.hpp"

//**********************************************************************************
//**********************************************************************************
// METHOD : project2sphere (1)
// DESC   : Project Cartesian coords to sphere, specified by rad argument, and
//          express in Cartesian coords.  Necessary for 2d grids.
// ARGS   : tmesh: Vector of triangles (vertices), modified to contain 
//                 projected coordinates
// RETURNS: none.
//**********************************************************************************
template<typename T>
void GGridIcos::project2sphere(GTVector<GTriangle<T>> &tmesh, T rad)
{
	GEOFLOW_TRACE();
  GString serr = "GridIcos::project2sphere (1): ";

  T          r, xlat, xlong;
  GTPoint<T> v;

  for ( GSIZET i=0; i<tmesh_.size(); i++ ) { // loop over all triangles in tmesh
    for ( GSIZET j=0; j<3; j++ ) { // guaranteed to have 3 points each
      v = *tmesh[i].v[j];
      r = v.norm();
      xlat  = asin(v.x3/r);
      xlong = atan2(v.x2,v.x1);
      xlong = xlong < 0.0 ? 2.0*PI+xlong : xlong;
      v.x1 = rad*cos(xlat)*cos(xlong);
      v.x2 = rad*cos(xlat)*sin(xlong);
      v.x3 = rad*sin(xlat);
      *tmesh[i].v[j] = v;
    }
  }

} // end of method project2sphere (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : project2sphere (2)
// DESC   : Project Cartesian coords to sphere, specified by rad argument.
//          Necessary for 2d grids.
// ARGS   : plist: Vector of points, modified to contain 
//                 projected coordinates. Must be 3d points.
// RETURNS: none.
//**********************************************************************************
template<typename T>
void GGridIcos::project2sphere(GTVector<GTPoint<T>> &plist, T rad)
{
	GEOFLOW_TRACE();
  GString serr = "GridIcos::project2sphere (2): ";

  T          r, xlat, xlong;
  GTPoint<T> v;

  for ( GSIZET i=0; i<plist.size(); i++ ) { // loop over all points
    v = plist[i];
    r = v.norm();
    xlat  = asin(v.x3/r);
    xlong = atan2(v.x2,v.x1);
    xlong = xlong < 0.0 ? 2.0*PI+xlong : xlong;
    v.x1 = rad*cos(xlat)*cos(xlong);
    v.x2 = rad*cos(xlat)*sin(xlong);
    v.x3 = rad*sin(xlat);
    plist[i] = v;
  }

} // end of method project2sphere (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : project2sphere (3)
// DESC   : Project Cartesian coords to sphere, specified by rad argument.
//          Necessary for 2d grids.
// ARGS   : plist: Vector of vectors, one vector for each dimension, modified to 
//                 contain projected coordinates. Must be 3d points.
// RETURNS: none.
//**********************************************************************************
template<typename T>
void GGridIcos::project2sphere(GTVector<GTVector<T>> &plist, T rad)
{
	GEOFLOW_TRACE();
  GString serr = "GridIcos::project2sphere (3): ";

  T r, xlat, xlong, x, y, z;

  for ( GSIZET i=0; i<plist[0].size(); i++ ) { // loop over all points
    x = plist[0][i]; y = plist[1][i]; z = plist[2][i];
    r = sqrt(x*x + y*y + z*z);
    xlat  = asin(z/r);
    xlong = atan2(y,x);
    xlong = xlong < 0.0 ? 2.0*PI+xlong : xlong;
    plist[0][i] = rad*cos(xlat)*cos(xlong);
    plist[1][i] = rad*cos(xlat)*sin(xlong);
    plist[2][i] = rad*sin(xlat);
  }

} // end of method project2sphere (3)


//**********************************************************************************
//**********************************************************************************
// METHOD : spherical2xyz (1)
// DESC   : Transform from spherical-polar to Cartesian coords
// ARGS   : plist: Vector of points, representing spherical coordinates 
//                 (r, lat, long), to be converted to (x, y, z), in-place
// RETURNS: none.
//**********************************************************************************
template<typename T>
void GGridIcos::spherical2xyz(GTVector<GTPoint<T>*> &plist)
{
	GEOFLOW_TRACE();
  GString serr = "GridIcos::spherical2xyz(1): ";

  T r, xlat, xlong;

  for ( GSIZET i=0; i<plist.size(); i++ ) { // loop over all points
    r = plist[i]->x1; xlat = plist[i]->x2; xlong = plist[i]->x3;
    plist[i]->x1 = r*cos(xlat)*cos(xlong);
    plist[i]->x2 = r*cos(xlat)*sin(xlong);
    plist[i]->x3 = r*sin(xlat);
  }

} // end of method spherical2xyz (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : spherical2xyz (2)
// DESC   : Transform from spherical-polar to Cartesian coords
// ARGS   : plist: Vector of points, representing spherical coordinates 
//                 (r, lat, long), to be converted to (x, y, z), in-place
// RETURNS: none.
//**********************************************************************************
template<typename T>
void GGridIcos::spherical2xyz(GTVector<GTPoint<T>> &plist)
{
	GEOFLOW_TRACE();
  GString serr = "GridIcos::spherical2xyz(2): ";

  T r, xlat, xlong;

  for ( GSIZET i=0; i<plist.size(); i++ ) { // loop over all points
    r = plist[i].x1; xlat = plist[i].x2; xlong = plist[i].x3;
    plist[i].x1 = r*cos(xlat)*cos(xlong);
    plist[i].x2 = r*cos(xlat)*sin(xlong);
    plist[i].x3 = r*sin(xlat);
  }

} // end of method spherical2xyz (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : spherical2xyz (3)
// DESC   : Transform from spherical-polar to Cartesian coords
// ARGS   : plist: Vector of points, representing spherical coordinates 
//                 (r, lat, long), to be converted to (x, y, z), in-place
// RETURNS: none.
//**********************************************************************************
template<typename T>
void GGridIcos::spherical2xyz(GTVector<GTVector<T>> &plist)
{
	GEOFLOW_TRACE();
  GString serr = "GridIcos::spherical2xyz(3): ";
  assert(plist.size() >= 3 && "Invalid dimensionality on input array");

  T r, xlat, xlong;

  for ( GSIZET i=0; i<plist[0].size(); i++ ) { // loop over all points
    r = plist[0][i]; xlat = plist[1][i]; xlong = plist[2][i];
    plist[0][i] = r*cos(xlat)*cos(xlong);
    plist[1][i] = r*cos(xlat)*sin(xlong);
    plist[2][i] = r*sin(xlat);
  }

} // end of method spherical2xyz (3)


//**********************************************************************************
//**********************************************************************************
// METHOD : xyz2spherical (1)
// DESC   : Transform from Cartesian coords to spherical-polar
// ARGS   : plist: Vector of points, representing Cartesian coordinates 
//                 (x,y,z) to be converted to (r, lat, long), in-place
// RETURNS: none.
//**********************************************************************************
template<typename T>
void GGridIcos::xyz2spherical(GTVector<GTPoint<T>*> &plist)
{
	GEOFLOW_TRACE();
  GString serr = "GridIcos::xyz2spherica(1): ";

  T r, x, y, z;

  for ( GSIZET i=0; i<plist.size(); i++ ) { // loop over all points
   x = plist[i]->x1; y = plist[i]->x2; z = plist[i]->x3;
   r = sqrt(x*x + y*y + z*z);
   plist[i]->x1 = r;
   plist[i]->x2 = asin(z/r);
   plist[i]->x3 = atan2(y,x);
   plist[i]->x3 = plist[i]->x3 < 0.0 ? 2.0*PI+plist[i]->x3 : plist[i]->x3;
  }

} // end of method xyz2spherical (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : xyz2spherical (2)
// DESC   : Transform from Cartesian coords to spherical-polar
// ARGS   : plist: Vector of points, representing Cartesian coordinates 
//                 (x,y,z) to be converted to (r, lat, long), in-place
// RETURNS: none.
//**********************************************************************************
template<typename T>
void GGridIcos::xyz2spherical(GTVector<GTVector<T>> &plist)
{
	GEOFLOW_TRACE();
  GString serr = "GridIcos::xyz2spherica(2): ";

  T r, x, y, z;

  for ( GSIZET i=0; i<plist.size(); i++ ) { // loop over all points
   x = plist[i][0]; y = plist[i][1]; z = plist[i][2];
   r = sqrt(x*x + y*y + z*z);
   plist[i][0] = r;
   plist[i][1] = asin(z/r);
   plist[i][2] = atan2(y,x);
   plist[i][2] = plist[i][2] < 0.0 ? 2.0*PI+plist[i][2] : plist[i][2];
  }

} // end of method xyz2spherical (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : xyz2spherical (3)
// DESC   : Transform from Cartesian coords to spherical-polar
// ARGS   : plist: Vector of points, representing Cartesian coordinates 
//                 (x,y,z) to be converted to (r, lat, long), in-place
// RETURNS: none.
//**********************************************************************************
template<typename T>
void GGridIcos::xyz2spherical(GTVector<GTPoint<T>> &plist)
{
	GEOFLOW_TRACE();
  GString serr = "GridIcos::xyz2spherica(3): ";

  T r, x, y, z;

  for ( GSIZET i=0; i<plist.size(); i++ ) { // loop over all points
   x = plist[i].x1; y = plist[i].x2; z = plist[i].x3;
   r = sqrt(x*x + y*y + z*z);
   plist[i].x1 = r;
   plist[i].x2 = asin(z/r);
   plist[i].x3 = atan2(y,x);
   plist[i].x3 = plist[i].x3 < 0.0 ? 2.0*PI+plist[i].x3 : plist[i].x3;
  }

} // end of method xyz2spherical (3)


//**********************************************************************************
//**********************************************************************************
// METHOD : cart2gnomonic
// DESC   : Transform Cartesian coords to gnomonic space
// ARGS   : clist: Vector of Cartesian points
//          rad  : radius of sphere
//          xlatc,
//          xlongc: lat and long of reference vector (usually a centroid)
//          glist: converted gnomonic points
//          
// RETURNS: none.
//**********************************************************************************
template<typename T>
void GGridIcos::cart2gnomonic(GTVector<GTPoint<T>> &clist, T rad, T xlatc, T xlongc, GTVector<GTPoint<T>> &glist)
{
	GEOFLOW_TRACE();
  GString serr = "GridIcos::cart2gnomonic: ";


  // Gnomonic transform given by (th=lat, and phi = long; thc, phic refer
  // to refernce lat and long):
  // x = rad [ cos(th) sin(phi-phic) ] / 
  //     [sin(thc)sin(th) + cos(thc)cos(th)cos(phi-phic)]
  // y = rad [ cos(thc)sin(th) - sin(thc)cos(th)cos(phi-phic) ] / 
  //     [sin(thc)sin(th) + cos(thc)cos(th)cos(th-thc)]
  // This can be rewritten as:
  // x = rad tan(longp), y = rad tan(latp) sec(longp)
  // where
  // latp  = arcsin[cos(thc)sin(th) - sin(thc)cos(th)cos(phi-phic)]
  // longp = atan [ [ cos(th) sin(phi-phic) / 
  //                  [sin(thc)sin(th) + cos(thc)cos(th)cos(phi-phic)] ] ]
  // 

  T xlatp, xlongp;
  T den, r, xlat, xlong;
  T eps = 10.0*std::numeric_limits<GFTYPE>::epsilon();
  for ( GSIZET i=0; i<clist.size(); i++ ) { // loop over all points
    r      = clist[i].norm();
    xlat   = asin(clist[i].x3/r);
    xlong  = atan2(clist[i].x2,clist[i].x1);
    xlong  = xlong < 0.0 ? 2.0*PI+xlong : xlong;
#if 0
    den    = sin(xlatc)*sin(xlat) + cos(xlatc)*cos(xlat)*cos(xlong-xlongc);  
    xlatp  = asin( cos(xlatc)*sin(xlat) - sin(xlatc)*cos(xlat)*cos(xlong-xlongc) );
    xlongp = atan2( cos(xlat)*sin(xlong-xlongc),den );

    glist[i].x1 = rad*tan(xlongp);
    glist[i].x2 = rad*tan(xlatp)*sec(xlongp);
#else
    den    = sin(xlatc)*sin(xlat) + cos(xlatc)*cos(xlat)*cos(xlong-xlongc); 
    den    = fabs(den) < eps ? 0.0 : 1.0/den;
    glist[i].x1 = rad*cos(xlat)*sin(xlong-xlongc)*den;
    glist[i].x2 = rad*( cos(xlatc)*sin(xlat) - sin(xlatc)*cos(xlat)*cos(xlong-xlongc) ) * den;
#endif
  }

} // end of method cart2gnomonic


//**********************************************************************************
//**********************************************************************************
// METHOD : gnomonic2cart (1)
// DESC   : Transform gnomonic coords to Cartesian space
// ARGS   : glist: Vector of gnomonic coords
//          rad  : radius of sphere
//          xlatc,
//          xlongc: lat and long of reference vector (usually a centroid)
//          clist: converted Cartesian coords
//          
// RETURNS: none.
//**********************************************************************************
template<typename T>
void GGridIcos::gnomonic2cart(GTVector<GTVector<T>> &glist, T rad, T xlatc, T xlongc, GTVector<GTVector<T>> &clist)
{
	GEOFLOW_TRACE();
  GString serr = "GridIcos::gnomonic2cart (1): ";
  assert(glist.size() >= 2 && clist.size() >= 3 && "Incompaible coordinate dimensions");


  // Gnomonic transform given by (th=lat, and phi = long; thc, phic refer
  // to refernce lat and long):
  //   x = rad [ cos(th) sin(phi-phic) ] / 
  //       [sin(thc)sin(th) + cos(thc)cos(th)cos(phi-phic)]
  //   y = rad [ cos(thc)sin(th) - sin(thc)cos(th)cos(phi-phic) ] / 
  //       [sin(thc)sin(th) + cos(thc)cos(th)cos(th-thc)]
  // 
  // Reverse this here... Solving for tan(th), and cos(phi-phic), arrive at:
  //   th  = asin( cos(b)sin(thc) + y*sin(b)*cos(thc)/rho )
  //   phi = phic + atan( x*sin(b) / ( rho*cos(thc)*cos(b) - y*sin(thc)*sin(b) ) ) 
  // where
  //   rho = sqrt(x^2 + y^2)
  //   b   = atan(rho)  
  //
  // (From Wolfram Research)

  T X, Y;
  T beta, rho, sign, x, xlat, xlong, y;
  T eps = 10.0*std::numeric_limits<GFTYPE>::epsilon();
  for ( GSIZET i=0; i<glist[0].size(); i++ ) { // loop over all points
    x      = glist[0][i];
    y      = glist[1][i];


    if ( fabs(x) < eps && fabs(y) < eps ) {
      xlat = xlatc;
      xlong = xlongc;
    } else {
      rho    = sqrt(x*x + y*y);
      beta   = atan(rho);
      Y      = cos(beta)*sin(xlatc) + (y*sin(beta)*cos(xlatc))/rho;
      sign   = copysign(1.0, Y);
      Y      = sign* MIN(fabs(Y),1.0);
      xlat   = asin(Y);
      Y      = x*sin(beta);
      X      = rho*cos(xlatc)*cos(beta)-y*sin(xlatc)*sin(beta);
      xlong  = xlongc + atan2(Y, X);
    }

    // Convert to spherical-polar to  Cart. coordinates:
    clist[0][i] = rad*cos(xlat)*cos(xlong);
    clist[1][i] = rad*cos(xlat)*sin(xlong);
    clist[2][i] = rad*sin(xlat);
  }

} // end of method gnomonic2cart (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : gnomonic2cart (2)
// DESC   : Transform gnomonic coords to Cartesian space
// ARGS   : glist: Vector of gnomonic coords
//          rad  : radius of sphere
//          xlatc,
//          xlongc: lat and long of reference vector (usually a centroid)
//          clist: converted Cartesian coords
//          
// RETURNS: none.
//**********************************************************************************
template<typename T>
void GGridIcos::gnomonic2cart(GTVector<GTPoint<T>> &glist, T rad, T xlatc, T xlongc, GTVector<GTPoint<T>> &clist)
{
	GEOFLOW_TRACE();
  GString serr = "GridIcos::gnomonic2cart(2): ";


  // Gnomonic transform given by (th=lat, and phi = long; thc, phic refer
  // to refernce lat and long):
  //   x = rad [ cos(th) sin(phi-phic) ] / 
  //       [sin(thc)sin(th) + cos(thc)cos(th)cos(phi-phic)]
  //   y = rad [ cos(thc)sin(th) - sin(thc)cos(th)cos(phi-phic) ] / 
  //       [sin(thc)sin(th) + cos(thc)cos(th)cos(th-thc)]
  // Reverse this here... Solving for tan(th), and cos(phi-phic), arrive at:
  //   th  = asin( cos(b)sin(thc) + y*sin(b)*cos(thc)/rho )
  //   phi = phic + atan( x*sin(b) / ( rho*cos(thc)*cos(b) - y*sin(thc)*sin(b) ) ) 
  // where
  //   rho = sqrt(x^2 + y^2)
  //   b   = atan(rho)  
  // (From Wolfram Research)

  T X, Y;
  T beta, rho, sign, x, xlat, xlong, y;
  T eps = 10.0*std::numeric_limits<GFTYPE>::epsilon();
  for ( GSIZET i=0; i<clist.size(); i++ ) { // loop over all points
    x      = glist[i][0];
    y      = glist[i][1];

    if ( fabs(x) < eps && fabs(y) < eps ) {
      xlat = xlatc;
      xlong = xlongc;
    } else {
      rho    = sqrt(x*x + y*y);
      beta   = atan(rho);
      Y      = cos(beta)*sin(xlatc) + (y*sin(beta)*cos(xlatc))/rho;
      sign   = copysign(1.0, Y);
      Y      = sign* MIN(fabs(Y),1.0);
      xlat   = asin(Y);
      Y      = x*sin(beta);
      X      = rho*cos(xlatc)*cos(beta)-y*sin(xlatc)*sin(beta);
      xlong  = xlongc + atan2(Y, X);
    }
    // Convert to spherical-polar to  Cart. coordinates:
    clist[i][0] = rad*cos(xlat)*cos(xlong);
    clist[i][1] = rad*cos(xlat)*sin(xlong);
    clist[i][2] = rad*sin(xlat);
  }

} // end of method gnomonic2cart(2)


//**********************************************************************************
//**********************************************************************************
// METHOD : reorderverts2d
// DESC   : Reorder specified elem vertices to be consistent with
//          shape functions
// ARGS   : uverts : list of unordered 3-vertices
//          tverts : temp array of points equal in number to uverts
//          isort  : array of sorting indirection indices
//          overts : array of ordered 3-vertices, returned
// RETURNS: none.
//**********************************************************************************
template<typename T>
void GGridIcos::reorderverts2d(GTVector<GTPoint<T>> &uverts, GTVector<GTPoint<T>> &tverts, GTVector<GSIZET> &isort, 
                               GTVector<GTPoint<T>> &overts)
{
	GEOFLOW_TRACE();
  GString serr = "GridIcos::reorderverts2d: ";

  assert(uverts.size() == 4 
      && overts.size() == 4
      && "Incorrect number of vertices");

  T                xlatc, xlongc;
  GTPoint<T>       c(3);
  GTVector<T>      x(4);
  GTVector<GSIZET> Ixy(4);
  
  for ( auto j=0; j<tverts.size(); j++ ) tverts[j].resize(2);

  c = (uverts[0] + uverts[1] + uverts[2] + uverts[3])*0.25; // elem centroid
  xlatc  = asin(c.x3); // reference lat/long
  xlongc = atan2(c.x2,c.x1);
  xlongc = xlongc < 0.0 ? 2*PI+xlongc : xlongc;

  // Convert to gnomonic coords so that we
  // can examine truly 2d planar coords when
  // re-ordering:
  cart2gnomonic<T>(uverts, 1.0, xlatc, xlongc, tverts); // gnomonic vertices of quads

  isort.resize(4);
  for ( auto i=0; i<tverts.size(); i++ ) { 
    x[i] = tverts[i].x1;
  }

  x.sortincreasing(Ixy);

  // Do 'left' -hand vertices:
  if ( tverts[Ixy[0]].x2 < tverts[Ixy[1]].x2 ) {
    isort[0] = Ixy[0];
    isort[3] = Ixy[1];
  } else {
    isort[0] = Ixy[1];
    isort[3] = Ixy[0];
  }

  // Do 'right' -hand vertices:
  if ( tverts[Ixy[2]].x2 < tverts[Ixy[3]].x2 ) {
    isort[1] = Ixy[2];
    isort[2] = Ixy[3];
  } else {
    isort[1] = Ixy[3];
    isort[2] = Ixy[2];
  }

  for ( auto j=0; j<4; j++ ) overts[j] = uverts[isort[j]];
  
} // end of method reorderverts2d



//**********************************************************************************
//**********************************************************************************
// METHOD : copycast (1)
// DESC   : copy 'from' array to 'to' array, while casting. Also, make sure
//          sizes are correct.
// ARGS   : from : 'from' array
//          to   : 'to' array; may be of different type than 'from'
// RETURNS: none.
//**********************************************************************************
template<typename TF, typename TT>
void GGridIcos::copycast(GTVector<GTVector<TF>> &from, GTVector<GTVector<TT>> &to)
{
	GEOFLOW_TRACE();
  assert(to.size() == from.size() && "Incompatible dimensions");
  for ( auto j=0; j<to.size(); j++ ) {
    to[j].resize(from[j].size()); 
    for ( auto i=0; i<to[j].size(); i++ ) {
      to[j][i] = static_cast<TT>(from[j][i]);
    }
  }
  
} // end of method copycast (1)



//**********************************************************************************
//**********************************************************************************
// METHOD : copycast (2)
// DESC   : copy 'from' array to 'to' array, while casting. Also, make sure
//          sizes are correct.
// ARGS   : from : 'from' array
//          to   : 'to' array; may be of different type than 'from'
// RETURNS: none.
//**********************************************************************************
template<typename TF, typename TT>
void GGridIcos::copycast(GTVector<GTVector<TF>*> &from, GTVector<GTVector<TT>*> &to)
{
	GEOFLOW_TRACE();
  assert(to.size() == from.size() && "Incompatible dimensions");
  for ( auto j=0; j<to.size(); j++ ) {
    to[j]->resize(from[j]->size()); 
    for ( auto i=0; i<to[j]->size(); i++ ) {
      (*to[j])[i] = static_cast<TT>((*from[j])[i]);
    }
  }
  
} // end of method copycast (2)



//**********************************************************************************
//**********************************************************************************
// METHOD : copycast (3)
// DESC   : copy 'form' point to 'to' point, while casting
// ARGS   : from : 'from' point
//          to   : 'to' point; may be of different type than 'from'
// RETURNS: none.
//**********************************************************************************
template<typename TF, typename TT>
void GGridIcos::copycast(GTPoint<TF> &from, GTPoint<TT> &to)
{
	GEOFLOW_TRACE();
  for ( auto j=0; j<to.size(); j++ ) {
    to[j] = static_cast<TT>(from[j]);
  }
  
} // end of method copycast (3)



//**********************************************************************************
//**********************************************************************************
// METHOD : interleave
// DESC   : Utility routine to interleave 2 rows of points 
// ARGS   : R0   : row of points above
//          R1   : row of points below
//          I    : 'row index' (0 ... iLevel+1 of R1 
//          Rz   : list of vertices (points) interleaved, so that following
//                 ordering represents triangles in 'Lagrangian refinement':
//                 (Rz(0), Rz(1), Rz(2)); (Rz(1) Rz(2) Rz(3)) , ...
// RETURNS: none.
//**********************************************************************************
template<typename T>
void GGridIcos::interleave(GTVector<GTPoint<T>> &R0, GTVector<GTPoint<T>> &R1,
                   GINT I, GTVector<GTPoint<T>> &Rz)
{
	GEOFLOW_TRACE();
  GString serr = "GridIcos::interleave: ";

  // Interlaeave R0 and R1:
  for ( GSIZET j=0; j<Rz.size(); j++ ) {
    if ( j%2 == 0 ) Rz   [j] = R1[j/2];
    else            Rz   [j] = R0[(j-1)/2];
  }

} // end of method interleave


//**********************************************************************************
//**********************************************************************************
// METHOD : lagvert
// DESC   : Utility routine to compute 'Lagrangian-refined' vertices from 'base' 
//          vertices and vertex indices. Not vectorized.
//          Given bse vertices, find vertices at 'row index' I
// ARGS   : a,b,c: base vertices
//          I    : 'row index' (0 ... iLevel+1
//          R    : list of vertices (points) at I
// RETURNS: none.
//**********************************************************************************
template<typename T>
void GGridIcos::lagvert(GTPoint<T>&a, GTPoint<T> &b, GTPoint<T> &c,
                   GINT I, GTVector<GTPoint<T>> &R)
{
	GEOFLOW_TRACE();
  GString serr = "GridIcos::lagvert: ";

  T   xI, xJ;

  GTPoint<T> rL(3);
  GTPoint<T> rR(3);

  R.resizem(I+1);

  T fact = 1.0/static_cast<T>(nrows_+1);

  // Build 'rail' points on L and R:
  xI = static_cast<T>(I);
  rL = a + ( (b - a) * (xI * fact) );
  rR = a + ( (c - a) * (xI * fact) );

  // Compute R vertices based on refinement indices:
  fact = I > 0 ? 1.0/static_cast<T>(I) : 1.0;
  for ( GSIZET j=0; j<I+1; j++ ) {
    xJ = static_cast<T>(j);
    R[j] = rL + (rR - rL)*(xJ*fact);
  }

} // end of method lagvert


//**********************************************************************************
//**********************************************************************************
// METHOD : order_latlong2d
// DESC   : Order 2d vertices on exit s.t. they roughly define a 'box' 
//          in spherical coords. Used only for quad elements.
// ARGS   : verts : Array of vertices, re-ordered on exit
// RETURNS: none.
//**********************************************************************************
template<typename T>
void GGridIcos::order_latlong2d(GTVector<GTPoint<T>> &verts)
{
	GEOFLOW_TRACE();
  assert(verts.size() == 4 && "4 vertices must be provided");

  GString              serr = "GGridIcos::order_latlong2d: ";
  GTVector<GSIZET>     isortlon(4);
  GTVector<T>  lon(4), lat(4);
  GTVector<GTPoint<T>> cverts(4);       // copy of input verts
  GTVector<GTPoint<T>> sverts(4);       // sverts in sph coords

  cverts = verts;
  sverts = verts;
  xyz2spherical<T>(sverts); // convert verts to latlon

  // Isolate lat, lon:
  for ( GSIZET j=0; j<4; j++ ) {
    lat[j] = sverts[j].x2;
    lon[j] = sverts[j].x3;
  }

  // Sort lon in increasing order:
  lon.sortincreasing(isortlon);

  // Check vertices near 0-2pi axis:
  if ( fabs(lon[isortlon[0]] - lon[isortlon[3]]) < PI ) {
    for ( GSIZET j=0; j<4; j++ ) {
      if ( lon[j] > 1.5*PI && lon[j] <=2.0*PI ) lon[j] -= 2.0*PI;
    }
  }

  lon.sortincreasing(isortlon);

  // Find 2 points with smallest lon, set v0 and v3
  // based on lat to define 'box':
  if ( lat[isortlon[0]] < lat[isortlon[1]] ) {
    verts[0] = cverts[isortlon[0]];
    verts[3] = cverts[isortlon[1]];
  }
  else {
    verts[0] = cverts[isortlon[1]];
    verts[3] = cverts[isortlon[0]];
  }
  
  // Find 2 points with largest lon, set v1 and v2
  // based on lat to define 'box':
  if ( lat[isortlon[2]] < lat[isortlon[3]] ) {
    verts[1] = cverts[isortlon[2]];
    verts[2] = cverts[isortlon[3]];
  }
  else {
    verts[1] = cverts[isortlon[3]];
    verts[2] = cverts[isortlon[2]];
  }

#if 0
  cout << serr << " on entry: verts=" << cverts << endl;
  cout << serr << " on exit : verts=" << verts << endl;
#endif

} // end, method order_latlong2d


//**********************************************************************************
//**********************************************************************************
// METHOD : order_triangles
// DESC   : Order triangle vertices
// ARGS   : tmesh: Array of vertices, re-ordered on exit
// RETURNS: none.
//**********************************************************************************
template<typename T>
void GGridIcos::order_triangles(GTVector<GTriangle<T>> &tmesh)
{
	GEOFLOW_TRACE();
  GString              serr = "GGridIcos::order_triangles: ";
  GTVector<GSIZET>     isortlon(3);
  GTVector<T>          lon(3), lat(3);
  GTVector<GTPoint<T>> cverts(3);       // copy of input verts
  GTVector<GTPoint<T>> sverts(3);       // sverts in sph coords

  iup_.resize(tmesh.size());
  iup_ = 0;
  for ( GSIZET i=0; i<tmesh.size(); i++ ) {
    for ( GSIZET j=0; j<3; j++ ) cverts[j] = *tmesh[i].v[j];
    sverts = cverts;
    xyz2spherical<T>(sverts); // convert verts to latlon
    for ( GSIZET j=0; j<3; j++ ) {
      lat[j] = sverts[j].x2;
      lon[j] = sverts[j].x3;
    }
    lon.sortincreasing(isortlon);

    // Check vertices near 0-2pi axis; if triangle
    // spans it, subtract 2pi to make longitude negative:
    if ( fabs(lon[isortlon[0]] - lon[isortlon[2]]) < PI ) {
      for ( GSIZET j=0; j<3; j++ ) {
        if ( lon[j] > 1.5*PI && lon[j] <= 2.0*PI ) lon[j] -= 2.0*PI;
      }
    }
    lon.sortincreasing(isortlon);
    
    if ( lat[isortlon[1]] > lat[isortlon[0]]
      && lat[isortlon[1]] > lat[isortlon[2]] ) { // pointing upwards
     *tmesh[i].v[0] =  cverts[isortlon[0]];
     *tmesh[i].v[1] =  cverts[isortlon[2]];
     *tmesh[i].v[2] =  cverts[isortlon[1]];
      iup_[i] = 1;
    }

    if ( lat[isortlon[1]] < lat[isortlon[0]]
      && lat[isortlon[1]] < lat[isortlon[2]] ) { // pointing downwards
     *tmesh[i].v[0] =  cverts[isortlon[1]];
     *tmesh[i].v[1] =  cverts[isortlon[2]];
     *tmesh[i].v[2] =  cverts[isortlon[0]];
    }
  }

} // end, method order_triangles


