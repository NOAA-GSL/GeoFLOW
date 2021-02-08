//==================================================================================
// Module       : gmtk.cpp
// Date         : 1/31/18 (DLR)
// Description  : Math TooKit: namespace of C routines for various
//                math functions
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

#include "ggrid.hpp"
#include "ggrid_box.hpp"
#include "ggrid_icos.hpp"
#include "gmtk.hpp"
#include "tbox/tracer.hpp"


using namespace std;

//#if !defined(_GMTK_GLOBAL_DATA)
//  #define _GMTK_GLOBAL_DATA
  GINT szMatCache_ = _G_MAT_CACHE_SIZE;
  GINT szVecCache_ = _G_VEC_CACHE_SIZE;
//#endif


namespace GMTK 
{

//**********************************************************************************
//**********************************************************************************
// METHOD : curl 
// DESC   : Compute curl component, idir, of input vector field
//          
// ARGS   : grid : grid
//          u    : input vector field. Must have >= GDIM components.
//          idir : curl component to compute. Must be appropriate for 
//                 problem dimension.
//          tmp  : tmp vector; must be of at least length 2.
//          curlc: result
// RETURNS: none.
//**********************************************************************************
template<>
void curl(GGrid &grid, const GTVector<GTVector<GFTYPE>*> &u, const GINT idir, 
          GTVector<GTVector<GFTYPE>*> &tmp, GTVector<GFTYPE> &curlc)
{
  GEOFLOW_TRACE();

  assert(tmp.size() >= 2 && "Insufficient temp space");

  // Handle 1c cases in 2d or 3d:
  if  ( u.size() < 2 ) {
     curlc = 0.0; 
  }

  // Handle 2-d case:
  else if ( GDIM == 2 && u.size() >= GDIM && grid.gtype() != GE_2DEMBEDDED ) {
    switch (idir) {
      case 1:
      case 2:
        curlc = 0.0;
        break;
      case 3:
        grid.deriv(*u[1], 1, *tmp[0], curlc);
        grid.deriv(*u[0], 2, *tmp[0], *tmp[1]);
        curlc -= *tmp[1];
        break;
      default:
        assert( FALSE && "Invalid component specified");
        break;
    }
  }

  // Handle 2.5-d or 2d-3c case:
  else if ( GDIM == 2 && u.size() > GDIM && grid.gtype() != GE_2DEMBEDDED ) {
    switch (idir) {
      case 1:
        grid.deriv(*u[2], 2, *tmp[0], curlc);
        curlc *= -1.0;
        break;
      case 2:
        grid.deriv(*u[2], 1, *tmp[0], curlc);
        break;
      case 3:
        grid.deriv(*u[1], 1, *tmp[0], curlc);
        grid.deriv(*u[0], 2, *tmp[0], *tmp[1]);
        curlc -= *tmp[1];
        break;
      default:
        assert( FALSE && "Invalid component specified");
        break;
    }
  }

  // Handle 2d-2c regular types:
  else if ( GDIM == 2  && u.size() == 2 && grid.gtype() == GE_REGULAR ) {
    switch (idir) {
      case 1:
      case 2:
        curlc = 0.0;
        break;
      case 3:
        grid.deriv(*u[1], 1, *tmp[0], curlc);
        grid.deriv(*u[0], 2, *tmp[0], *tmp[1]);
        curlc -= *tmp[1];
        break;
      default:
        assert( FALSE && "Invalid component specified");
        break;
    }
  }

  // Handle 3d-3c or embedded cases:
  else if ( GDIM == 3 || grid.gtype() == GE_2DEMBEDDED ) {
    switch (idir) {
      case 1:
        grid.deriv(*u[1], 3, *tmp[0], curlc);
        grid.deriv(*u[2], 2, *tmp[0], *tmp[1]);
        curlc -= *tmp[1];
        break;
      case 2:
        grid.deriv(*u[2], 1, *tmp[0], curlc);
        grid.deriv(*u[0], 3, *tmp[0], *tmp[1]);
        curlc -= *tmp[1];
        break;
      case 3:
        grid.deriv(*u[1], 1, *tmp[0], curlc);
        grid.deriv(*u[0], 2, *tmp[0], *tmp[1]);
        curlc -= *tmp[1];
        break;
    }
  }
  else {
    assert(FALSE && "Curl cannot be computed");
  }


  return;

} // end of method curl


//**********************************************************************************
//**********************************************************************************
// METHOD : grad
// DESC   : Compute gradient component, idir, of input vector field.
//          Note: Don't really need this, as it's just another way
//                to refer to the Cartesian 'deriv' method in GGrid
//          
// ARGS   : grid : grid
//          u    : input (scalar) field. 
//          idir : gradient component to compute. Must be appropriate for 
//                 problem dimension.
//          tmp  : tmp vector; must be of at least length 1.
//          gradc: result
// RETURNS: none.
//**********************************************************************************
template<>
void grad(GGrid &grid, GTVector<GFTYPE> &u, const GINT idir, 
          GTVector<GTVector<GFTYPE>*> &tmp, GTVector<GFTYPE> &gradc)
{
  GEOFLOW_TRACE();
  assert ( idir >0 && idir <=3 && "Invalid compoment specified");

  grid.deriv(u, idir, *tmp[0], gradc);

} // end of method grad


//**********************************************************************************
//**********************************************************************************
// METHOD : constrain2sphere(1)
// DESC   : Project/constrain input 3-vector to sphere:
//          -  -     -                           -   -  -
//          |vx|     |(r^2-x^2)   -xy     -xz    |   |vx|
//        P |vy| =   |   -yx    (r^2-y^2) -yz    |   |vy|
//          |vz|     |   -zx      -zy   (r^2-z^2)|   |vz|
//          -  -     -                           -   -  -
//
//        This is derived from a Lagrange multiplier constraint
//        that requires all vectors, v, to be normal to radial
//        vector, s.t. x.v = 0. 
//          
// ARGS   : grid : Grid. If not of the correct type, nothing is done
//          v    : Array of vector components to be constrained
//          Pv   : Projected vector. Must be different from v
// RETURNS: none
//**********************************************************************************
template<>
void constrain2sphere(GGrid &grid, const GTVector<GTVector<GFTYPE>*> &v, GTVector<GTVector<GFTYPE>*> &Pv)
{
  GEOFLOW_TRACE();
  if ( grid.gtype() != GE_2DEMBEDDED || v.size() != 3 ) return;

  assert( v.size() >= 3 && "Incompatible dimensionality");

  GSIZET nxy = grid.ndof();
  GFTYPE r2, x, y, z;
  GTVector<GTVector<GFTYPE>> *xnodes = &grid.xNodes();

  for ( GSIZET j=0; j<nxy; j++ ) {
    x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
    r2 = x*x + y*y + z*z;
    (*Pv[0])[j] =  (*v[0])[j]*(r2-x*x) - (*v[1])[j]*x*y      - (*v[2])[j]*x*z;
    (*Pv[1])[j] = -(*v[0])[j]*y*x      + (*v[1])[j]*(r2-y*y) - (*v[2])[j]*y*z;
    (*Pv[2])[j] = -(*v[0])[j]*z*x      - (*v[1])[j]*z*y      + (*v[2])[j]*(r2-z*z);
   }

} // end of method constrain2sphere (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : constrain2sphere(2)
// DESC   : Project/constrain input 3-vector to sphere:
//          -  -     -                           -   -  -
//          |vx|     |(r^2-x^2)   -xy     -xz    |   |vx|
//        P |vy| =   |   -yx    (r^2-x^2) -yz    |   |vy|
//          |vz|     |   -zx      -zy   (r^2-z^2)|   |vz|
//          -  -     -                           -   -  -
//
//        This is derived from a Lagrange multiplier constraint
//        that requires all vectors, v, to be normal to radial
//        vector, s.t. x.v = 0. 
//          
// ARGS   : grid : Grid. If not of the correct type, nothing is done
//          v    : Array of vector components to be constrained, 
//                 modified on output to be projected components
// RETURNS: none
//**********************************************************************************
template<>
void constrain2sphere(GGrid &grid, GTVector<GTVector<GFTYPE>*> &v)
{
  GEOFLOW_TRACE();
  if ( grid.gtype() != GE_2DEMBEDDED || v.size() < 3 ) return;

  GSIZET nxy = grid.ndof();
  GFTYPE ri, r2, x, y, z;
  GTVector<GTVector<GFTYPE>> *xnodes = &grid.xNodes();

  x = (*xnodes)[0][0]; y = (*xnodes)[1][0]; z = (*xnodes)[2][0];
  r2 = x*x + y*y + z*z;
  ri = 1.0/r2;
  for ( GSIZET j=0; j<nxy; j++ ) {
    x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
    (*v[0])[j] = ( (*v[0])[j]*(r2-x*x) - (*v[1])[j]*x*y      - (*v[2])[j]*x*z )*ri;
    (*v[1])[j] = (-(*v[0])[j]*y*x      + (*v[1])[j]*(r2-y*y) - (*v[2])[j]*y*z )*ri;
    (*v[2])[j] = (-(*v[0])[j]*z*x      - (*v[1])[j]*z*y      + (*v[2])[j]*(r2-z*z) )*ri;
   }

} // end of method constrain2sphere (2)



//**********************************************************************************
//**********************************************************************************
// METHOD : vsphere2cart 
// DESC   : Convert vector field from spherical coords to Cartesian.
//
// ARGS   : grid : Grid. If not of the correct type, nothing is done
//          vsph : Array of vector components. If we have GE_2DEMBEDDED grid,
//                 there must be at least 2 components, and only the first 2
//                 are used, and assumed to be latitudual, and longitudinal
//                 respectively. If grid is a 3D spherical grid, then
//                 vector components are assumed to be (r, lat, long).
//                 If grid is REGULAR, then this transformation cannot be done.
//          vtype: Vector type of spherical coords GVECTYPE_(PHYS, CONTRAVAR, COVAR)
//          vcart: Cartesian vector field. Must have at least 3 components, and
//                 are returned as (x, y, z).
// RETURNS: none
//**********************************************************************************
template<>
void vsphere2cart(GGrid &grid, const GTVector<GTVector<GFTYPE>*> &vsph, GVectorType vtype, GTVector<GTVector<GFTYPE>*> &vcart)
{
  GEOFLOW_TRACE();
if      ( GDIM == 2 && grid.gtype() == GE_2DEMBEDDED ) {
assert( vsph.size() >= 2 && "GE_2DEMBEDDED grid requires 2 spherical components");
}
else if ( grid.gtype() == GE_DEFORMED ) {
assert( vsph.size() >= 3 && "GE_DEFORMED grid requires 3 spherical components");
}
else if ( grid.gtype() != GE_REGULAR ) {
assert( FALSE && "GE_REGULAR grid will not allow this transformation");
}


GSIZET           nxy = grid.ndof();
GFTYPE           x, y, z, tiny;
GFTYPE           phi, r, theta;
GFTYPE           gpp, gtt;
GFTYPE           vthcontra, vphicontra;
GTVector<GTVector<GFTYPE>> *xnodes = &grid.xNodes();

tiny = numeric_limits<GFTYPE>::epsilon();

//   v_i_cart = vtheta dx_i/dtheta + vphi dx_i/dphi
// where
//   vtheta, vhi are _contravariant_ (upper-index) components
// Note: Metric is orthogonal:
//   g_ij = (1, h_theta^2, h_phi^2) = (1, r^2, (r cos(theta))^2 )
if ( grid.gtype() == GE_2DEMBEDDED ) {
if ( vtype == GVECTYPE_PHYS ) { // vsph are physical components
for ( GSIZET j=0; j<nxy; j++ ) {
x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
r     = sqrt(x*x + y*y + z*z);
theta = asin(z/r);
phi   = atan2(y,x);
(*vcart[0])[j] = -(*vsph[0])[j]*  sin(theta)*cos(phi) 
	       -  (*vsph[1])[j]*             sin(phi);
(*vcart[1])[j] = -(*vsph[0])[j]*  sin(theta)*sin(phi) 
	       +  (*vsph[1])[j]*             cos(phi);
(*vcart[2])[j] =  (*vsph[0])[j]*  cos(theta);
}
}
else if ( vtype == GVECTYPE_COVAR ) { // vsph are covar. components
for ( GSIZET j=0; j<nxy; j++ ) {
x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
r     = sqrt(x*x + y*y + z*z);
theta = asin(z/r);
phi   = atan2(y,x);
vthcontra  = (*vsph[0])[j];
vphicontra = (*vsph[1])[j];
gtt        = r*r; gpp = pow(r*cos(theta),2);
vthcontra  = (*vsph[0])[j]/(gtt+tiny);
vphicontra = (*vsph[1])[j]/(gpp+tiny);
(*vcart[0])[j] = -vthcontra *r*sin(theta)*cos(phi) 
	       -  vphicontra*r*cos(theta)*sin(phi);
(*vcart[1])[j] = -vthcontra *r*sin(theta)*sin(phi) 
	       +  vphicontra*r*cos(theta)*cos(phi);
(*vcart[2])[j] =  vthcontra *r*cos(theta);
}
}
else if ( vtype == GVECTYPE_CONTRAVAR ) { // vsph are contravar. components
for ( GSIZET j=0; j<nxy; j++ ) {
x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
r     = sqrt(x*x + y*y + z*z);
theta = asin(z/r);
phi   = atan2(y,x);
vthcontra  = (*vsph[0])[j];
vphicontra = (*vsph[1])[j];
(*vcart[0])[j] = -vthcontra *r*sin(theta)*cos(phi) 
	       -  vphicontra*r*cos(theta)*sin(phi);
(*vcart[1])[j] = -vthcontra *r*sin(theta)*sin(phi) 
	       +  vphicontra*r*cos(theta)*cos(phi);
(*vcart[2])[j] =  vthcontra *r*cos(theta);
}
}
}

if ( grid.gtype() == GE_DEFORMED ) {
if      ( vtype == GVECTYPE_PHYS ) {
for ( GSIZET j=0; j<nxy; j++ ) {
x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
r     = sqrt(x*x + y*y + z*z);
theta = asin(z/r);
phi   = atan2(y,x);
gtt        = r; gpp = r*cos(theta);
vthcontra  = (*vsph[1])[j]/(gtt+tiny);
vphicontra = (*vsph[2])[j]/(gpp+tiny);
(*vcart[0])[j] =  (*vsph[0])[j]*  cos(theta)*cos(phi)
	       -  (*vsph[1])[j]*r*sin(theta)*cos(phi) 
	       +  (*vsph[2])[j]*r*cos(theta)*sin(phi);
(*vcart[1])[j] =  (*vsph[0])[j]*  cos(theta)*sin(phi)
	       -  (*vsph[1])[j]*r*sin(theta)*sin(phi) 
	       +  (*vsph[2])[j]*r*cos(theta)*cos(phi);
(*vcart[2])[j] =  (*vsph[0])[j]*  sin(theta)
	       +  (*vsph[1])[j]*r*cos(theta);
}
}
else if ( vtype == GVECTYPE_COVAR ) {
for ( GSIZET j=0; j<nxy; j++ ) {
x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
r     = sqrt(x*x + y*y + z*z);
theta = asin(z/r);
phi   = atan2(y,x);
gtt        = r*r; gpp = pow(r*cos(theta),2);
vthcontra  = (*vsph[1])[j]/(gtt+tiny);
vphicontra = (*vsph[2])[j]/(gpp+tiny);
(*vcart[0])[j] =  (*vsph[0])[j]*  cos(theta)*cos(phi)
	       -  (*vsph[1])[j]*r*sin(theta)*cos(phi) 
	       +  (*vsph[2])[j]*r*cos(theta)*sin(phi);
(*vcart[1])[j] =  (*vsph[0])[j]*  cos(theta)*sin(phi)
	       -  (*vsph[1])[j]*r*sin(theta)*sin(phi) 
	       +  (*vsph[2])[j]*r*cos(theta)*cos(phi);
(*vcart[2])[j] =  (*vsph[0])[j]*  sin(theta)
	       +  (*vsph[1])[j]*r*cos(theta);
}
}
else if ( vtype == GVECTYPE_CONTRAVAR ) {
for ( GSIZET j=0; j<nxy; j++ ) {
x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
r     = sqrt(x*x + y*y + z*z);
theta = asin(z/r);
phi   = atan2(y,x);
vthcontra  = (*vsph[0])[j];
vphicontra = (*vsph[1])[j];
(*vcart[0])[j] =  (*vsph[0])[j]*  cos(theta)*cos(phi)
	       -  (*vsph[1])[j]*r*sin(theta)*cos(phi) 
	       +  (*vsph[2])[j]*r*cos(theta)*sin(phi);
(*vcart[1])[j] =  (*vsph[0])[j]*  cos(theta)*sin(phi)
	       -  (*vsph[1])[j]*r*sin(theta)*sin(phi) 
	       +  (*vsph[2])[j]*r*cos(theta)*cos(phi);
(*vcart[2])[j] =  (*vsph[0])[j]*  sin(theta)
	       +  (*vsph[1])[j]*r*cos(theta);
}
}
}

} // end of method vsphere2cart


//**********************************************************************************
//**********************************************************************************
// METHOD : vcart2sphere
// DESC   : Convert vector from Cartesian coords to spherical coords.
//
// ARGS   : grid : Grid object
//          vcart: Cartesian vector field. Must have at least 3 components, and
//                 are returned as (x, y, z).
//          vtype: Vector type in spherical coords GVECTYPE_(PHYS, CONTRAVAR, COVAR)
//          vsph : Array of vector components. If we have GE_2DEMBEDDED grid,
//                 there must be at least 2 components, and only the first 2
//                 are used, and assumed to be latitudinal, and longitudinal
//                 respectively. If grid is a 3D spherical grid, then
//                 vector components are assumed to be (r, lat, long).
//                 If grid is REGULAR, then this transformation cannot be done.
// RETURNS: none
//**********************************************************************************
template<>
void vcart2sphere(GGrid &grid, const GTVector<GTVector<GFTYPE>*> &vcart, GVectorType vtype, GTVector<GTVector<GFTYPE>*> &vsph)
{
  GEOFLOW_TRACE();
assert( vcart.size() >= 3 && "Transformation requires 3 Cartesian components");
if      ( GDIM == 2 && grid.gtype() == GE_2DEMBEDDED ) {
assert( vsph.size() >= 2 && "GE_2DEMBEDDED grid requires 2 spherical components");
}
else if ( grid.gtype() == GE_DEFORMED ) {
assert( vsph.size() >= 3 && "GE_DEFORMED grid requires 3 spherical components");
}
else if ( grid.gtype() != GE_REGULAR ) {
assert( FALSE && "GE_REGULAR grid will not allow this transformation");
}

/*

GSIZET           nxy = grid.ndof();
GFTYPE           x, y, z, tiny;
GFTYPE           phi, r, theta;
GFTYPE           gpp, gtt;
GFTYPE           vthcontra, vphicontra;
GTVector<GTVector<GFTYPE>> *xnodes = &grid.xNodes();

tiny = numeric_limits<GFTYPE>::epsilon();

//   v_i_cart = vtheta dx_i/dtheta + vphi dx_i/dphi
// where
//   vtheta, vhi are _contravariant_ (upper-index) components
// Note: Metric is orthogonal:
//   g_ij = (1, h_theta^2, h_phi^2) = (1, r^2, (r cos(theta))^2 )
if ( grid.gtype() == GE_2DEMBEDDED ) {
if ( vtype == GVECTYPE_PHYS ) { // vsph are physical components
for ( GSIZET j=0; j<nxy; j++ ) {
x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
r     = sqrt(x*x + y*y + z*z);
theta = asin(z/r);
phi   = atan2(y,x);
(*vcart[0])[j] = -(*vsph[0])[j]*  sin(theta)*cos(phi) 
	       -  (*vsph[1])[j]*             sin(phi);
(*vcart[1])[j] = -(*vsph[0])[j]*  sin(theta)*sin(phi) 
	       +  (*vsph[1])[j]*             cos(phi);
(*vcart[2])[j] =  (*vsph[0])[j]*  sin(theta);
}
}
else if ( vtype == GVECTYPE_COVAR ) { // vsph are covar. components
for ( GSIZET j=0; j<nxy; j++ ) {
x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
r     = sqrt(x*x + y*y + z*z);
theta = asin(z/r);
phi   = atan2(y,x);
vthcontra  = (*vsph[0])[j];
vphicontra = (*vsph[1])[j];
gtt        = r*r; gpp = pow(r*cos(theta),2);
vthcontra  = (*vsph[0])[j]/(gtt+tiny);
vphicontra = (*vsph[1])[j]/(gpp+tiny);
(*vcart[0])[j] = -vthcontra *r*sin(theta)*cos(phi) 
	       -  vphicontra*r*cos(theta)*sin(phi);
(*vcart[1])[j] = -vthcontra *r*sin(theta)*sin(phi) 
	       +  vphicontra*r*cos(theta)*cos(phi);
(*vcart[2])[j] =  vthcontra *r*cos(theta);
}
}
else if ( vtype == GVECTYPE_CONTRAVAR ) { // vsph are contravar. components
for ( GSIZET j=0; j<nxy; j++ ) {
x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
r     = sqrt(x*x + y*y + z*z);
theta = asin(z/r);
phi   = atan2(y,x);
vthcontra  = (*vsph[0])[j];
vphicontra = (*vsph[1])[j];
(*vcart[0])[j] = -vthcontra *r*sin(theta)*cos(phi) 
	       -  vphicontra*r*cos(theta)*sin(phi);
(*vcart[1])[j] = -vthcontra *r*sin(theta)*sin(phi) 
	       +  vphicontra*r*cos(theta)*cos(phi);
(*vcart[2])[j] =  vthcontra *r*cos(theta);
}
}
}

if ( grid.gtype() == GE_DEFORMED ) {
if      ( vtype == GVECTYPE_PHYS ) {
for ( GSIZET j=0; j<nxy; j++ ) {
x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
r     = sqrt(x*x + y*y + z*z);
theta = asin(z/r);
phi   = atan2(y,x);
gtt        = r; gpp = r*cos(theta);
vthcontra  = (*vsph[1])[j]/(gtt+tiny);
vphicontra = (*vsph[2])[j]/(gpp+tiny);
(*vcart[0])[j] =  (*vsph[0])[j]*  cos(theta)*cos(phi)
	       -  (*vsph[1])[j]*r*sin(theta)*cos(phi) 
	       +  (*vsph[2])[j]*r*cos(theta)*sin(phi);
(*vcart[1])[j] =  (*vsph[0])[j]*  cos(theta)*sin(phi)
	       -  (*vsph[1])[j]*r*sin(theta)*sin(phi) 
	       +  (*vsph[2])[j]*r*cos(theta)*cos(phi);
(*vcart[2])[j] =  (*vsph[0])[j]*  sin(theta)
	       +  (*vsph[1])[j]*r*cos(theta);
}
}
else if ( vtype == GVECTYPE_COVAR ) {
for ( GSIZET j=0; j<nxy; j++ ) {
x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
r     = sqrt(x*x + y*y + z*z);
theta = asin(z/r);
phi   = atan2(y,x);
gtt        = r*r; gpp = pow(r*cos(theta),2);
vthcontra  = (*vsph[1])[j]/(gtt+tiny);
vphicontra = (*vsph[2])[j]/(gpp+tiny);
(*vcart[0])[j] =  (*vsph[0])[j]*  cos(theta)*cos(phi)
	       -  (*vsph[1])[j]*r*sin(theta)*cos(phi) 
	       +  (*vsph[2])[j]*r*cos(theta)*sin(phi);
(*vcart[1])[j] =  (*vsph[0])[j]*  cos(theta)*sin(phi)
	       -  (*vsph[1])[j]*r*sin(theta)*sin(phi) 
	       +  (*vsph[2])[j]*r*cos(theta)*cos(phi);
(*vcart[2])[j] =  (*vsph[0])[j]*  sin(theta)
	       +  (*vsph[1])[j]*r*cos(theta);
}
}
else if ( vtype == GVECTYPE_CONTRAVAR ) {
for ( GSIZET j=0; j<nxy; j++ ) {
x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
r     = sqrt(x*x + y*y + z*z);
theta = asin(z/r);
phi   = atan2(y,x);
vthcontra  = (*vsph[0])[j];
vphicontra = (*vsph[1])[j];
(*vcart[0])[j] =  (*vsph[0])[j]*  cos(theta)*cos(phi)
	       -  (*vsph[1])[j]*r*sin(theta)*cos(phi) 
	       +  (*vsph[2])[j]*r*cos(theta)*sin(phi);
(*vcart[1])[j] =  (*vsph[0])[j]*  cos(theta)*sin(phi)
	       -  (*vsph[1])[j]*r*sin(theta)*sin(phi) 
	       +  (*vsph[2])[j]*r*cos(theta)*cos(phi);
(*vcart[2])[j] =  (*vsph[0])[j]*  sin(theta)
	       +  (*vsph[1])[j]*r*cos(theta);
}
}
}
*/

} // end of method vcart2sphere


//**********************************************************************************
//**********************************************************************************
// METHOD : cart2latlon 
// DESC   : Convert Cartesian position vectors to lat long
//
// ARGS   : grid   : Grid. If not of the correct type, nothing is done
//          cart   : x, y, z coords; only first 3 vectors are read
//          latlon : corresponding vectors of (lat,lon); only first 2 vectors 
//                   are written
// RETURNS: none
//**********************************************************************************
template<>
void cart2latlon(const GTVector<GTVector<GFTYPE>*> &cart, GTVector<GTVector<GFTYPE>*> &latlon)
{
  GEOFLOW_TRACE();
assert( cart.size() >= 3 && latlon.size() >= 2 && "Must have correct array sizes");


GSIZET           nxy = cart[0]->size();
GFTYPE           x, y, z;
GFTYPE           phi, r, theta;

for ( GSIZET j=0; j<nxy; j++ ) {
x = (*cart[0])[j]; y = (*cart[1])[j]; z = (*cart[2])[j];
r     = sqrt(x*x + y*y + z*z);
theta = asin(z/r);
phi   = atan2(y,x);
(*latlon[0])[j] = theta; // lat
(*latlon[1])[j] = phi;   // lon
}

} // end of method cart2latlon


//**********************************************************************************
//**********************************************************************************
// METHOD : rcart2sphere
// DESC   : Convert Cartesian position vectors to radius-lat-longcoords
//
// ARGS     cart    : x, y, z coords; only first 3 vectors are read
//          rlatlon : corresponding vectors of radius/lat/lon; only first 3 
//                   vectors are written
// RETURNS: none
//**********************************************************************************
template<>
void rcart2sphere(const GTVector<GTVector<GFTYPE>*> &cart, GTVector<GTVector<GFTYPE>*> &rlatlon)
{
  GEOFLOW_TRACE();
assert( cart.size() >= 3 && rlatlon.size() >= 3 && "Must have correct array sizes");


GSIZET           nxy = cart[0]->size();
GFTYPE           x, y, z;
GFTYPE           phi, r, theta;

for ( GSIZET j=0; j<nxy; j++ ) {
x = (*cart[0])[j]; y = (*cart[1])[j]; z = (*cart[2])[j];
r     = sqrt(x*x + y*y + z*z);
theta = asin(z/r);
phi   = atan2(y,x);
(*rlatlon[0])[j] = r;
(*rlatlon[1])[j] = theta;
(*rlatlon[2])[j] = phi;
}

} // end of method rcart2sphere


//**********************************************************************************
//**********************************************************************************
// METHOD : energy
// DESC   : 
//             Compute volume-integrated mean of energy from input 
//             vector field:
//               energy = 0.5 * Int _u_^2 dV / Int dV
//          
// ARGS   : 
//          grid    : grid object
//          u       : vector field; entire field used to compute energy
//          tmp     : tmp vector of length at least 2, each
//                    of same length as u
//          isglobal: do global reduction
//          ismax   : do max norm, instead of L2
// RETURNS: GFTYPE energy
//**********************************************************************************
template<>
GFTYPE energy(GGrid &grid, const GTVector<GTVector<GFTYPE>*> & u, GTVector<GTVector<GFTYPE>*> &tmp, GBOOL isglobal, GBOOL ismax)
{
    GEOFLOW_TRACE();
GDOUBLE                     ener, local;
GC_COMM                     comm = grid.get_comm();

// Find _u_^2 = Sum_l u_l ^2
*tmp[1] = *u[0]; tmp[1]->rpow(2);
for ( GINT l=1; l<u.size(); l++ ) {
*tmp[0] = *u[l]; tmp[0]->rpow(2);
*tmp[1] += *tmp[0];
}

if ( ismax ) {
ener =  0.5*tmp[1]->amax();
if ( isglobal ) {
local = ener;
GComm::Allreduce(&local, &ener, 1, T2GCDatatype<GDOUBLE>() , GC_OP_MAX, comm);
}
}
else {
ener  = static_cast<GDOUBLE>(grid.integrate(*tmp[1], *tmp[0], isglobal));
ener *=  0.5*static_cast<GDOUBLE>(grid.ivolume());
}

return static_cast<GFTYPE>(ener);

} // end of method energy


//**********************************************************************************
//**********************************************************************************
// METHOD : enstrophy
// DESC   : 
//             Compute volume-integrated mean 
//                 1/2 Int |curl u |^2 dV / Int dV
//          
// ARGS   : 
//          grid    : grid object
//          u       : vector field; entire field used to compute energy
//          tmp     : tmp vector of length at least 4, each
//                    of same length as u
//          isglobal: do global reduction
//          ismax   : if TRUE, then compute max of integrand, and return, 
//                    instead of computing mean
// RETURNS: GFTYPE enstrophy
//**********************************************************************************
template<>
GFTYPE enstrophy(GGrid &grid, const GTVector<GTVector<GFTYPE>*> & u, GTVector<GTVector<GFTYPE>*> &tmp, GBOOL isglobal, GBOOL ismax)
{
    GEOFLOW_TRACE();
assert(tmp.size() >= 4 && "Insufficient temp space");


GINT                        ibeg, iend;
GDOUBLE                     enst, local;
GC_COMM                     comm = grid.get_comm();
GTVector<GFTYPE>           *cc;
GTVector<GTVector<GFTYPE>*> utmp(3);

utmp[0] = tmp[0];
utmp[1] = tmp[1];
cc      = tmp[2];

*tmp[3] = 0.0;
if ( u.size() == 3 ) {
for ( GINT l=0; l<u.size(); l++ ) {
GMTK::curl(grid, u, l+1, utmp, *cc);
cc->rpow(2);
*tmp[3] += *cc;
}
}
else if ( u.size() == 2 ) {
GMTK::curl(grid, u, 3, utmp, *cc);
cc->rpow(2);
*tmp[3] += *cc;
}

if ( ismax ) {
enst =  0.5*static_cast<GDOUBLE>(tmp[3]->amax());
if ( isglobal ) {
local = enst;
GComm::Allreduce(&local, &enst, 1, T2GCDatatype<GDOUBLE>() , GC_OP_MAX, comm);
}
}
else {
enst  = static_cast<GDOUBLE>(grid.integrate(*tmp[3], *tmp[0], isglobal));
enst *=  0.5*grid.ivolume();
}

return static_cast<GFTYPE>(enst);

} // end of method enstrophy



//**********************************************************************************
//**********************************************************************************
// METHOD : helicity
// DESC   : 
//             Compute volume-integrated mean 
//                 Int |curl u \dot u| dV / Int dV
//          
// ARGS   : 
//          grid    : grid object
//          u       : vector field; entire field used to compute energy
//          tmp     : tmp vector of length at least 4, each
//                    of same length as u
//          isglobal: do global reduction
//          ismax   : if TRUE, then compute abs max of integrand, and return, 
//                    instead of computing mean
// RETURNS: GFTYPE helicity
//**********************************************************************************
template<>
GFTYPE helicity(GGrid &grid, const GTVector<GTVector<GFTYPE>*> & u, GTVector<GTVector<GFTYPE>*> &tmp, GBOOL isglobal, GBOOL ismax)
{
    GEOFLOW_TRACE();
assert(tmp.size() >= 4 && "Insufficient temp space");


GDOUBLE                     hel, local;
GC_COMM                     comm = grid.get_comm();
GTVector<GFTYPE>           *cc;
GTVector<GTVector<GFTYPE>*> utmp(3);

utmp[0] = tmp[0];
utmp[1] = tmp[1];
cc      = tmp[2];

*tmp[3] = 0.0;
if ( u.size() == 3 ) {
for ( GINT l=0; l<3; l++ ) {
GMTK::curl(grid, u, l+1, utmp, *cc);
if ( u.size() > l ) cc->pointProd(*u[l]);
*tmp[3] += *cc;
}
}

if ( ismax ) {
hel =  static_cast<GDOUBLE>(tmp[3]->amax());
if ( isglobal ) {
local = hel;
GComm::Allreduce(&local, &hel, 1, T2GCDatatype<GDOUBLE>() , GC_OP_MAX, comm);
}
}
else {
hel  = static_cast<GDOUBLE>(grid.integrate(*tmp[3], *tmp[0], isglobal));
hel *=  static_cast<GDOUBLE>(grid.ivolume());
}

return static_cast<GFTYPE>(hel);

} // end of method helicity 



//**********************************************************************************
//**********************************************************************************
// METHOD : relhelicity
// DESC   : 
//             Compute volume-integrated mean 
//                 Int |curl u \dot u| /(|u| {curl u|) dV / Int dV
//          
// ARGS   : 
//          grid    : grid object
//          u       : vector field; entire field used to compute energy
//          tmp     : tmp vector of length at least 5, each
//                    of same length as u
//          isglobal: do global reduction
//          ismax   : if TRUE, then compute abs max of integrand, and return, 
//                    instead of computing mean
// RETURNS: GFTYPE helicity
//**********************************************************************************
template<>
GFTYPE relhelicity(GGrid &grid, const GTVector<GTVector<GFTYPE>*> & u, GTVector<GTVector<GFTYPE>*> &tmp, GBOOL isglobal, GBOOL ismax)
{
    GEOFLOW_TRACE();
assert(tmp.size() >= 5 && "Insufficient temp space");


GDOUBLE                     local, rhel;
GC_COMM                     comm = grid.get_comm();
GTVector<GFTYPE>           *cc;
GTVector<GTVector<GFTYPE>*> utmp(3);

utmp[0] = tmp[0];
utmp[1] = tmp[1];
cc      = tmp[2];

*tmp[3] = 0.0;
*tmp[4] = 0.0;

// Compute u. curl u:
if ( u.size() == 3 ) {
for ( GINT l=0; l<3; l++ ) {
GMTK::curl(grid, u, l+1, utmp, *cc);
cc->pointProd(*u[l]);
*tmp[3] += *cc;
}
}

// Compute |curl u|:
if ( u.size() == 3 ) {
for ( GINT l=0; l<3; l++ ) {
GMTK::curl(grid, u, l+1, utmp, *cc);
cc->rpow(2);
*tmp[4] += *cc;
}
}
else if ( u.size() == 2 ) {
GMTK::curl(grid, u, 3, utmp, *cc);
cc->rpow(2);
*tmp[4] += *cc;
} 
tmp[4]->rpow(0.5);


// Compute |u|:
*tmp[0] = 0.0;
for ( GINT l=0; l<u.size(); l++ ) {
*tmp[1] = *u[l];
tmp[1]->rpow(2);
*tmp[0] += *tmp[1];
}
tmp[0]->rpow(0.5);


// Compute u. (curl u) /|u| |curl u| integrand:
GFTYPE tiny = 100.0*numeric_limits<GFTYPE>::epsilon();

tmp[0]->pointProd(*tmp[4]); // compute |u| |curl u|
for ( GSIZET k=0; k<utmp[0]->size(); k++ ) {
// (*tmp[1])[k] = fabs((*tmp[0])[k]) <= tiny ? 0.0 : (*tmp[3])[k]/(*tmp[0])[k];  
if ( fabs((*tmp[0])[k]) <= tiny ) {
(*tmp[1])[k] = 0.0;
}
else {
(*tmp[1])[k] = (*tmp[3])[k]/(*tmp[0])[k];
}
}


if ( ismax ) {
rhel =  static_cast<GDOUBLE>(tmp[1]->amax());
if ( isglobal ) {
local = rhel;
GComm::Allreduce(&local, &rhel, 1, T2GCDatatype<GDOUBLE>() , GC_OP_MAX, comm);
}
}
else {
rhel  = static_cast<GDOUBLE>(grid.integrate(*tmp[1], *tmp[0], isglobal));
rhel *= static_cast<GDOUBLE>(grid.ivolume());
}

return static_cast<GFTYPE>(rhel);

} // end of method relhelicity 



//**********************************************************************************
//**********************************************************************************
// METHOD : energyinj
// DESC   : 
//             Compute volume-integrated mean energy injection
//                 Int |u \dot f| dV / Int dV
//          
// ARGS   : 
//          grid    : grid object
//          u       : velocity field; entire field used to compute energy
//          uf      : forcing field; each must be non-NULL on entry; if
//                 any are NULL, or uf.size == 0, then return value is 0.
//          tmp     : tmp vector of length at least 2, each
//                    of same length as u
//          isglobal: do global reduction
//          ismax   : if TRUE, then compute abs max of integrand, and return, 
//                    instead of computing mean
// RETURNS: GFTYPE energy injection rate
//**********************************************************************************
template<>
GFTYPE energyinj(GGrid &grid, const GTVector<GTVector<GFTYPE>*> &u,  const GTVector<GTVector<GFTYPE>*> &uf, GTVector<GTVector<GFTYPE>*> &tmp, GBOOL isglobal, GBOOL ismax)
{
  GEOFLOW_TRACE();
if ( uf.size() == 0 ) return 0.0;

GBOOL bnull = FALSE;
for ( GINT l=0; l<u.size(); l++ ) bnull = bnull || uf[l] == NULLPTR;

if ( bnull ) return 0.0;

assert(tmp.size() >= 2 && "Insufficient temp space");


GDOUBLE                     einj, local;
GC_COMM                     comm = grid.get_comm();


u[0]->pointProd(*uf[0], *tmp[0]);
for ( GINT l=1; l<u.size(); l++ ) {
  assert(uf[l]!= NULLPTR && "NULL force not allowed");
u[l]->pointProd(*uf[l], *tmp[1]);
*tmp[0] += *tmp[1];
}

if ( ismax ) {
einj =  static_cast<GDOUBLE>(tmp[0]->amax());

if ( isglobal ) {
local = einj;
GComm::Allreduce(&local, &einj, 1, T2GCDatatype<GDOUBLE>() , GC_OP_MAX, comm);
}

}
else {
einj  = static_cast<GDOUBLE>(grid.integrate(*tmp[0], *tmp[1], isglobal));
einj *=  static_cast<GDOUBLE>(grid.ivolume());
}

return static_cast<GFTYPE>(einj);

} // end of method energyinj



//**********************************************************************************
//**********************************************************************************
// METHOD : domathop
// DESC   : 
//             Carry out math (mainly differential) operation, based on 
//             specified string description.
//
//             Return in appropriate component of output array
//          
// ARGS   : 
//          grid    : grid object
//          uin     : input state. May be vector or scalar, in an
//                    order and number appropriate for specified operation
//          sop     : string operation
//          tmp     : tmp vector of length at least 1, each
//                    of same length as uin
//          uout    : output array
//          iuout   : which indices of uout contain valid output data; 
//                    sized to reflect valid number of components
// RETURNS: none
//**********************************************************************************
template<>
void domathop(GGrid &grid, const GTVector<GTVector<GFTYPE>*> &uin,  const GString sop, GTVector<GTVector<GFTYPE>*> &utmp, GTVector<GTVector<GFTYPE>*> &uout, GTVector<GINT> &iuout)
{
  GEOFLOW_TRACE();
  GINT                        ib, nxy;
  GTVector<GTVector<GFTYPE>*> tmp(3);

  if      ( "div"  == sop ) { // operates on a vector field...
    // produces a scalar field:
    nxy = grid.gtype() == GE_2DEMBEDDED ? 3 : GDIM;
    assert(utmp .size() >= 1   && "Insufficient temp space");
    assert(uin  .size() >= nxy && "Insufficient no. input components");
    assert(uout .size() >= 1   && "Incorrect no. output components");
    GMTK::grad(grid, *uin[0], 1, utmp, *uout[0]);
    iuout.resize(1); iuout[0] = 0;
    for ( auto j=1; j<nxy; j++ ) {
      GMTK::grad(grid, *uin[j], j+1, utmp, *utmp[utmp.size()-1]);
      *uout[0] += *utmp[utmp.size()-1];
    }
  }
  else if ( "grad" == sop ) {  // operates on a scalar field...
    // produces a vector field:
  nxy = grid.gtype() == GE_2DEMBEDDED ? 3 : GDIM;
    assert(utmp .size() >= 1   && "Insufficient temp space");
    assert(uin  .size() >= 1   && "Incorrect no. input components");
    assert(uout .size() >= nxy && "Insufficient no. output components");
    iuout.resize(nxy); 
      for ( auto j=0; j<nxy; j++ ) {
      GMTK::grad(grid, *uin[0], j+1, utmp, *uout[j]);
      iuout[j] = j; 
    }
  }
  else if ( "gradmag" == sop ) {  // operates on a scalar field...
    // produces a scalar field:
    nxy = grid.gtype() == GE_2DEMBEDDED ? 3 : GDIM;
    assert(utmp .size() >= 2   && "Insufficient temp space");
    assert(uin  .size() >= 1   && "Incorrect no. input components");
    assert(uout .size() >= 1   && "Insufficient no. output components");
    iuout.resize(1); iuout[0] = 0; 
    tmp[0] = utmp[0]; tmp[1] = utmp[1];
    GMTK::grad(grid, *uin[0], 1, tmp, *uout[0]);
    uout[0]->rpow(2);
    for ( auto j=1; j<nxy; j++ ) {
      GMTK::grad(grid, *uin[0], j+1, tmp, *tmp[1]);
      tmp[1]->rpow(2);
      *uout[0] += *tmp[1];
    }
    uout[0]->rpow(0.5);
  }
  else if ( "curl" == sop ) { // operates on a vector field...
    // produces a vector field:
    nxy = grid.gtype() == GE_2DEMBEDDED ? 3 : GDIM;
    ib = 1; 
    if ( GDIM == 2 && grid.gtype() != GE_2DEMBEDDED ) { nxy = 1; ib = 3; }
    assert(utmp .size() >= 2   && "Insufficient temp space");
    assert(uin  .size() >= nxy && "Insufficient no. input components");
    assert(uout .size() >= nxy && "Insufficient no. output components");
    iuout.resize(nxy); 
    for ( auto j=0; j<nxy; j++ ) {
      GMTK::curl(grid, uin, j+ib, utmp, *uout[j]);
      iuout[j] = j; 
    }
  }
  else if ( "curlmag" == sop ) { // operates on a vector field...
    // produces a scalar field:
    nxy = grid.gtype() == GE_2DEMBEDDED ? 3 : GDIM;
    ib = 1; 
    if ( GDIM == 2 && grid.gtype() != GE_2DEMBEDDED ) { nxy = 1; ib = 3; }
    assert(utmp .size() >= 3   && "Insufficient temp space");
    assert(uin  .size() >= nxy && "Insufficient no. input components");
    assert(uout .size() >= 1   && "Insufficient no. output components");
    iuout.resize(1); iuout[0] = 0; 
    tmp[0] = utmp[0]; tmp[1] = utmp[1]; tmp[2] = utmp[2];
    GMTK::curl(grid, uin, ib, utmp, *uout[0]);
    uout[0]->rpow(2);
    for ( auto j=1; j<nxy; j++ ) {
      GMTK::curl(grid, uin, j+ib, tmp, *tmp[2]);
      tmp[2]->rpow(2);
      *uout[0] += *tmp[2];
    }
    uout[0]->rpow(0.5);
  }
  else if ( "vmag" == sop ) { // operates on a vector field...
    // produces a scalar field:
    nxy = grid.gtype() == GE_2DEMBEDDED ? 3 : GDIM;
    assert(utmp .size() >= 1   && "Insufficient temp space");
    assert(uin  .size() >= nxy && "Insufficient no. input components");
    assert(uout .size() >= 1   && "Insufficient no. output components");
    iuout.resize(1); iuout[0] = 0; 
    tmp[0] = utmp[0]; ;
    *uout[0] = *uin[0]; uout[0]->rpow(2);
    for ( auto j=1; j<nxy; j++ ) {
      *tmp[0] = *uin[j]; tmp[0]->rpow(2);
      *uout[0] += *tmp[0];
    }
    uout[0]->rpow(0.5);
  }
  else if ( "lapderivs" == sop ) { // operates on a scalar field
    // produces a vector field of each termin Laplacnin:
    nxy = grid.gtype() == GE_2DEMBEDDED ? 3 : GDIM;
    assert(utmp .size() >= 3   && "Insufficient temp space");
    assert(uin  .size() >= 1   && "Insufficient no. input components");
    assert(uout .size() >= nxy && "Insufficient no. output components");
    for ( auto j=0; j<3; j++ ) tmp[j] = utmp[j];
    iuout.resize(nxy); 
    for ( auto j=0; j<nxy; j++ ) {
      GMTK::grad(grid, *uin[0], j+1, tmp , *tmp[2]);
      GMTK::grad(grid, *tmp[2], j+1, utmp, *uout[j]);
      iuout[j] = j; 
    }
  }
  else {
  assert(FALSE && "Invalid math operation");
  }

} // end of method domathop 


} // end, namespace GMTK
