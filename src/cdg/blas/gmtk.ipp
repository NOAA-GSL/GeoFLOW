//==================================================================================
// Module       : gmtk.ipp
// Date         : 1/31/18 (DLR)
// Description  : Math TooKit: namespace of C routines for various
//                math functions
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================


#include "tbox/tracer.hpp"


namespace GMTK 
{



//**********************************************************************************
//**********************************************************************************
// METHOD : fact
// DESC   : Compute factorial of specified number
//
// ARGS   : n  : whole number
// RETURNS: T factorial
//**********************************************************************************
template<typename T>
T fact(T n)
{
	GEOFLOW_TRACE();
   T Tzero = static_cast<T>(0);
   T Tone  = static_cast<T>(1);
   if ( n == Tzero || n == Tone ) {
     return Tone;
   }
   else {
     return n*fact(n-Tone);
   }

} // end, method fact


//**********************************************************************************
//**********************************************************************************
// METHOD : Plm_cart
// DESC   : Compute assosiated Legendre polynomials at each grid point
//          specified by xnodes in Carteswian coordinates.
//              Note: Adapted from Numerical Recipes
//
// ARGS   : l,m,  : orbital ang mom. quantum number, and azimuthal quantum number
//          xnodes: Cartesian coordinate arrays
//          plm   : P^l_m for each grid point
// RETURNS: none
//**********************************************************************************
template<typename T>
void Plm_cart(GINT l, GINT m, GTVector<GTVector<T>> &xnodes, GTVector<T> &plm)
{
	GEOFLOW_TRACE();
  T fact;
  T pmm, pmmp1, pll, somx2;
  T colat, r;
  T xlmm, xlpm;
  T x, y, z;

  assert( m >= 0 && m <= l );

  for ( auto j=0; j<xnodes[0].size(); j++ ) {
    x     = xnodes[0][j]; y = xnodes[1][j]; z = xnodes[2][j];
    r     = sqrt(x*x +y*y + z*z);

    colat = acos(z/r); // Note: this is co-latitude
    x     = cos(colat);
    assert( fabs(x) <= 1.0  );

    pmm = 1.0; // compute Pm_m 
    if ( m > 0 ) {
      somx2 = sqrt(fabs((1.0-x)*(1.0+x))); 
      fact  = 1.0;
      for ( auto i=1; i<=m; i++) {
        pmm  *= -fact*somx2;
        fact += 2.0; 
      }
    }
    if (l == m) {
      plm[j] = pmm; 
      continue;
    }
    else {   // compute P^m_m+1    
      pmmp1 = x*(2*m+1)*pmm; 
      if (l == (m+1)) {
        plm[j] = pmmp1; 
        continue;
      }
      else { // compute P^m_l; l > m+1
        for ( auto ll=m+2; ll<=l; ll++) {
          xlpm  = ll + m - 1;
          xlmm  = ll - m;
          pll   = (x*(2*ll-1)*pmmp1 - xlpm*pmm)/xlmm; 
          pmm   = pmmp1;
          pmmp1 = pll;
        }
        plm[j] = pll; 
        continue;
      }
    }
  } // end, coord loop

} // end, method Plm_cart


//**********************************************************************************
//**********************************************************************************
// METHOD : Ylm_cart
// DESC   : Compute spherical harmonics at each grid point
//          specified by xnodes in Carteswian coordinates.
//          Computed as:
//              Ylm  =  (-1)^m sqrt{ (2l+1)(l-m)! / [4pi (l+m)!] } X
//                      P_lm exp(i m phi)
//          If iri = 0, return Ylm; else returns Ylm^*. If |m| > l, set Ylm=0
//          
//   
// ARGS   : l,m,  : orbital ang mom. quantum number, and azimuthal quantum number
//          xnodes: Cartesian coordinate arrays
//          iri   : if 0, return Ylm, else, return complex conjugate, Ylm^*
//          ylm_r: real comp. of Y^l_m for each grid point
//          ylm_i: imag. comp. of Y^l_m for each grid point
// RETURNS: none
//**********************************************************************************
template<typename T>
void Ylm_cart(GINT l, GINT m, GTVector<GTVector<T>> &xnodes, GINT iri, GTVector<T> &ylm_r, GTVector<T> &ylm_i)
{
	GEOFLOW_TRACE();
  T phi;
  T x, y, z;
  GDOUBLE rfact, xl, xm;
  GDOUBLE xf1, xf2;

  assert( iri >= 0 );
  assert( l >= 0 );

  if ( abs(m) > l ) {
    ylm_r = 0.0;
    ylm_i = 0.0;
    return;
  }
  
  xl = l; xm = m;
  xf1 = GMTK::fact<GDOUBLE>(xl-abs(xm));
  xf2 = GMTK::fact<GDOUBLE>(xl+abs(xm));
  rfact = sqrt( fabs(2.0*xl+1.0)/(4.0*PI) * (xf1 / xf2) ); 

  GMTK::Plm_cart<T>(l, abs(m), xnodes, ylm_r);
  ylm_i = ylm_r;
  for ( auto j=0; j<xnodes[0].size(); j++ ) {
    x     = xnodes[0][j]; y = xnodes[1][j]; z = xnodes[2][j];
    phi   = atan2(y,x);
    ylm_r[j] *= pow(-1.0,m)*rfact*cos(xm*phi);
    ylm_i[j] *= pow(-1.0,m)*rfact*sin(xm*phi);
  } // end, coord loop

  if ( iri > 0 ) { // create complex conjugate
    ylm_i *= -1.0;
  }

} // end, method Ylm_cart


//**********************************************************************************
//**********************************************************************************
// METHOD : dYlm_cart
// DESC   : Compute 1-derivatives of spherical harmonics at each grid point
//          specified by xnodes in Carteswian coordinates.
//          Cmputed as:
//              dYlm/dtheta = m cot theta Ylm + 
//                            sqrt[(l-m)(l+m+1)] e^(-i phi) Yl(m+1)
//              dYlm/dphi   = i m Ylm 
//          If iri = 0, return dYlm; else return dYlm^*. 
//          
//   
// ARGS   : l,m,  : orbital ang mom. quantum number, and azimuthal quantum number
//          xnodes: Cartesian coordinate arrays
//          idir  : If 1, take d/d\theta; if 2, take d/d\phi
//          iri   : if 0, return dYlm, else, return complex conjugate, dYlm^*
//          cexcl : exclude colat < cexcl and (PI-colat) < cexcl from computation 
//                  (set to 0). Expressed in degrees.
//          tmp   : tmp arrays. 2 are required.
//          dylm_r: real comp. of dY^l_m for each grid point
//          dylm_i: imag. comp. of dY^l_m for each grid point
// RETURNS: none
//**********************************************************************************
template<typename T>
void dYlm_cart(GINT l, GINT m, GTVector<GTVector<T>> &xnodes, GINT idir,  GINT iri, T cexcl,  GTVector<GTVector<T>*> &tmp, GTVector<T> &dylm_r, GTVector<T> &dylm_i)
{
	GEOFLOW_TRACE();
  T colat, phi, r, rexcl;
  T x, y, z;
  T xl, xm;
  T cotc, rfact, rfact1, rfact2;

  rexcl = cexcl * PI/180.0;

  assert( iri >= 0 );
  assert( idir == 1 || idir == 2 );

  xl    = l; xm = m;

  if ( idir == 1 ) {

#if 0
    // Compute: 
    // dYlm/dtheta = m cot theta Ylm + 
    //               sqrt[(l-m)(l+m+1)] e^(-i phi) Yl(m+1)
    assert( tmp.size() >= 2 );
    rfact = sqrt( fabs((xl-xm)*(xl+xm+1.0)) );
    Ylm_cart<T>(l, m  , xnodes, 0, dylm_r  , dylm_i); // Ylm
    Ylm_cart<T>(l, m+1, xnodes, 0, *tmp[0], *tmp[1]); // Y^l_m+1
    for ( auto j=0; j<xnodes[0].size(); j++ ) {
      x     = xnodes[0][j]; y = xnodes[1][j]; z = xnodes[2][j];
      r     = sqrt(x*x + y*y + z*z);
      colat = acos(z/r);
      if ( colat < rexcl || (PI-colat) < rexcl ) {
        dylm_r[j] = dylm_i[j] = 0.0;
        continue;
      }
      cotc  = cos(colat)/sin(colat);
      phi   = atan2(y,x);
      dylm_r[j] = xm*cotc*dylm_r[j] 
                + rfact*( cos(phi)*(*tmp[0])[j] + sin(phi)*(*tmp[1])[j] );
      dylm_i[j] = xm*cotc*dylm_i[j]
                - rfact*( cos(phi)*(*tmp[1])[j] - sin(phi)*(*tmp[0])[j] );
    } // end, coord loop
#else
    // Compute: 
    // dYlm/dtheta = exp(-i phi) sqrt((l-m)(l+m+1)) Y^l_m+1
    //             - exp (i phi) sqrt*(l+m)(l-m+1)) Y^l_m-1
    assert( tmp.size() >= 2 );
#if 0
    rfact1 = 0.5*sqrt( fabs((xl-xm)*(xl+xm+1.0)) );
    rfact2 = 0.5*sqrt( fabs((xl+xm)*(xl-xm+1.0)) );
#else
    rfact1 = sqrt( fabs((xl-xm)*(xl+xm+1.0)) );
    rfact2 = sqrt( fabs((xl+xm)*(xl-xm+1.0)) );
#endif
    Ylm_cart<T>(l, m+1, xnodes, 0, *tmp[0], *tmp[1]); // Y^l_m+1
    for ( auto j=0; j<xnodes[0].size(); j++ ) {
      x     = xnodes[0][j]; y = xnodes[1][j]; z = xnodes[2][j];
      r     = sqrt(x*x + y*y + z*z);
      phi   = atan2(y,x);
      dylm_r[j] = rfact1*(cos(phi)*(*tmp[0])[j] + sin(phi)*(*tmp[1])[j]);
      dylm_i[j] = rfact1*(cos(phi)*(*tmp[1])[j] - sin(phi)*(*tmp[0])[j]);
    } // end, coord loop

    Ylm_cart<T>(l, m-1, xnodes, 0, *tmp[0], *tmp[1]); // Y^l_m-1
    for ( auto j=0; j<xnodes[0].size(); j++ ) {
      x     = xnodes[0][j]; y = xnodes[1][j]; z = xnodes[2][j];
      r     = sqrt(x*x + y*y + z*z);
      phi   = atan2(y,x);
      dylm_r[j] -= rfact2*(cos(phi)*(*tmp[0])[j] - sin(phi)*(*tmp[1])[j]);
      dylm_i[j] -= rfact2*(cos(phi)*(*tmp[1])[j] + sin(phi)*(*tmp[0])[j]);
    } // end, coord loop
#endif

  }
  else if ( idir == 2 ) {
    // Compute:
    // dYlm/dphi   = i m Ylm 
    Ylm_cart<T>(l, m, xnodes, 0, *tmp[0], *tmp[1]);   // Ylm
    for ( auto j=0; j<xnodes[0].size(); j++ ) {
      dylm_r[j] = -xm * (*tmp[1])[j];
      dylm_i[j] =  xm * (*tmp[0])[j];
    } // end, coord loop
  

  }

  if ( iri > 0 ) { // create complex conjugate
    dylm_i *= -1.0;
  }

} // end, method dYlm_cart


//**********************************************************************************
//**********************************************************************************
// METHOD : ddYlm_cart
// DESC   : Compute 2-derivatives of spherical harmonics at each grid point
//          specified by xnodes in Carteswian coordinates.
//          If iri = 0, return ddYlm; else return ddYlm^*. 
//          
//   
// ARGS   : l,m,  : orbital ang mom. quantum number, and azimuthal quantum number
//          xnodes: Cartesian coordinate arrays
//          idir  : If 1, take d^2/d\theta^2; if 2, take d^2/d\phi^2
//          iri   : if 0, return dYlm, else, return complex conjugate, dYlm^*
//          cexcl : exclude colat < cexcl and (PI-colat) < cexcl from computation 
//                  (set to 0). Expressed in degrees.
//          tmp   : tmp arrays. 2 are required.
//          dylm_r: real comp. of dY^l_m for each grid point
//          dylm_i: imag. comp. of dY^l_m for each grid point
// RETURNS: none
//**********************************************************************************
template<typename T>
void ddYlm_cart(GINT l, GINT m, GTVector<GTVector<T>> &xnodes, GINT idir, GINT iri, T cexcl,  GTVector<GTVector<T>*> &tmp, GTVector<T> &dylm_r, GTVector<T> &dylm_i)
{
	GEOFLOW_TRACE();
  T colat, phi, r, rexcl;
  T x, y, z;
  T xl, xm;
  T cotc, csc, rfact;

  assert( iri >= 0 );
  assert( idir == 1 || idir == 2 );
  assert( tmp.size() >= 2 );

  xl   = l; xm = m;

 
  if ( idir == 1 ) {

    Ylm_cart<T>(l, m  , xnodes, 0, dylm_r  , dylm_i);   // Ylm
    Ylm_cart<T>(l, m+1, xnodes, 0, *tmp[0], *tmp[1]); // Y^l_m+1
    rfact = sqrt( fabs((xl-xm)*(xl+xm-1.0)) ) * (2.0*xm+1.0);
    for ( auto j=0; j<xnodes[0].size(); j++ ) {
      x     = xnodes[0][j]; y = xnodes[1][j]; z = xnodes[2][j];
      r     = sqrt(x*x + y*y + z*z);
      colat = acos(z/r);
      if ( colat < rexcl || (PI-colat) < rexcl ) {
        dylm_r[j] = dylm_i[j] = 0.0;
        continue;
      }
      cotc  = cos(colat)/sin(colat);
      csc   = 1.0 / sin(colat);
      phi   = atan2(y,x);

      dylm_r[j] = xm*( xm*cotc*cotc - csc*csc )*dylm_r[j] 
                + rfact*cotc*( cos(phi)*(*tmp[0])[j] + sin(phi)*(*tmp[1])[j] );
      dylm_i[j] = xm*( xm*cotc*cotc - csc*csc )*dylm_i[j]
                + rfact*cotc*( cos(phi)*(*tmp[1])[j] - sin(phi)*(*tmp[0])[j] );

    } // end, coord loop

    Ylm_cart<T>(l, m+2, xnodes, 0, *tmp[0], *tmp[1]); // Y^l_m+2
    rfact = sqrt( fabs((xl-xm)*(xl-xm-1.0)*(xm+xl+2.0)*(xm+xl+1.0)) );
    for ( auto j=0; j<xnodes[0].size(); j++ ) {
      x     = xnodes[0][j]; y = xnodes[1][j]; z = xnodes[2][j];
      phi   = atan2(y,x);
      dylm_r[j] += rfact*( cos(2.0*phi)*(*tmp[0])[j] + sin(2.0*phi)*(*tmp[1])[j] );
      dylm_i[j] += rfact*( cos(2.0*phi)*(*tmp[1])[j] - sin(2.0*phi)*(*tmp[0])[j] );
    } // end, coord loop

  }
  else if ( idir == 2 ) {

    Ylm_cart<T>(l, m, xnodes, 0, *tmp[0], *tmp[1]);   // Ylm
    for ( auto j=0; j<xnodes[0].size(); j++ ) {
      dylm_r[j] = -xm*xm * (*tmp[0])[j];
      dylm_i[j] = -xm*xm * (*tmp[1])[j];
    } // end, coord loop
  

  }

  if ( iri > 0 ) { // create complex conjugate
    dylm_i *= -1.0;
  }

} // end, method ddYlm_cart


//**********************************************************************************
//**********************************************************************************
// METHOD : rYlm_cart
// DESC   : Compute real spherical harmonics at each grid point
//          specified by xnodes in Carteswian coordinates.
//          Computed as:
//              rYlm  =  (-1)^m/sqrt(2) (Ylm + Ylm^*); m > 0
//                                       Yl0         ; m = 0
//                     -i(-1)^m/sqrt(2) (Ylm - Ylm^*); m < 0
//           Adapted from 
//           Miguel A. Blanco, et al., 
//             J. Mol. Struct. Theochem 419:19 (1997)
//   
// ARGS   : l,m,  : orbital ang mom. quantum number, and azimuthal quantum number
//          xnodes: Cartesian coordinate arrays
//          tmp   : tmp array the size of ylm
//          rylm  : real comp. of Y^l_m for each grid point
// RETURNS: none
//**********************************************************************************
template<typename T>
void rYlm_cart(GINT l, GINT m, GTVector<GTVector<T>> &xnodes, GTVector<T> &tmp, GTVector<T> &rylm)
{
	GEOFLOW_TRACE();
  T rfact;


  if ( abs(m) > l ) {
    rylm = 0.0;
    return;
  }

  if ( m > 0 ) { 
    rfact = pow(-1.0,m)/sqrt(2.0);
    GMTK::Ylm_cart<T>(l, m, xnodes, 0, rylm, tmp);
    rylm *= 2.0*rfact;
  }
  else if ( m < 0 ) {
    rfact = pow(-1.0,m)/sqrt(2.0);
    GMTK::Ylm_cart<T>(l, m, xnodes, 0, tmp, rylm);
    rylm *= 2.0*rfact;
  }
  else {
    GMTK::Ylm_cart<T>(l, m, xnodes, 0, rylm, tmp);
  }

} // end, method rYlm_cart


//**********************************************************************************
//**********************************************************************************
// METHOD : drYlm_cart
// DESC   : Compute idir--derivative real spherical harmonics at each grid point
//          specified by xnodes in Carteswian coordinates.
//          Computed as:
//             drYlm  =  (-1)^m/sqrt(2) (d_idir Ylm + d_idir Ylm^*); m > 0
//                                       d_idir Yl0                ; m = 0
//                     -i(-1)^m/sqrt(2) (d_idir Ylm - d_idir Ylm^*); m < 0
//          
//   
// ARGS   : l,m,  : orbital ang mom. quantum number, and azimuthal quantum number
//          xnodes: Cartesian coordinate arrays
//          idir  : If 1, computes d/d\theta; if 2, computes d/d\phi
//          cexcl : exclude colat < cexcl and (PI-colat) < cexcl from computation 
//                  (set to 0). Expressed in degrees.
//          tmp   : tmp arrays. 3 are required.
//          dylm  : deriv of real comp. of Y^l_m for each grid point
// RETURNS: none
//**********************************************************************************
template<typename T>
void drYlm_cart(GINT l, GINT m, GTVector<GTVector<T>> &xnodes, GINT idir, T cexcl, GTVector<GTVector<T>*>  &tmp, GTVector<T> &drylm)
{
	GEOFLOW_TRACE();
  T rfact;
  GTVector<T> *dr, *di;  // real and imaginary comps of derivative
  GTVector<GTVector<T>> mytmp(2);

  if ( abs(m) > l ) {
    drylm = 0.0;
    return;
  }

  assert( idir == 1 || idir == 2 );
  assert( tmp.size() >= 3 );

  if ( m > 0 ) { 
    rfact = pow(-1.0,m)/sqrt(2.0);
    di = tmp[2];
    GMTK::dYlm_cart<T>(l, m, xnodes, idir, 0, cexcl, tmp, drylm, *di);
    drylm *= 2.0*rfact;
  }
  else if ( m < 0 ) {
    rfact = pow(-1.0,m)/sqrt(2.0);
    dr = tmp[2];
    GMTK::dYlm_cart<T>(l, m, xnodes, idir, 0, cexcl, tmp, *dr, drylm);
    drylm *= 2.0*rfact;
  }
  else {
    di = tmp[2];
    GMTK::dYlm_cart<T>(l, m, xnodes, idir, 0, cexcl, tmp, drylm, *di);
  }

} // end, method drYlm_cart


//**********************************************************************************
//**********************************************************************************
// METHOD : ddrYlm_cart
// DESC   : Compute idir--2nd derivative of real spherical harmonics at 
//          each grid point specified by xnodes in Carteswian coordinates.
//          Computed as:
//            ddrYlm  =  (-1)^m/sqrt(2) (d^2_idir Ylm + d^2_idir Ylm^*); m > 0
//                                       d^2_idir Yl0                  ; m = 0
//                     -i(-1)^m/sqrt(2) (d^2_idir Ylm - d^2_idir Ylm^*); m < 0
//          
//   
// ARGS   : l,m,  : orbital ang mom. quantum number, and azimuthal quantum number
//          xnodes: Cartesian coordinate arrays
//          idir  : If 1, computes d/d\theta; if 2, computes d/d\phi
//          cexcl : exclude colat < cexcl and (PI-colat) < cexcl from computation 
//                  (set to 0). Expressed in degrees.
//          tmp   : tmp arrays. 3 are required.
//          dylm  : deriv. of real comp. of Y^l_m for each grid point
// RETURNS: none
//**********************************************************************************
template<typename T>
void ddrYlm_cart(GINT l, GINT m, GTVector<GTVector<T>> &xnodes, GINT idir, T cexcl, GTVector<GTVector<T>*>  &tmp, GTVector<T> &drylm)
{
	GEOFLOW_TRACE();
  T rfact;
  GTVector<T> *dr, *di;  // real and imaginary comps of derivative
  GTVector<GTVector<T>> mytmp(2);

  if ( abs(m) > l ) {
    drylm = 0.0;
    return;
  }

  assert( idir == 1 || idir == 2 );
  assert( tmp.size() >= 3 );

  if ( m > 0 ) { 
    rfact = pow(-1.0,m)/sqrt(2.0);
    di = tmp[2];
    GMTK::ddYlm_cart<T>(l, m, xnodes, idir, 0, cexcl, tmp, drylm, *di);
    drylm *= 2.0*rfact;
  }
  else if ( m < 0 ) {
    rfact = pow(-1.0,m)/sqrt(2.0);
    dr = tmp[2];
    GMTK::ddYlm_cart<T>(l, m, xnodes, idir, 0, cexcl, tmp, *dr, drylm);
    drylm *= 2.0*rfact;
  }
  else {
    di = tmp[2];
    GMTK::ddYlm_cart<T>(l, m, xnodes, idir, 0, cexcl, tmp, drylm, *di);
  }

} // end, method ddrYlm_cart


//**********************************************************************************
//**********************************************************************************
// METHOD : Rx3
// DESC   : Compute finite rotation of amount alpha about x axis, and return
//          in the input arrays
//             x_new  = Rx x_old,
//          where 
//               | 1    0    0     |
//          Rx = | 0   cos a sin a |
//               | 0  -sin a cos a |
//          Note: This is a 'passive' rotation of coord system. User must
//                late alpha --> -alpha to get an active rotation of a vector
//          
// ARGS   : alpha : angle of rotation
//          y, z  : arrays of vector components
// RETURNS: none.
//**********************************************************************************
template<typename T>
void Rx3(T alpha, GTVector<T> &y, GTVector<T> &z)
{
	GEOFLOW_TRACE();
  assert( z.size() == y.size() &&  "Incompatible vectors");

  for ( auto j=0; j<y.size(); j++ ) { // cycle over all vector elems
    y[j]  =  cos(alpha)*y[j] + sin(alpha)*z[j];
    z[j]  = -sin(alpha)*y[j] + cos(alpha)*z[j];
  }

} // end of method Rx3 


//**********************************************************************************
//**********************************************************************************
// METHOD : Ry3
// DESC   : Compute finite rotation of amount alpha about y axis, and return
//          in the input arrays
//             x_new  = Ry x_old,
//          where 
//               | cos a    0    sin a |
//          Rx = | 0        1    0     |
//               | -sin a   0    cos a |
//          Note: This is a 'passive' rotation of coord system. User must
//                late alpha --> -alpha to get an active rotation of a vector
//          
// ARGS   : alpha : angle of rotation
//          x, z  : arrays of vector components
// RETURNS: none.
//**********************************************************************************
template<typename T>
void Ry3(T alpha, GTVector<T> &x, GTVector<T> &z)
{
	GEOFLOW_TRACE();
  assert( z.size() == x.size() &&  "Incompatible vectors");

  for ( auto j=0; j<x.size(); j++ ) { // cycle over all vector elems
    x[j]  =  cos(alpha)*x[j] + sin(alpha)*z[j];
    z[j]  = -sin(alpha)*x[j] + cos(alpha)*z[j];
  }

} // end of method Ry3 


//**********************************************************************************
//**********************************************************************************
// METHOD : Rz3
// DESC   : Compute finite rotation of amount alpha about z axis, and return
//          in the input arrays
//             x_new  = Rz x_old,
//          where 
//               |  cos a sin a  0 |
//          Rz = | -sin a cos a  0 |
//               | 0    0        1 |
//          Note: This is a 'passive' rotation of coord system. User must
//                late alpha --> -alpha to get an active rotation of a vector
//          
// ARGS   : alpha  : angle of rotation
//          x, y   : arrays of vector components
// RETURNS: none.
//**********************************************************************************
template<typename T>
void Rz3(T alpha, GTVector<T> &x, GTVector<T> &y)
{
	GEOFLOW_TRACE();
  assert( y.size() == x.size() &&  "Incompatible vectors");

  for ( auto j=0; j<x.size(); j++ ) { // cycle over all vector elems
    x[j]  =  cos(alpha)*x[j] + sin(alpha)*y[j];
    y[j]  = -sin(alpha)*x[j] + cos(alpha)*y[j];
  }

} // end of method Rz3 


//**********************************************************************************
//**********************************************************************************
// METHOD : cross_prod_k (1)
// DESC   : Compute cross/vector product with hat(k)
//             C = A X hat(k)
//          
// ARGS   : A    : array of pointers to vectors; must each have at least 3
//                 elements for 3-d vector product. All vector elements must
//                 have the same length
//          iind : pointer to indirection array, used if non_NULLPTR.
//          nind : number of indices in iind. This is the number of cross products
//                 that will be computed, if iind is non-NULLPTR
//          isgn : multiplies product by sign(isgn)
//          C    : array of pointers to cross produc vectors; must have
//                 at least 3 elements for 3-d vector products. All vector
//                 elements must have the same length, of at least length nind if
//                 iind != NULLPTR, and of length of x[?], y[?] if iind==NULLPTR.
// RETURNS: GTVector & 
//**********************************************************************************
template<typename T>
void cross_prod_k(GTVector<GTVector<T>*> &A, GINT *iind, GINT nind, GINT isgn, GTVector<GTVector<T>*> &C)
{
	GEOFLOW_TRACE();
  assert( A.size() >= 2 && C.size() >= 2 &&  "Incompatible dimensionality");

  GSIZET n;
  T      fact = isgn < 0 ? -1.0 : 1.0;
  if ( iind != NULLPTR ) {
    for ( auto k=0; k<nind; k++ ) { // cycle over all coord pairs
      n = iind[k];
      C[0][k] =  (*A[1])[n]*fact;
      C[1][k] = -(*A[0])[n]*fact;
      if ( C.size() > 2 ) C[2][n] = 0.0;
    }
  }
  else {
    for ( auto n=0; n<A[0]->size(); n++ ) { // cycle over all coord pairs
      C[0][n] =  (*A[1])[n]*fact;
      C[1][n] = -(*A[0])[n]*fact;
      if ( C.size() > 2 ) C[2][n] = 0.0;
    }
  }

} // end of method cross_prod_k (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : cross_prod_k (2)
// DESC   : Compute cross/vector product with hat(k)
//             C = A X hat(k)
//          
// ARGS   : Ai   : Vector components x, y, z; must each have same no. elements.
//          iind : pointer to indirection array, used if non_NULLPTR.
//          nind : number of indices in iind. This is the number of cross products
//                 that will be computed, if iind is non-NULLPTR
//          isgn : multiplies product by sign(isgn)
//          Ci   : Vector components for solution, must each have the same no. elements
//                 as in Ai, Bi, unless iind != NULLPTR, in which case, Ci must each
//                 have at least nind elements.
// RETURNS: GTVector & 
//**********************************************************************************
template<typename T>
void cross_prod_k(GTVector<T> &Ax, GTVector<T> &Ay, 
                  GINT *iind, GINT nind, GINT isgn, 
                  GTVector<T> &Cx, GTVector<T> &Cy)
{
	GEOFLOW_TRACE();
  GSIZET n;
  T      fact = isgn < 0 ? -1.0 : 1.0;
  if ( iind != NULLPTR ) {
    for ( auto k=0; k<nind; k++ ) { // cycle over all coord pairs
      n = iind[k];
      Cx[k] =  Ay[n]*fact;
      Cy[k] = -Ax[n]*fact;
    }
  }
  else {
    for ( auto n=0; n<Ax.size(); n++ ) { // cycle over all coord pairs
      Cx[n] =  Ay[n]*fact;
      Cy[n] = -Ax[n]*fact;
    }
  }

} // end of method cross_prod_k (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : cross_prod (1)
// DESC   : compute cross/vector product 
//             C = A X B
//          
// ARGS   : A    : array of pointers to vectors; must each have at least 3
//                 elements for 3-d vector product. All vector elements must
//                 have the same length
//          B    : array of pointers to vectors; must each have at least 3
//                 elements for 3-d vector product. All vector elements must
//                 have the same length
//          iind : pointer to indirection array, used if non_NULLPTR.
//          nind : number of indices in iind. This is the number of cross products
//                 that will be computed, if iind is non-NULLPTR
//          C    : array of pointers to cross produc vectors; must have
//                 at least 3 elements for 3-d vector products. All vector
//                 elements must have the same length, of at least length nind if
//                 iind != NULLPTR, and of length of x[?], y[?] if iind==NULLPTR.
// RETURNS: GTVector & 
//**********************************************************************************
template<typename T>
void cross_prod(GTVector<GTVector<T>*> &A, GTVector<GTVector<T>*> &B, 
                GINT *iind, GINT nind, GTVector<GTVector<T>*> &C)
{
	GEOFLOW_TRACE();
  assert( A.size() >= 3 && B.size() && C.size() >= 3 && "Incompatible input vectors");

  GSIZET n;
  T      x1, y1, z1;
  T      x2, y2, z2;

  if ( iind != NULLPTR ) {
    for ( auto k=0; k<nind; k++ ) { // cycle over all coord pairs
      n = iind[k];
      x1 = (*A[0])[n]; y1 = (*A[1])[n]; z1 = (*A[2])[n];
      x2 = (*B[0])[n]; y2 = (*B[1])[n]; z2 = (*B[2])[n];
      C[0][k] = y1*z2 - z1*y2; 
      C[1][k] = z1*x2 - z2*x1; 
      C[2][k] = x1*y2 - x2*y1;
    }
  }
  else {
    for ( auto n=0; n<A[0]->size(); n++ ) { // cycle over all coord pairs
      x1 = (*A[0])[n]; y1 = (*A[1])[n]; z1 = (*A[2])[n];
      x2 = (*B[0])[n]; y2 = (*B[1])[n]; z2 = (*B[2])[n];
      C[0][n] = y1*z2 - z1*y2; 
      C[1][n] = z1*x2 - z2*x1; 
      C[2][n] = x1*y2 - x2*y1;
    }
  }

} // end of method cross_prod (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : cross_prod (2)
// DESC   : compute cross/vector product 
//             C = A X B
//          
// ARGS   : Ai   : Vector components x, y, z; must each have same no. elements.
// ARGS   : Bi   : Vector components x, y, z; must each have same no. elements, as Ai
//                 have the same length
//          iind : pointer to indirection array, used if non_NULLPTR.
//          nind : number of indices in iind. This is the number of cross products
//                 that will be computed, if iind is non-NULLPTR
//          Ci   : Vector components for solution, must each have the same no. elements
//                 as in Ai, Bi, unless iind != NULLPTR, in which case, Ci must each
//                 have at least nind elements.
// RETURNS: GTVector & 
//**********************************************************************************
template<typename T>
void cross_prod(GTVector<T> &Ax, GTVector<T> &Ay, GTVector<T> &Az,
                GTVector<T> &Bx, GTVector<T> &By, GTVector<T> &Bz,
                GINT *iind, GINT nind, 
                GTVector<T> &Cx, GTVector<T> &Cy, GTVector<T> &Cz)
{
	GEOFLOW_TRACE();
  GSIZET n;
  T      x1, y1, z1;
  T      x2, y2, z2;

  if ( iind != NULLPTR ) {
    for ( auto k=0; k<nind; k++ ) { // cycle over all coord pairs
      n = iind[k];
      x1 = Ax[n]; y1 = Ay[n]; z1 = Az[n];
      x2 = Bx[n]; y2 = By[n]; z2 = Bz[n];
      Cx[k] = y1*z2 - z1*y2; 
      Cy[k] = z1*x2 - z2*x1; 
      Cz[k] = x1*y2 - x2*y1;
    }
  }
  else {
    for ( auto n=0; n<Ax.size(); n++ ) { // cycle over all coord pairs
      x1 = Ax[n]; y1 = Ay[n]; z1 = Az[n];
      x2 = Bx[n]; y2 = By[n]; z2 = Bz[n];
      Cx[n] = y1*z2 - z1*y2; 
      Cy[n] = z1*x2 - z2*x1; 
      Cz[n] = x1*y2 - x2*y1;
    }
  }

} // end of method cross_prod (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : cross_prod (3)
// DESC   : compute idir component of cross/vector product 
//             C = A X B
//          
// ARGS   : A    : array of pointers to vectors; must each have at least 3
//                 elements for 3-d vector product. All vector elements must
//                 have the same length if non-NULL. May have NULL components, 
//                 which are treated as 0
//          B    : array of pointers to vectors; must each have at least 3
//                 elements for 3-d vector product. All vector elements must
//                 have the same length. All components must be non-NULL
//          iind : pointer to indirection array, used if non_NULLPTR.
//          nind : number of indices in iind. This is the number of cross products
//                 that will be computed, if iind is non-NULLPTR
//          idir : which component to compute
//          C    : idir component
// RETURNS: GTVector & 
//**********************************************************************************
template<typename T>
void cross_prod(GTVector<GTVector<T>*> &A, GTVector<GTVector<T>*> &B, 
                GINT idir, GTVector<T> &C)
{
	GEOFLOW_TRACE();
  assert( (A.size() == 2 || A.size() == 3) && "Incompatible dimensionality");
  assert( A.size() == B.size() && "Incompatible input vectors");

  T      x1, y1, z1;
  T      x2, y2, z2;

  if ( A.size() == 2 ) { // 2d vector
    if ( idir != 3 ) { // only a 3-component
      C = 0.0;
      return;
    }
    if ( A[0] != NULLPTR && A[1] != NULLPTR ) {
      for ( auto n=0; n<A[0]->size(); n++ ) { // cycle over all elemsnts
        x1 = (*A[0])[n]; y1 = (*A[1])[n]; 
        x2 = (*B[0])[n]; y2 = (*B[1])[n]; 
        C[n] = x1*y2 - x2*y1;
      }
    }
    else if ( A[0] != NULLPTR && A[1] == NULLPTR ) {
      for ( auto n=0; n<A[0]->size(); n++ ) { // cycle over all elemsnts
        x1 = (*A[0])[n]; 
        y2 = (*B[1])[n]; 
        C[n] = x1*y2;
      }
    }
    else if ( A[0] == NULLPTR && A[1] != NULLPTR ) {
      for ( auto n=0; n<A[1]->size(); n++ ) { // cycle over all elemsnts
        y1 = (*A[1])[n]; 
        x2 = (*B[0])[n]; 
        C[n] = -x2*y1;
      }
    }
    else {
      assert(FALSE);
    }
  }
  else { // 3d vectors
    switch ( idir ) {
    case 1:
      if ( A[1] != NULLPTR && A[2] != NULLPTR ) {
        for ( auto n=0; n<A[1]->size(); n++ ) { // cycle over all coord pairs
          y1 = (*A[1])[n]; z1 = (*A[2])[n];
          y2 = (*B[1])[n]; z2 = (*B[2])[n];
          C[n] = y1*z2 - z1*y2; 
        }
      }
      else if ( A[1] != NULLPTR && A[2] == NULLPTR ) {
        for ( auto n=0; n<A[0]->size(); n++ ) { // cycle over all coord pairs
          y1 = (*A[1])[n]; 
          z2 = (*B[2])[n];
          C[n] = y1*z2;
        }
      }
      else if ( A[1] == NULLPTR && A[2] != NULLPTR ) {
        for ( auto n=0; n<A[2]->size(); n++ ) { // cycle over all coord pairs
          z1 = (*A[2])[n];
          y2 = (*B[1])[n]; 
          C[n] =  -z1*y2; 
        }
      }
      else {
        assert(FALSE);
      }
      break;
    case 2:
      if ( A[0] != NULLPTR && A[2] != NULLPTR ) {
        for ( auto n=0; n<A[0]->size(); n++ ) { // cycle over all coord pairs
          x1 = (*A[0])[n]; z1 = (*A[2])[n];
          x2 = (*B[0])[n]; z2 = (*B[2])[n];
          C[n] = z1*x2 - z2*x1; 
        }
      }
      else if ( A[0] != NULLPTR && A[2] == NULLPTR ) {
        for ( auto n=0; n<A[0]->size(); n++ ) { // cycle over all coord pairs
          x1 = (*A[0])[n]; 
          z2 = (*B[2])[n]; 
          C[n] = -z2*x1; 
        }
      }
      else if ( A[0] == NULLPTR && A[2] != NULLPTR ) {
        for ( auto n=0; n<A[2]->size(); n++ ) { // cycle over all coord pairs
          z1 = (*A[2])[n];
          x2 = (*B[0])[n]; 
          C[n] = z1*x2;
        }
      }
      else {
        assert(FALSE);
      }
      break;
    case 3:
      if ( A[0] != NULLPTR && A[1] != NULLPTR ) {
        for ( auto n=0; n<A[0]->size(); n++ ) { // cycle over all coord pairs
          x1 = (*A[0])[n]; y1 = (*A[1])[n];
          x2 = (*B[0])[n]; y2 = (*B[1])[n];
          C[n] = x1*y2 - x2*y1;
        }
      }
      else if ( A[0] != NULLPTR && A[1] == NULLPTR ) {
        for ( auto n=0; n<A[0]->size(); n++ ) { // cycle over all coord pairs
          x1 = (*A[0])[n]; 
          y2 = (*B[1])[n];
          C[n] = x1*y2;
        }
      }
      else if ( A[0] == NULLPTR && A[1] != NULLPTR ) {
        for ( auto n=0; n<A[1]->size(); n++ ) { // cycle over all coord pairs
          y1 = (*A[1])[n];
          x2 = (*B[0])[n]; 
          C[n] = -x2*y1;
        }
      }
      else {
        assert(FALSE);
      }
      break;
    default:
      assert(FALSE);
    }
  }

} // end of method cross_prod (3)


//**********************************************************************************
//**********************************************************************************
// METHOD : cross_prod_s 
// DESC   : Compute idir component of cross/vector product 
//             C = A X B
//          where A's components are constants
//          
// ARGS   : A    : Vector of constants, of same length as B.
//          B    : array of pointers to vectors; must each have at least 3
//                 elements for 3-d vector product. All vector elements must
//                 have the same length. All components must be non-NULL
//          iind : pointer to indirection array, used if non_NULLPTR.
//          nind : number of indices in iind. This is the number of cross products
//                 that will be computed, if iind is non-NULLPTR
//          idir : which component to compute
//          C    : idir component
// RETURNS: GTVector & 
//**********************************************************************************
template<typename T>
void cross_prod_s(GTVector<T> &A, GTVector<GTVector<T>*> &B, 
                GINT idir, GTVector<T> &C)
{
	GEOFLOW_TRACE();
  assert( A.size() == B.size() && "Incompatible input vectors");
  assert(idir > 0 && idir < A.size());

  T      x1, y1, z1;
  T      x2, y2, z2;

  if ( A.size() == 2 ) { // 2d vectors
    if ( idir != 3 ) { // only a 3-component
      C = 0.0;
      return;
    }
    for ( auto n=0; n<B[0]->size(); n++ ) { // cycle over all elements
      x1 = A[0]; y1 = A[1]; 
      x2 = (*B[0])[n]; y2 = (*B[1])[n]; 
      C[n] = x1*y2 - x2*y1;
    }
  }
  else { // 3d vectors
    switch ( idir ) {
    case 1:
      for ( auto n=0; n<B[1]->size(); n++ ) { // cycle over all coord pairs
        y1 = A[1]      ; z1 = A[2];
        y2 = (*B[1])[n]; z2 = (*B[2])[n];
        C[n] = y1*z2 - z1*y2; 
      }
      break;
    case 2:
      for ( auto n=0; n<B[0]->size(); n++ ) { // cycle over all coord pairs
        x1 = A[0]      ; z1 = A[2];
        x2 = (*B[0])[n]; z2 = (*B[2])[n];
        C[n] = z1*x2 - z2*x1; 
      }
      break;
    case 3:
      for ( auto n=0; n<B[0]->size(); n++ ) { // cycle over all coord pairs
        x1 = A[0]      ; y1 = A[1];
        x2 = (*B[0])[n]; y2 = (*B[1])[n];
        C[n] = x1*y2 - x2*y1;
      }
      break;
    default:
      assert(FALSE);
    }
  }

} // end of method cross_prod_s (3)


//**********************************************************************************
//**********************************************************************************
// METHOD : normalize_euclidean
// DESC   : 
//             Compute Euclidean norm of each 'point', and return, overriding
//             input vectors
//          
// ARGS   : x    : array of pointers to vectors; must each have at least 3
//                 elements for 3-d vector product. All vector elements must
//                 have the same length
//          iind : pointer to indirection array, used if non_NULLPTR.
//          nind : number of indices in iind. This is the number of cross products
//                 that will be computed, if iind is non-NULLPTR
//          x0   : normalization constant
// RETURNS: GTVector & 
//**********************************************************************************
template<typename T>
void normalize_euclidean(GTVector<GTVector<T>*> &x, GINT *iind, GINT nind, T x0)
{
	GEOFLOW_TRACE();
  GSIZET n;
  T      xn;

  // NOTE: The better way to do this, especially for long lists is use
  //       outer loop over coord dimensions, and store progressive sum
  //       over each for each tuple. But this requires vector storage for
  //       each coordinate for all tuples.
  xn = 0.0;
  if ( iind != NULLPTR ) {
    for ( auto k=0; k<nind; k++ ) { // cycle over all n-tuples
      n = iind[k];
      xn = 0.0;
      for ( auto l=0; l<x.size(); l++ ) xn += (*x[l])[n]*(*x[l])[n];
      xn = x0/pow(xn,0.5);
      for ( auto l=0; l<x.size(); l++ ) (*x[l])[n] *= xn;
    }
  }
  else {
    for ( auto n=0; n<x[0]->size(); n++ ) { // cycle over all n-tuples
      xn = 0.0;
      for ( auto l=0; l<x.size(); l++ ) xn += (*x[l])[n]*(*x[l])[n];
      xn = x0/pow(xn,0.5);
      for ( auto l=0; l<x.size(); l++ ) (*x[l])[n] *= xn;
    }
  }

} // end of method normalize_euclidean


//**********************************************************************************
//**********************************************************************************
// METHOD :    paxy (1)
// DESC   : 
//             compute z = axy
//          
// ARGS   : z : return vector
//          x : vector factor
//          a : const multiplying x
//          y : vector factor
//          
// RETURNS: none
//**********************************************************************************
template<typename T>
void paxy(GTVector<T> &z, const GTVector<T> &x, T a, const GTVector<T> &y) 
{
	GEOFLOW_TRACE();
  if ( y.size() > 1 ) {
    for ( auto j=0; j<x.size(); j++ ) { 
      z[j] = a*x[j]*y[j];
    }
  }
  else { // to make consistent with GTVector
    for ( auto j=0; j<x.size(); j++ ) { 
      z[j] = a*x[j]*y[0];
    }
  }
} // end of method paxy (1)


//**********************************************************************************
//**********************************************************************************
// METHOD :    paxy (2)
// DESC   : 
//             compute x = axy
//          
// ARGS   : 
//          x : vector factor, returned
//          a : const multiplying x
//          y : vector factor
//          
// RETURNS: none
//**********************************************************************************
template<typename T>
void paxy(GTVector<T> &x, T a, const GTVector<T> &y) 
{
	GEOFLOW_TRACE();
  if ( y.size() > 1 ) {
    for ( auto j=0; j<x.size(); j++ ) { 
      x[j] = a*x[j]*y[j];
    }
  }
  else { // to make consistent with GTVector
    for ( auto j=0; j<x.size(); j++ ) { 
      x[j] = a*x[j]*y[0];
    }
  }
} // end of method paxy (2)


//**********************************************************************************
//**********************************************************************************
// METHOD :    saxpby (1)
// DESC   : 
//             compute x = ax + by
//          
// ARGS   : x : vector , updated
//          a : const multiplying x
//          y : vector, must be same size as x
//          b : const multiplying y
//          
// RETURNS: GTVector & 
//**********************************************************************************
template<typename T>
void saxpby(GTVector<T> &x, T a, GTVector<T> &y, T b) 
{
	GEOFLOW_TRACE();
  if ( y.size() > 1 ) {
    for ( auto j=0; j<x.size(); j++ ) { 
      x[j] = a*x[j] + b*y[j];
    }
  }
  else { // to make consistent with GTVector
    for ( auto j=0; j<x.size(); j++ ) { 
      x[j] = a*x[j] + b*y[0];
    }
  }
} // end of method saxpby (1)


//**********************************************************************************
//**********************************************************************************
// METHOD :    saxpby (2)
// DESC   : 
//             compute z = ax + by
//          
// ARGS   : z : vector, updated
//          x : summand vector
//          a : const multiplying x
//          y : vector, must be same size as x
//          b : const multiplying y
//          
// RETURNS: GTVector & 
//**********************************************************************************
template<typename T>
void saxpby(GTVector<T> &z, GTVector<T> &x, T a, GTVector<T> &y, T b) 
{
	GEOFLOW_TRACE();
  if ( y.size() > 1 ) {
    for ( auto j=0; j<x.size(); j++ ) { 
      z[j] = a*x[j] + b*y[j];
    }
  }
  else { // to make consistent with GTVector
    for ( auto j=0; j<x.size(); j++ ) { 
      z[j] = a*x[j] + b*y[0];
    }
  }
} // end of method saxpby (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : normalizeL2
// DESC   :
//             L2 Normalize input field, u --> c u,  s.t.
//               0.5 * Int (c u)^2  dV / Int dV = E0
//
// ARGS   :
//          grid : grid object
//          u    : array of pointers to vectors; must each have at least 3
//                 elements for 3-d vector product. All vector elements must
//                 have the same length. Is normalized on exit.
//          tmp  : tmp vector of length at least 2, each
//                 of same length as x
//          E0   : normalization value
// RETURNS: none.
//**********************************************************************************
template<typename Grid, typename T>
void normalizeL2(Grid &grid, GTVector<GTVector<T>*> &u, GTVector<GTVector<T>*> &tmp, T E0)
{
	GEOFLOW_TRACE();
  GSIZET  n;
  GDOUBLE xn, xint;

  // xinit gives 0.5  Int _u_^2 dV / Int dV:
  xint = static_cast<GDOUBLE>(GMTK::energy<Grid,T>(grid, u, tmp, TRUE, FALSE));
  xn   = sqrt(E0 / xint);
  for ( GINT l=0; l<u.size(); l++ ) {
    *u[l] *= xn;
  }


} // end of method normalizeL2


//**********************************************************************************
//**********************************************************************************
// METHOD : zero (1)
// DESC   :
//          Set values < std::numeric_limits<>::epsilon()
//          identically to 0
//
// ARGS   :
//          u    : field
// RETURNS: none.
//**********************************************************************************
template<typename T>
void zero(GTVector<T> &u)
{
	GEOFLOW_TRACE();
  T tiny=std::numeric_limits<T>::epsilon();

  for ( auto j=0; j<u.size(); j++ ) {
    if ( fabs(u[j]) < tiny ) u[j] = 0.0;
  }


} // end of method zero (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : zero (2)
// DESC   :
//          Set values < std::numeric_limits<>::epsilon()
//          identically to 0
//
// ARGS   :
//          u    : field
// RETURNS: none.
//**********************************************************************************
template<typename T>
void zero(GTVector<GTVector<T>*> &u)
{
	GEOFLOW_TRACE();
  T tiny=std::numeric_limits<T>::epsilon();

  for ( auto i=0; i<u.size(); i++ ) {
    for ( auto j=0; j<u[i]->size(); j++ ) {
      if ( fabs((*u[i])[j]) < tiny ) (*u[i])[j] = 0.0;
    }
  }


} // end of method zero (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : zero (3)
// DESC   :
//          Set values < std::numeric_limits<>::epsilon()
//          identically to 0
//
// ARGS   :
//          u    : field
// RETURNS: none.
//**********************************************************************************
template<typename T>
void zero(GTVector<GTVector<T>> &u)
{
	GEOFLOW_TRACE();
  T tiny=std::numeric_limits<T>::epsilon();

  for ( auto i=0; i<u.size(); i++ ) {
    for ( auto j=0; j<u[i].size(); j++ ) {
      if ( fabs(u[i][j]) < tiny ) u[i][j] = 0.0;
    }
  }


} // end of method zero (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : maxbyelem
// DESC   : Find max of q on each element, return in max array
//
// ARGS   : grid: Grid object
//          q   : field over all elems.
//          max : max of q on each element. Allocated if necessary.
// RETURNS: none.
//**********************************************************************************
template<typename Grid, typename T>
void maxbyelem(Grid &grid, GTVector<T> &q, GTVector<T> &max)
{
	GEOFLOW_TRACE();
  GSIZET     ibeg, iend;
  typename Grid::GElemList *elems = &grid.elems();

  max.resizem(grid.nelems());

  for ( auto e=0; e<grid.nelems(); e++ ) {
    ibeg = (*elems)[e]->igbeg(); iend = (*elems)[e]->igend();
    q.range(ibeg, iend); // restrict global vecs to local range
    max[e] = q.amax();
  } // end, elem loop
  q.range_reset();

} // end of method maxbyelem


//**********************************************************************************
//**********************************************************************************
// METHOD : minbyelem
// DESC   : Find min of q on each element, return in min array
//
// ARGS   : grid: Grid object
//          q   : field over all elems.
//          min : min of q on each element. Allocated if necessary.
// RETURNS: none.
//**********************************************************************************
template<typename Grid, typename T>
void minbyelem(Grid &grid, GTVector<T> &q, GTVector<T> &min)
{
	GEOFLOW_TRACE();
  GSIZET ibeg, iend;
  typename Grid::GElemList *elems = &grid.elems();

  min.resizem(grid.nelems());

  for ( auto e=0; e<grid.nelems(); e++ ) {
    ibeg = (*elems)[e]->igbeg(); iend = (*elems)[e]->igend();
    q.range(ibeg, iend); // restrict global vecs to local range
    min[e] = q.amin();
  } // end, elem loop
  q.range_reset();

} // end of method maxbyelem


//**********************************************************************************
//**********************************************************************************
// METHOD : I2_X_D1 (1)
// DESC   : Apply tensor product operator to vector:
//            y = I2 X D1 u
//          where I2 is the 2-direction's identity, and D1 is
//          the deriv matrix in the 1-direction
// ARGS   : D1  : 1-direction (dense) operator
//          u   : operand vector; of size N1 X N2
//          N1-2: dimensions of u if interpreted as a matrix
//          y   : return vector result; must be at least of size
//                N1 X N2
// RETURNS: none
//**********************************************************************************
template<typename T>
void I2_X_D1(GTMatrix<T> &D1,
             GTVector<T> &u, GSIZET N1, GSIZET N2, GTVector<T> &y)
{
  GEOFLOW_TRACE_RENAME("I2_X_D1(5 args)");
  ASSERT_MSG((u.size() >= N1*N2 && y.size() >= N1*N2), "GMTK::I2_X_D1 (1) incompatible size");
  
  GSIZET ND1, ND2;
  ND1 = D1.size(1);
  ND2 = D1.size(2);

  // Compute y = I2_X_D1 u:
  if      ( std::is_same<T,GFLOAT>::value ) {
    fmxm((GFLOAT*)y.data(), (GFLOAT*)(D1.data().data()), &ND1, &ND2, (GFLOAT*)(u.data()), &N1, &N2, &szMatCache_);
  }
  else if ( std::is_same<T,GDOUBLE>::value ) {
    dmxm((GDOUBLE*)y.data(), (GDOUBLE*)(D1.data().data()), &ND1, &ND2, (GDOUBLE*)(u.data()), &N1, &N2, &szMatCache_);
  }
  else if ( std::is_same<T,GQUAD>::value ) {
    qmxm((GQUAD*)y.data(), (GQUAD*)(D1.data().data()), &ND1, &ND2, (GQUAD*)(u.data()), &N1, &N2, &szMatCache_);
  }
  else {
    assert(FALSE);
  }

} // end of method I2_X_D1 (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : I2_X_D1  (2)
// DESC   : Apply tensor product operator to vector:
//            y = I2 X D1 u
//          where I2 is the 2-direction's identity, and D1 is
//          the deriv matrix in the 1-direction, and u is the
//          the full state component
// ARGS   : D1  : 1-direction (dense) operator
//          u   : operand vector consisting of Ne 'elements' 
//                each of size ize N1 X N2
//          N1-2: dimensions of u elemewnts if interpreted as matrices
//          Ne  : number of 'elements' in u
//          y   : return vector result; must be at least of size
//                N1 X N2
// RETURNS: none
//**********************************************************************************
template<typename T>
void I2_X_D1(GTMatrix<T> &D1,
             GTVector<T> &u, GSIZET N1, GSIZET N2, GSIZET Ne, GTVector<T> &y)
{
  GEOFLOW_TRACE_RENAME("I2_X_D1(6 args)");
  ASSERT_MSG((u.size() >= N1*N2 && y.size() >= N1*N2), "GMTK::I2_X_D1 (2) incompatible size");
  GSIZET ND1, ND2, Nu;

  ND1 = D1.size(1);
  ND2 = D1.size(2);
  Nu = N2 * Ne;

  // Compute y = I2_X_D1 u:
  if      ( std::is_same<T,GFLOAT>::value ) {
    fmxm((GFLOAT*)y.data(), (GFLOAT*)(D1.data().data()), &ND1, &ND2, (GFLOAT*)u.data(), &N1, &Nu, &szMatCache_);
  }
  else if ( std::is_same<T,GDOUBLE>::value ) {
    dmxm((GDOUBLE*)y.data(), (GDOUBLE*)(D1.data().data()), &ND1, &ND2, (GDOUBLE*)u.data(), &N1, &Nu, &szMatCache_);
  }
  else if ( std::is_same<T,GQUAD>::value ) {
    qmxm((GQUAD*)y.data(), (GQUAD*)(D1.data().data()), &ND1, &ND2, (GQUAD*)u.data(), &N1, &Nu, &szMatCache_);
  }
  else {
    assert(FALSE);
  }

} // end of method I2_X_D1 (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : I2_X_D1  (3)
// DESC   : Apply tensor product operator to vector:
//            y = I2 X D1 u
//          where I2 is the 2-direction's identity, and D1 is
//          the deriv matrix in the 1-direction, and u is the
//          the full state component
// ARGS   : D1   : 1-direction (dense) operator
//          u    : operand vector consisting of Ne 'elements' 
//                 each of size ize N1 X N2
//          N1-2 : dimensions of u elemewnts if interpreted as matrices
//          Ne   : number of 'elements' in u
//          cudat: cuMatBlockDat structure data
//          y    : return vector result; must be at least of size
//                 N1 X N2
// RETURNS: none
//**********************************************************************************
template<typename T>
void I2_X_D1(GTMatrix<T> &D1,
             GTVector<T> &u, GSIZET N1, GSIZET N2, GSIZET Ne, GCBLAS::cuMatBlockDat &cudat, GTVector<T> &y)
{
  GEOFLOW_TRACE_RENAME("I2_X_D1(7 args)");
  ASSERT_MSG((u.size() >= N1*N2 && y.size() >= N1*N2), "GMTK::I2_X_D1 (3) incompatible size");

  GSIZET ND1, ND2, Nu;
  GINT   M, N, K, lda, ldb, ldc;

  ND1 = D1.size(1);
  ND2 = D1.size(2);
  Nu = N2 * Ne;

#if defined(GEOFLOW_USE_CBLAS) || defined(GEOFLOW_USE_CUBLAS)

  M   = ND1;
  N   = N2*Ne;
  K   = ND2;
  lda = M;
  ldb = K;
  ldc = M;
  GCBLAS::gemm<T>(cudat.hcublas, GCBLAS::CblasColMajor, GCBLAS::CblasNoTrans, GCBLAS::CblasNoTrans,
                    M, N, K, 1.0, (T*)(D1.data().data()), lda, (T*)(u.data()), ldb, 0.0, (T*)y.data(), ldc);
#else

  if      ( std::is_same<T,GFLOAT>::value ) {
    fmxm((GFLOAT*)y.data(), (GFLOAT*)(D1.data().data()), &ND1, &ND2, (GFLOAT*)u.data(), &N1, &Nu, &szMatCache_);
  }
  else if ( std::is_same<T,GDOUBLE>::value ) {
    dmxm((GDOUBLE*)y.data(), (GDOUBLE*)(D1.data().data()), &ND1, &ND2, (GDOUBLE*)u.data(), &N1, &Nu, &szMatCache_);
  }
  else if ( std::is_same<T,GQUAD>::value ) {
    qmxm((GQUAD*)y.data(), (GQUAD*)(D1.data().data()), &ND1, &ND2, (GQUAD*)u.data(), &N1, &Nu, &szMatCache_);
  }
  else {
    assert(FALSE);
  }

#endif

} // end of method I2_X_D1 (3)


//**********************************************************************************
//**********************************************************************************
// METHOD : D2_X_D1 
// DESC   : Apply tensor product operator to vector:
//            y = D2 X D1 u
// ARGS   : D1  : 1-direction (dense) operator 
//          D2T : transpose of 2-direction (dense) operator
//          u   : operand vector; must be at least of size
//                D1.size(2) x D2.size(2) = D1.size(2) x D2T.size(1)
//          tmp : temp space; may be re-allocated, but only if 
//                current size < required (must be at least 
//                D1.size(1) x D2.size(2) = D1.size(1) x D2T.size(1)
//                before reallocation is done).
//          y   : return vector result; must be at least of size
//                D1.size(1) x D2.size(1) = D1.size(1) x D2T.size(2)
// RETURNS: none
//**********************************************************************************
template <typename T>
void D2_X_D1(GTMatrix<T> &D1, GTMatrix<T>  &D2T, 
             GTVector<T> &u, GTVector<T> &tmp, GTVector<T> &y)
{
	GEOFLOW_TRACE();
  
  GSIZET   N11, N12, N21, N22;
  N11 = D1 .size(1);
  N12 = D1 .size(2);
  N21 = D2T.size(1);
  N22 = D2T.size(2);
  ASSERT_MSG((u.size() >= N12*N21 && y.size() >= N11*N22), "GMTK::D2_X_D1 incompatible size");

  // Compute y = D2_X_D1 u as: y = D1 U D2T, where U is u 
  // considered in matrix form:

  // Resize tmp only if its current size is less than required:
  tmp.resizem(N11*N21);

  if      ( std::is_same<T,GFLOAT>::value ) {
    // tmp = I2_X_D1 * u == D1 U (in mat form):
    fmxm((GFLOAT*)tmp.data(), (GFLOAT*)D1.data().data(), &N11, &N12, (GFLOAT*)u.data(), &N12, &N21, &szMatCache_);

    // y = D2_X_I1 * tmp == TMP D2T (in mat form):
    fmxm((GFLOAT*)y.data(), (GFLOAT*)tmp.data(), &N11, &N12, (GFLOAT*)D2T.data().data(), &N21, &N22, &szMatCache_);
  }
  else if ( std::is_same<T,GDOUBLE>::value ) {
    // tmp = I2_X_D1 * u == D1 U (in mat form):
    dmxm((GDOUBLE*)tmp.data(), (GDOUBLE*)D1.data().data(), &N11, &N12, (GDOUBLE*)u.data(), &N12, &N21, &szMatCache_);

    // y = D2_X_I1 * tmp == TMP D2T (in mat form):
    dmxm((GDOUBLE*)y.data(), (GDOUBLE*)tmp.data(), &N11, &N12, (GDOUBLE*)D2T.data().data(), &N21, &N22, &szMatCache_);

  }
  else if ( std::is_same<T,GQUAD>::value ) {
    // tmp = I2_X_D1 * u == D1 U (in mat form):
    qmxm((GQUAD*)tmp.data(), (GQUAD*)D1.data().data(), &N11, &N12, (GQUAD*)u.data(), &N12, &N21, &szMatCache_);

    // y = D2_X_I1 * tmp == TMP D2T (in mat form):
    qmxm((GQUAD*)y.data(), (GQUAD*)tmp.data(), &N11, &N12, (GQUAD*)D2T.data().data(), &N21, &N22, &szMatCache_);

  }
  else {
    assert(FALSE);
  }


} // end of method D2_X_D1


//**********************************************************************************
//**********************************************************************************
// METHOD : Dg2_X_D1
// DESC   : Apply tensor product operator to vector:
//            y = diag(D2) X D1 u
//          where D2 is specified as a vector, and D1 is 
//          the deriv matrix in the 1-direction
// ARGS   : D1  : 1-direction (dense) operator 
//          Dg2 : diag(D2), as a vector
//          u   : operand vector; must be at least of size
//                D1.size(2) x D2.size(2) = D1.size(2) x D2T.size(1)
//          tmp : temp space; may be re-allocated, but only if 
//                current size < required (must be at least 
//                D1.size(1) x D2.size(2) = D1.size(1) x D2T.size(1)
//                before reallocation is done).
//          y   : return vector result; must be at least of size
//                D1.size(1) x D2.size(1) = D1.size(1) x D2T.size(2)
// RETURNS: none
//**********************************************************************************
template<typename T> 
void Dg2_X_D1(GTMatrix<T> &D1, GTVector<T> &Dg2, GTVector<T> &u, 
              GTVector<T> &tmp, GTVector<T> &y)
{
	GEOFLOW_TRACE();
  
  GSIZET   N11, N12, N2;
  N11 = D1.size(1);
  N12 = D1.size(2);
  N2  = Dg2.size();
  ASSERT_MSG((u.size() >= N11*N2 && y.size() >= N11*N2), "GMTK::Dg_X_D1 incompatible size");
 

  // Resize tmp only if its current size is less than required:
  tmp.resizem(N11*N2);

  // Compute y = Dg2_X_D1 u as (Dg2 X I) (I X D1) U:

  if      ( std::is_same<T,GFLOAT>::value ) {
    // tmp = I2_X_D1 * u == D1 U (in mat form):
    fmxm((GFLOAT*)tmp.data(), (GFLOAT*)D1.data().data(), &N11, &N12, (GFLOAT*)u.data(), &N12, &N2, &szMatCache_);

    // y = Dg2_X_I1 * tmp == TMP diag(D2T) = TMP diag(D2)  (in mat form):
    fmxDm((GFLOAT*)y.data(), (GFLOAT*)tmp.data(), &N11, &N2, (GFLOAT*)Dg2.data(), &N2, &szMatCache_);
  }
  else if ( std::is_same<T,GDOUBLE>::value ) {
    // tmp = I2_X_D1 * u == D1 U (in mat form):
    dmxm((GDOUBLE*)tmp.data(), (GDOUBLE*)D1.data().data(), &N11, &N12, (GDOUBLE*)u.data(), &N12, &N2, &szMatCache_);

    // y = Dg2_X_I1 * tmp == TMP diag(D2T) = TMP diag(D2)  (in mat form):
    dmxDm((GDOUBLE*)y.data(), (GDOUBLE*)tmp.data(), &N11, &N2, (GDOUBLE*)Dg2.data(), &N2, &szMatCache_);
  }
  else if ( std::is_same<T,GQUAD>::value ) {
    // tmp = I2_X_D1 * u == D1 U (in mat form):
    qmxm((GQUAD*)tmp.data(), (GQUAD*)D1.data().data(), &N11, &N12, (GQUAD*)u.data(), &N12, &N2, &szMatCache_);

    // y = Dg2_X_I1 * tmp == TMP diag(D2T) = TMP diag(D2)  (in mat form):
    qmxDm((GQUAD*)y.data(), (GQUAD*)tmp.data(), &N11, &N2, (GQUAD*)Dg2.data(), &N2, &szMatCache_);
  }
  else {
    assert(FALSE);
  }


} // end of method Dg2_X_D1


//**********************************************************************************
//**********************************************************************************
// METHOD : D2_X_I1  (1)
// DESC   : Apply tensor product operator to vector:
//            y = D2 X I1 u
//          where I1 is the 1-direction's identity, and D2
//          the deriv matrix in the 2-direction
// ARGS   : D2T : 2-direction (dense) operator transpose 
//          u   : operand vector; must be at least of size
//                D1.size(2) x D2.size(2) = D1.size(2) x D2T.size(1)
//          N1-2: dimensions of u if interpreted as matrix
//          y   : return vector result; must be at least of size
//                D1.size(1) x D2.size(1) = D1.size(1) x D2T.size(2)
// RETURNS: none
//**********************************************************************************
template<typename T>
void D2_X_I1(GTMatrix<T> &D2T, 
              GTVector<T> &u, GSIZET N1, GSIZET N2, GTVector<T> &y)
{
	GEOFLOW_TRACE();
  ASSERT_MSG((u.size() >= N1*N2 && y.size() >= N1*N2), "GMTK::D2_X_I1 (1) incompatible size");

  GSIZET N21, N22;
  N21 = D2T.size(1);
  N22 = D2T.size(2);

  // Compute y = I2_X_D1 u = u * D2T:
  if      ( std::is_same<T,GFLOAT>::value ) {
    fmxm((GFLOAT*)y.data(), (GFLOAT*)u.data(), &N1, &N2, (GFLOAT*)D2T.data().data(), &N21, &N22, &szMatCache_);
  }
  else if ( std::is_same<T,GDOUBLE>::value ) {
    dmxm((GDOUBLE*)y.data(), (GDOUBLE*)u.data(), &N1, &N2, (GDOUBLE*)D2T.data().data(), &N21, &N22, &szMatCache_);
  }
  else if ( std::is_same<T,GQUAD>::value ) {
    qmxm((GQUAD*)y.data(), (GQUAD*)u.data(), &N1, &N2, (GQUAD*)D2T.data().data(), &N21, &N22, &szMatCache_);
  }
  else {
    assert(FALSE);
  }

} // end of method D2_X_I1 (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : D2_X_I1 (2)
// DESC   : Apply tensor product operator to vector:
//            y = D2 X I1 u
//          where I1 is the 1-direction's identity, and D2
//          the deriv matrix in the 2-direction
// ARGS   : D2T : 2-direction (dense) operator transpose 
//          u   : operand vector consisting of Ne 'elements' 
//                each of size ize N1 X N2
//          N1-2: dimensions of u elemewnts if interpreted as matrices
//          Ne  : number of 'elements' in u
//          y   : return vector result; must be at least of size
//                D1.size(1) x D2.size(1) = D1.size(1) x D2T.size(2)
// RETURNS: none
//**********************************************************************************
template<typename T>
void D2_X_I1(GTMatrix<T> &D2T, 
              GTVector<T> &u, GSIZET N1, GSIZET N2, GSIZET Ne, GTVector<T> &y)
{
	GEOFLOW_TRACE();
  GSIZET N21, N22, Nu;

  N21 = D2T.size(1);
  N22 = D2T.size(2);
  Nu = N1 * N2;
  ASSERT_MSG((u.size() >= N1*N2*Ne && y.size() >= N1*N2*Ne), "GMTK::D2_X_I1 (2) incompatible size");

  // Compute y = I2_X_D1 u = u * D2T:
  if      ( std::is_same<T,GFLOAT>::value ) {
    for ( auto i=0; i<Ne; i++ ) {
      fmxm((GFLOAT*)y.data()+i*Nu, (GFLOAT*)u.data()+i*Nu, &N1, &Nu, (GFLOAT*)D2T.data().data(), &N21, &N22, &szMatCache_);
    }
  }
  else if ( std::is_same<T,GDOUBLE>::value ) {
    for ( auto i=0; i<Ne; i++ ) {
      dmxm((GDOUBLE*)y.data()+i*Nu, (GDOUBLE*)u.data()+i*Nu, &N1, &Nu, (GDOUBLE*)D2T.data().data(), &N21, &N22, &szMatCache_);
    }
  }
  else if ( std::is_same<T,GQUAD>::value ) {
    for ( auto i=0; i<Ne; i++ ) {
      qmxm((GQUAD*)y.data()+i*Nu, (GQUAD*)u.data()+i*Nu, &N1, &Nu, (GQUAD*)D2T.data().data(), &N21, &N22, &szMatCache_);
    }
  }
  else {
    assert(FALSE);
  }

} // end of method D2_X_I1 (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : D2_X_I1 (3)
// DESC   : Apply tensor product operator to vector:
//            y = D2 X I1 u
//          where I1 is the 1-direction's identity, and D2
//          the deriv matrix in the 2-direction
// ARGS   : D2T  : 2-direction (dense) operator transpose 
//          u    : operand vector consisting of Ne 'elements' 
//                 each of size ize N1 X N2
//          N1-2 : dimensions of u elemewnts if interpreted as matrices
//          Ne   : number of 'elements' in u
//          cudat: cuMatBlockDat structure data
//          y    : return vector result; must be at least of size
//                  D1.size(1) x D2.size(1) = D1.size(1) x D2T.size(2)
// RETURNS: none
//**********************************************************************************
template<typename T>
void D2_X_I1(GTMatrix<T> &D2T, 
              GTVector<T> &u, GSIZET N1, GSIZET N2, GSIZET Ne, GCBLAS::cuMatBlockDat &cudat, GTVector<T> &y)
{
	GEOFLOW_TRACE();
  GSIZET N21, N22, Nu;
  GINT   M, N, K, lda, ldb, ldc;


  N21 = D2T.size(1);
  N22 = D2T.size(2);
  Nu = N1 * N2;
  ASSERT_MSG((u.size() >= N1*N2*Ne && y.size() >= N1*N2*Ne), "GMTK::D2_X_I1 (2) incompatible size");

// Compute y = I2_X_D1 u = u * D2T:
#if defined(GEOFLOW_USE_CBLAS) || defined(GEOFLOW_USE_CUBLAS)
  M   = N21;
  N   = N2;
  K   = N22;
  lda = M;
  ldb = K;
  ldc = M;
  
  GCBLAS::batched_gemm<T>(cudat, GCBLAS::CblasColMajor, GCBLAS::CblasNoTrans, GCBLAS::CblasNoTrans,
                           M, N, K, 1.0, (T*)(u.data()), lda, (T*)(D2T.data().data()), ldb, 0.0, (T*)y.data(), ldc);
#else

  if      ( std::is_same<T,GFLOAT>::value ) {
    for ( auto i=0; i<Ne; i++ ) {
      fmxm((GFLOAT*)y.data()+i*Nu, (GFLOAT*)u.data()+i*Nu, &N1, &Nu, (GFLOAT*)D2T.data().data(), &N21, &N22, &szMatCache_);
    }
  }
  else if ( std::is_same<T,GDOUBLE>::value ) {
    for ( auto i=0; i<Ne; i++ ) {
      dmxm((GDOUBLE*)y.data()+i*Nu, (GDOUBLE*)u.data()+i*Nu, &N1, &Nu, (GDOUBLE*)D2T.data().data(), &N21, &N22, &szMatCache_);
    }
  }
  else if ( std::is_same<T,GQUAD>::value ) {
    for ( auto i=0; i<Ne; i++ ) {
      qmxm((GQUAD*)y.data()+i*Nu, (GQUAD*)u.data()+i*Nu, &N1, &Nu, (GQUAD*)D2T.data().data(), &N21, &N22, &szMatCache_);
    }
  }
  else {
    assert(FALSE);
  }

#endif

} // end of method D2_X_I1 (3)


//**********************************************************************************
//**********************************************************************************
// METHOD : D2_X_Dg1 
// DESC   : Apply tensor product operator to vector:
//            y = D2 X diag(D1) u
//          where Dg  is the 1-direction's operator expressed as a vector, and D2
//          the dense matrix in the 2-direction
// ARGS   : Dg1 : diag(D1), as a vector
//          D2T : 2-direction (dense) operator transpose 
//          u   : operand vector; must be at least of size
//                D1.size(2) x D2.size(2) = D1.size(2) x D2T.size(1)
//          tmp : temp space; may be re-allocated, but only if 
//                current size < required (must be at least 
//                D1.size(1) x D2.size(2) = D1.size(1) x D2T.size(1)
//                before reallocation is done).
//          y   : return vector result; must be at least of size
//                D1.size(1) x D2.size(1) = D1.size(1) x D2T.size(2)
// RETURNS: none
//**********************************************************************************
template<typename T>
void D2_X_Dg1(GTVector<T> &Dg1, GTMatrix<T> &D2T, GTVector<T> &u, 
              GTVector<T> &tmp, GTVector<T> &y)
{
	GEOFLOW_TRACE();
  GSIZET   N1, N21, N22;

  N1  = Dg1.size();
  N21 = D2T.size(1);
  N22 = D2T.size(2);
  ASSERT_MSG((u.size() >= N1*N21 && y.size() >= N1*N21), "GMTK::D2_X_Dg1 incompatible size");

  // Compute y = D2_X_Dg u as: y = (D2_X_I)(I_X_Dg) U,
  //  where U is u considered in matrix form:
  tmp.resizem(N1*N21);
  if      ( std::is_same<T,GFLOAT>::value ) {
    // y = D2_X_I1 * tmp == TMP D2T (in mat form):
    fmxm((GFLOAT*)tmp.data(), (GFLOAT*)u.data(), &N1, &N21, (GFLOAT*)D2T.data().data(), &N21, &N22, &szMatCache_);

    // tmp = I2_X_D1 * u == D1 U (in mat form):
    fDmxm((GFLOAT*)y.data(), (GFLOAT*)Dg1.data(), &N1, (GFLOAT*)tmp.data(), &N1, &N21, &szMatCache_);
  }
  else if ( std::is_same<T,GDOUBLE>::value ) {
    // y = D2_X_I1 * tmp == TMP D2T (in mat form):
    dmxm((GDOUBLE*)tmp.data(), (GDOUBLE*)u.data(), &N1, &N21, (GDOUBLE*)D2T.data().data(), &N21, &N22, &szMatCache_);

    // tmp = I2_X_D1 * u == D1 U (in mat form):
    dDmxm((GDOUBLE*)y.data(), (GDOUBLE*)Dg1.data(), &N1, (GDOUBLE*)tmp.data(), &N1, &N21, &szMatCache_);

  }
  else if ( std::is_same<T,GQUAD>::value ) {
    // y = D2_X_I1 * tmp == TMP D2T (in mat form):
    qmxm((GQUAD*)tmp.data(), (GQUAD*)u.data(), &N1, &N21, (GQUAD*)D2T.data().data(), &N21, &N22, &szMatCache_);

    // tmp = I2_X_D1 * u == D1 U (in mat form):
    qDmxm((GQUAD*)y.data(), (GQUAD*)Dg1.data(), &N1, (GQUAD*)tmp.data(), &N1, &N21, &szMatCache_);
  }
  else {
    assert(FALSE);
  }


} // end of method D2_X_Dg1


//**********************************************************************************
//**********************************************************************************
// METHOD : D3_X_D2_X_D1 
// DESC   : Apply tensor product operator to vector:
//            y = D3 X D2 X D1 u
// ARGS   : D1  : 1-direction (dense) operator 
//          D2T : transpose of 2-direction (dense) operator
//          D3T : transpose of 3-direction (dense) operator
//          u   : operand vector; must be at least of size
//                D1.size(2) x D2.size(2) x D3.dm(2) = D1.size(2) x D2T.size(1) x D3T.size(1)
//          tmp : temp space; may be re-allocated, but only if 
//                current size < required (must be at least 
//                D1.size(1) x D2.size(2)  = D1.size(1) x D2T.size(1)
//                before reallocation is done).
//          y   : return vector result; must be at least of size
//                D1.size(1) x D2.size(1) x D3.size(1) = D1.size(1) x D2T.size(2) x D3T.size(2)
// RETURNS: none
//**********************************************************************************
template<typename T>
void D3_X_D2_X_D1(GTMatrix<T> &D1, GTMatrix<T>  &D2T, GTMatrix<T> &D3T,
                  GTVector<T> &u, GTVector<T> &tmp, GTVector<T> &y)
{
	GEOFLOW_TRACE();
  GSIZET   N11, N12, N21, N22, N31, N32;

  N11 = D1 .size(1);
  N12 = D1 .size(2);
  N21 = D2T.size(1);
  N22 = D2T.size(2);
  N31 = D3T.size(1);
  N32 = D3T.size(2);
  ASSERT_MSG((u.size() >= N12*N21*N31 && y.size() >= N11*N22*N32), "GMTK::D2_X_D1 incompatible size");

  // Compute y = D3_X_D2_X_D1 u as (D3 X I2 X I1)(I3 X D2 X I1)(I3 X I2 X D1) U: 

  GSIZET nxy = N21*N31;

  // Resize tmp only if its current size is less than required:
  tmp.resizem(nxy*N32);

  if      ( std::is_same<T,GFLOAT>::value ) {
    fmxm((GFLOAT*)y.data(), (GFLOAT*)D1.data().data(), &N11, &N12, (GFLOAT*)u.data(), &N12, &nxy, &szMatCache_);

  // tmp = I3_X_D2_X_I1 y:
    for ( auto k=0; k<N32; k++ ) { // do mxm op for each 'plane':
      fmxm((GFLOAT*)(tmp.data()+k*N11*N22), (GFLOAT*)(y.data()+k*N11*N22), &N11, &N22, (GFLOAT*)D2T.data().data(), &N21, &N22, &szMatCache_);
    }

    // y = D3 X I X I tmp:
    nxy = N11*N22;
    fmxm((GFLOAT*)y.data(), (GFLOAT*)tmp.data(), &nxy, &N32, (GFLOAT*)D3T.data().data(), &N31, &N32, &szMatCache_);
  }
  else if ( std::is_same<T,GDOUBLE>::value ) {
    dmxm((GDOUBLE*)y.data(), (GDOUBLE*)D1.data().data(), &N11, &N12, u.data(), &N12, &nxy, &szMatCache_);

  // tmp = I3_X_D2_X_I1 y:
    for ( auto k=0; k<N32; k++ ) { // do mxm op for each 'plane':
      dmxm((GDOUBLE*)(tmp.data()+k*N11*N22), (GDOUBLE*)(y.data()+k*N11*N22), &N11, &N22, (GDOUBLE*)D2T.data().data(), &N21, &N22, &szMatCache_);
    }

    // y = D3 X I X I tmp:
    nxy = N11*N22;
    dmxm((GDOUBLE*)y.data(), (GDOUBLE*)tmp.data(), &nxy, &N32, (GDOUBLE*)D3T.data().data(), &N31, &N32, &szMatCache_);
  }
  else if ( std::is_same<T,GQUAD>::value ) {
    qmxm((GQUAD*)y.data(), (GQUAD*)D1.data().data(), &N11, &N12, (GQUAD*)u.data(), &N12, &nxy, &szMatCache_);

  // tmp = I3_X_D2_X_I1 y:
    for ( auto k=0; k<N32; k++ ) { // do mxm op for each 'plane':
      qmxm((GQUAD*)(tmp.data()+k*N11*N22), (GQUAD*)(y.data()+k*N11*N22), &N11, &N22, (GQUAD*)D2T.data().data(), &N21, &N22, &szMatCache_);
    }

    // y = D3 X I X I tmp:
    nxy = N11*N22;
    qmxm((GQUAD*)y.data(), (GQUAD*)tmp.data(), &nxy, &N32, (GQUAD*)D3T.data().data(), &N31, &N32, &szMatCache_);
  }
  else {
    assert(FALSE);
  }


} // end of method D3_X_D2_X_D1



//**********************************************************************************
//**********************************************************************************
// METHOD : I3_X_I2_X_D1 (1)
// DESC   : Apply tensor product operator to vector:
//            y = I3 X I2 X D1 u
// ARGS   : D1      : 1-direction (dense) operator 
//          u       : operand vector; must be at least of size
//                    N1*N2*N3, with 1 changing most rapidly, then, 2, then 3
//          N1,N2,N3: coord dimensions of u, y
//          y   : return vector result of size >= N1 * N2 * N3
// RETURNS: none
//**********************************************************************************
template<typename T>
void I3_X_I2_X_D1(GTMatrix<T> &D1, GTVector<T> &u,
                  GSIZET N1, GSIZET N2, GSIZET N3,
                  GTVector<T> &y)
{
	GEOFLOW_TRACE();
  GSIZET  ND1, ND2, NYZ, NN;

  ND1 = D1.size(1);
  ND2 = D1.size(2);
  NYZ = N2*N3;
  NN  = N1*N2*N3;
  ASSERT_MSG(!(u.size() < NN || y.size() < NN), "GMTK::I3_X_I2_X_D1 (1): incompatible dimensions");

  if      ( std::is_same<T,GFLOAT>::value ) {
    fmxm((GFLOAT*)y.data(), (GFLOAT*)D1.data().data(), &ND1, &ND2, (GFLOAT*)u.data(), &N1, &NYZ, &szMatCache_);
  }
  else if ( std::is_same<T,GDOUBLE>::value ) {
    dmxm((GDOUBLE*)y.data(), (GDOUBLE*)D1.data().data(), &ND1, &ND2, (GDOUBLE*)u.data(), &N1, &NYZ, &szMatCache_);
  }
  else if ( std::is_same<T,GQUAD>::value ) {
    qmxm((GQUAD*)y.data(), (GQUAD*)D1.data().data(), &ND1, &ND2, (GQUAD*)u.data(), &N1, &NYZ, &szMatCache_);
  }
  else {
    assert(FALSE);
  }

} // end of method I3_X_I2_X_D1 (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : I3_X_I2_X_D1 (2)
// DESC   : Apply tensor product operator to vector:
//            y = I3 X I2 X D1 u
// ARGS   : D1      : 1-direction (dense) operator 
//          u       : operand vector; must be at least of size
//                    N1*N2*N3*Ne, with 1 changing most rapidly, then, 2, then 3,
//                    with Ne 'elements' each of size N1*X2*N3
//          N1-N3   : coord dimensions of u, y 'elements', if interpreted 
//                    as matrix dimensions
//          Ne      : number of 'elements' in u, y
//          y   : return vector result of size >= N1 * N2 * N3
// RETURNS: none
//**********************************************************************************
template<typename T>
void I3_X_I2_X_D1(GTMatrix<T> &D1, GTVector<T> &u,
                  GSIZET N1, GSIZET N2, GSIZET N3, GSIZET Ne,
                  GTVector<T> &y)
{
	GEOFLOW_TRACE();
  GSIZET  ND1, ND2, NYZ, NN, Nu;

  ND1 = D1.size(1);
  ND2 = D1.size(2);
  NYZ = N2*N3;
  NN  = N1*N2*N3*Ne;
  Nu = NYZ * Ne;
  ASSERT_MSG(!(u.size() < NN || y.size() < NN), "GMTK::I3_X_I2_X_D1 (2): u or y of incorrect size");
  ASSERT_MSG(!(N1 != ND2), "GMTK::I3_X_I2_X_D1 (2): incompatible dimensions");

  if      ( std::is_same<T,GFLOAT>::value ) {
    fmxm((GFLOAT*)y.data(), (GFLOAT*)D1.data().data(), &ND1, &ND2, (GFLOAT*)u.data(), &N1, &Nu, &szMatCache_);
  }
  else if ( std::is_same<T,GDOUBLE>::value ) {
    dmxm((GDOUBLE*)y.data(), (GDOUBLE*)D1.data().data(), &ND1, &ND2, (GDOUBLE*)u.data(), &N1, &Nu, &szMatCache_);
  }
  else if ( std::is_same<T,GQUAD>::value ) {
    qmxm((GQUAD*)y.data(), (GQUAD*)D1.data().data(), &ND1, &ND2, (GQUAD*)u.data(), &N1, &Nu, &szMatCache_);
  }
  else {
    assert(FALSE);
  }

} // end of method I3_X_I2_X_D1 (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : I3_X_D2_X_I1 (1)
// DESC   : Apply tensor product operator to vector:
//            y = I3 X D2 X I1 u
// ARGS   : D2T     : 2-direction (dense) operator transpose
//          u       : operand vector; must be at least of size
//                    N1*N2*N3, with 1 changing most rapidly, then, 2, then 3
//          N1,N2,N3: coord dimensions of u, y
//          y   : return vector result of size >= N1 * N2 * N3
// RETURNS: none
//**********************************************************************************
template<typename T>
void I3_X_D2_X_I1(GTMatrix<T> &D2T, GTVector<T> &u, 
                  GSIZET N1, GSIZET N2, GSIZET N3,
                  GTVector<T> &y)
{
	GEOFLOW_TRACE();
  GSIZET  ND1, ND2, NXY, NN;

  ND1 = D2T.size(1);
  ND2 = D2T.size(2);
  NXY = N1*N2;
  NN  = N1*N2*N3;
  ASSERT_MSG(!( u.size() < NN || y.size() < NN ), "GMTK::I3_X_D2_X_I1 (1): incompatible dimensions");

  if      ( std::is_same<T,GFLOAT>::value ) {
    for ( auto k=0; k<N3; k++ ) {
      fmxm((GFLOAT*)(y.data()+k*NXY), (GFLOAT*)(u.data()+k*NXY), &N1, &N2, (GFLOAT*)D2T.data().data(), 
           &ND1, &ND2, &szMatCache_);
    }
  }
  else if ( std::is_same<T,GDOUBLE>::value ) {
    for ( auto k=0; k<N3; k++ ) {
      dmxm((GDOUBLE*)(y.data()+k*NXY), (GDOUBLE*)(u.data()+k*NXY), &N1, &N2, (GDOUBLE*)D2T.data().data(), 
           &ND1, &ND2, &szMatCache_);
    }
  }
  else if ( std::is_same<T,GQUAD>::value ) {
    for ( auto k=0; k<N3; k++ ) {
      qmxm((GQUAD*)(y.data()+k*NXY), (GQUAD*)(u.data()+k*NXY), &N1, &N2, (GQUAD*)D2T.data().data(), 
           &ND1, &ND2, &szMatCache_);
    }
  }
  else {
    assert(FALSE);
  }


} // end of method I3_X_D2_X_I1 (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : I3_X_D2_X_I1 (2)
// DESC   : Apply tensor product operator to vector:
//            y = I3 X D2 X I1 u
// ARGS   : D2T     : 2-direction (dense) operator transpose
//          u       : operand vector; must be at least of size
//                    N1*N2*N3*Ne, with 1 changing most rapidly, then, 2, then 3,
//                    with Ne 'elements' each of size N1*X2*N3
//          N1-N3   : coord dimensions of u, y 'elements', if interpreted 
//                    as matrix dimensions
//          Ne      : number of 'elements' in u, y
//          y   : return vector result of size >= N1 * N2 * N3 * Ne
// RETURNS: none
//**********************************************************************************
template<typename T>
void I3_X_D2_X_I1(GTMatrix<T> &D2T, GTVector<T> &u, 
                  GSIZET N1, GSIZET N2, GSIZET N3, GSIZET Ne,
                  GTVector<T> &y)
{
	GEOFLOW_TRACE();
  GSIZET  ND1, ND2, NXY, NN, Nu;

  ND1 = D2T.size(1);
  ND2 = D2T.size(2);
  NXY = N1*N2;
  NN  = N1*N2*N3*Ne;
  Nu  = N1*N2*N3;
  ASSERT_MSG(!( u.size() < NN || y.size() < NN ), "GMTK::I3_X_D2_X_I1 (2): incompatible dimensions");

  if      ( std::is_same<T,GFLOAT>::value ) {
    for ( auto j=0; j<Ne; j++ ) {
      for ( auto k=0; k<N3; k++ ) {
        fmxm((GFLOAT*)(y.data()+k*NXY)+j*Nu, (GFLOAT*)(u.data()+k*NXY)+j*Nu, &N1, &N2, (GFLOAT*)D2T.data().data(), &ND1, &ND2, &szMatCache_);
      }
    }
  }
  else if ( std::is_same<T,GDOUBLE>::value ) {
    for ( auto j=0; j<Ne; j++ ) {
      for ( auto k=0; k<N3; k++ ) {
        dmxm((GDOUBLE*)y.data()+k*NXY+j*Nu, (GDOUBLE*)u.data()+k*NXY+j*Nu, &N1, &N2, (GDOUBLE*)D2T.data().data(), &ND1, &ND2, &szMatCache_);
      }
    }
  }
  else if ( std::is_same<T,GQUAD>::value ) {
    for ( auto j=0; j<Ne; j++ ) {
      for ( auto k=0; k<N3; k++ ) {
        qmxm((GQUAD*)(y.data()+k*NXY)+j*Nu, (GQUAD*)(u.data()+k*NXY)+j*Nu, &N1, &N2, (GQUAD*)D2T.data().data(), &ND1, &ND2, &szMatCache_);
      }
    }
  }
  else {
    assert(FALSE);
  }


} // end of method I3_X_D2_X_I1 (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : D3_X_I2_X_I1 (1)
// DESC   : Apply tensor product operator to vector:
//            y = D3 X I2 X I1 u
// ARGS   : D3T     : 3-direction (dense) operator transpose
//          u       : operand vector; must be at least of size
//                    N1*N2*N3, with 1 changing most rapidly, then, 2, then 3
//          N1,N2,N3: coord dimensions of u, y
//          y   : return vector result of size >= N1 * N2 * N3
// RETURNS: none
//**********************************************************************************
template<typename T>
void D3_X_I2_X_I1(GTMatrix<T> &D3T, GTVector<T> &u,
                  GSIZET N1, GSIZET N2, GSIZET N3, 
                  GTVector<T> &y)
{
	GEOFLOW_TRACE();
  GSIZET  ND1, ND2, NXY, NN;

  ND1 = D3T.size(1);
  ND2 = D3T.size(2);
  NXY = N1*N2;
  NN  = N1*N2*N3;
  ASSERT_MSG(!( u.size() < NN || y.size() < NN ), "GMTK::D3_X_I2_X_I1 (1): incompatible dimensions");

  if      ( std::is_same<T,GFLOAT>::value ) {
    fmxm((GFLOAT*)y.data(), (GFLOAT*)u.data(), &NXY, &N3, (GFLOAT*)D3T.data().data(), 
         &ND1, &ND2, &szMatCache_);
  }
  else if ( std::is_same<T,GDOUBLE>::value ) {
    dmxm((GDOUBLE*)y.data(), (GDOUBLE*)u.data(), &NXY, &N3, (GDOUBLE*)D3T.data().data(), 
         &ND1, &ND2, &szMatCache_);
  }
  else if ( std::is_same<T,GQUAD>::value ) {
    qmxm((GQUAD*)y.data(), (GQUAD*)u.data(), &NXY, &N3, (GQUAD*)D3T.data().data(), 
         &ND1, &ND2, &szMatCache_);
  }
  else {
    assert(FALSE);
  }

} // end of method D3_X_I2_X_I1 (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : D3_X_I2_X_I1 (2)
// DESC   : Apply tensor product operator to vector:
//            y = D3 X I2 X I1 u
// ARGS   : D3T     : 3-direction (dense) operator transpose
//          u       : operand vector; must be at least of size
//                    N1*N2*N3*Ne, with 1 changing most rapidly, then, 2, then 3,
//                    with Ne 'elements' each of size N1*X2*N3
//          N1-N3   : coord dimensions of u, y 'elements', if interpreted 
//                    as matrix dimensions
//          Ne      : number of 'elements' in u, y
//          y   : return vector result of size >= N1 * N2 * N3
// RETURNS: none
//**********************************************************************************
template<typename T>
void D3_X_I2_X_I1(GTMatrix<T> &D3T, GTVector<T> &u,
                  GSIZET N1, GSIZET N2, GSIZET N3, GSIZET Ne,
                  GTVector<T> &y)
{
	GEOFLOW_TRACE();
  GSIZET  ND1, ND2, NXY, NN, Nu;

  ND1 = D3T.size(1);
  ND2 = D3T.size(2);
  NXY = N1*N2;
  NN  = N1*N2*N3;
  Nu = N1*N2*N3;
  ASSERT_MSG(!( u.size() < NN || y.size() < NN ), "GMTK::D3_X_I2_X_I1 (2): incompatible dimensions");  

  if      ( std::is_same<T,GFLOAT>::value ) {
    for ( auto i=0; i<Ne; i++ ) {
      fmxm((GFLOAT*)y.data()+i*Nu, (GFLOAT*)u.data()+i*Nu, &NXY, &N3, (GFLOAT*)D3T.data().data(), &ND1, &ND2, &szMatCache_);
    }
  }
  else if ( std::is_same<T,GDOUBLE>::value ) {
    for ( auto i=0; i<Ne; i++ ) {
      dmxm((GDOUBLE*)y.data()+i*Nu, (GDOUBLE*)u.data()+i*Nu, &NXY, &N3, (GDOUBLE*)D3T.data().data(), &ND1, &ND2, &szMatCache_);
    }
  }
  else if ( std::is_same<T,GQUAD>::value ) {
    for ( auto i=0; i<Ne; i++ ) {
      qmxm((GQUAD*)y.data()+i*Nu, (GQUAD*)u.data()+i*Nu, &NXY, &N3, (GQUAD*)D3T.data().data(), 
           &ND1, &ND2, &szMatCache_);
    }
  }
  else {
    assert(FALSE);
  }

} // end of method D3_X_I2_X_I1 (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : Dg3_X_Dg2_X_D1 
// DESC   : Apply tensor product operator to vector:
//            y = diag(D3) X diag(D2) X D1 u
// ARGS   : D1   : GTMatrix of 1-operator
//          Dg1  : diag(D2)
//          Dg3  : diag(D3)
//          u    : GVector argument
//          tmp  : tmp space, resized here if necessary
//          y    : GVector result
// RETURNS: none
//**********************************************************************************
template<typename T>
void Dg3_X_Dg2_X_D1(GTMatrix<T> &D1, GTVector<T> &Dg2, GTVector<T> &Dg3,
                    GTVector<T> &u, GTVector<T> &tmp, GTVector<T> &y)
{
	GEOFLOW_TRACE();
  GSIZET N11, N12, N2, N3, NN, NXY, NYZ;

  N11 = D1.size(1);
  N12 = D1.size(2);
  N2  = Dg2.size();
  N3  = Dg3.size();
  NXY = N11*N2;
  NYZ = N2*N3;
  NN  = N11*N2*N3;
  ASSERT_MSG(!( u.size() < NN  || y.size() < NN ),"Dg3_X_D2_X_Dg1: incompatible vectors");

  tmp.resizem(NXY*N3);

  if      ( std::is_same<T,GFLOAT>::value ) {
    // tmp = I X I X D1 u:
    fmxm((GFLOAT*)y.data(), (GFLOAT*)D1.data().data(), &N11, &N12, (GFLOAT*)u.data(), &N11, &NYZ, &szMatCache_);

    // tmp1 = I X Diag(D2) X I tmp:
    for ( auto k=0; k<N3; k++ ) {
      fmxDm((GFLOAT*)(tmp.data()+k*NXY),  (GFLOAT*)(y.data()+k*NXY), &N11, &N12, (GFLOAT*)Dg2.data(), &N2, &szMatCache_);
    }

    // y = Dg3 X I X I tmp1:
    fmxDm((GFLOAT*)y.data(), (GFLOAT*)tmp.data(), &NXY, &N3, (GFLOAT*)Dg3.data(), &N3, &szMatCache_);
  }
  else if ( std::is_same<T,GDOUBLE>::value ) {
    // tmp = I X I X D1 u:
    dmxm((GDOUBLE*)y.data(), (GDOUBLE*)D1.data().data(), &N11, &N12, (GDOUBLE*)u.data(), &N11, &NYZ, &szMatCache_);

    // tmp1 = I X Diag(D2) X I tmp:
    for ( auto k=0; k<N3; k++ ) {
      dmxDm((GDOUBLE*)(tmp.data()+k*NXY),  (GDOUBLE*)(y.data()+k*NXY), &N11, &N12, (GDOUBLE*)Dg2.data(), &N2, &szMatCache_);
    }

    // y = Dg3 X I X I tmp1:
    dmxDm((GDOUBLE*)y.data(), (GDOUBLE*)tmp.data(), &NXY, &N3, (GDOUBLE*)Dg3.data(), &N3, &szMatCache_);
  }
  else if ( std::is_same<T,GQUAD>::value ) {
    // tmp = I X I X D1 u:
    qmxm((GQUAD*)y.data(), (GQUAD*)D1.data().data(), &N11, &N12, (GQUAD*)u.data(), &N11, &NYZ, &szMatCache_);

    // tmp1 = I X Diag(D2) X I tmp:
    for ( auto k=0; k<N3; k++ ) {
      qmxDm((GQUAD*)(tmp.data()+k*NXY),  (GQUAD*)(y.data()+k*NXY), &N11, &N12, (GQUAD*)Dg2.data(), &N2, &szMatCache_);
    }

    // y = Dg3 X I X I tmp1:
    qmxDm((GQUAD*)y.data(), (GQUAD*)tmp.data(), &NXY, &N3, (GQUAD*)Dg3.data(), &N3, &szMatCache_);
  }
  else {
    assert(FALSE);
  }

} // end, method Dg3_X_Dg2_X_D1


//**********************************************************************************
//**********************************************************************************
// METHOD : Dg3_X_D2_X_Dg1 
// DESC   : Apply tensor product operator to vector:
//            y = diag(D3) X D2 X diag(D1) u
// ARGS   : Dg1  : transpose of dense matrix, D1
//          D2T  : GTMatrix for 2-operator, transpose
//          Dg3  : diag(D3)
//          u    : GVector argument
//          tmp  : tmp space, resized here if necessary
//          y    : GVector result
// RETURNS: none
//**********************************************************************************
template<typename T>
void Dg3_X_D2_X_Dg1(GTVector<T> &Dg1, GTMatrix<T> &D2T, GTVector<T> &Dg3,
                    GTVector<T> &u, GTVector<T> &tmp, GTVector<T> &y)
{
	GEOFLOW_TRACE();
  GSIZET N1, N21, N22, N3, NN, NXY, NYZ;

  N1  = Dg1.size();
  N21 = D2T.size(1);
  N22 = D2T.size(2);
  N3  = Dg3.size();
  NXY = N1*N21;
  NYZ = N21*N3;
  NN  = N1*N21*N3;
  ASSERT_MSG(!( u.size() < NN  || y.size() < NN ), "Dg3_X_D2_X_Dg1: incompatible vectors");

  tmp.resizem(NXY*N3);

  if      ( std::is_same<T,GFLOAT>::value ) {
    // tmp = I X I X Diag(D1) u:
    fDmxm((GFLOAT*)y.data(), (GFLOAT*)Dg1.data(), &N1, (GFLOAT*)u.data(), &N1, &NYZ, &szMatCache_);

    // tmp1 = I X D2 X I tmp:
    for ( auto k=0; k<N3; k++ ) {
      fmxm((GFLOAT*)(tmp.data()+k*NXY), (GFLOAT*)(y.data()+k*NXY), &N1, &N21, (GFLOAT*)D2T.data().data(), &N21, &N22, &szMatCache_);
    }

    // y = Dg3 X I X I tmp1:
    fmxDm((GFLOAT*)y.data(),  (GFLOAT*)tmp.data(), &NXY, &N3, (GFLOAT*)Dg3.data(), &N3, &szMatCache_);
  }
  else if ( std::is_same<T,GDOUBLE>::value ) {
    // tmp = I X I X Diag(D1) u:
    dDmxm((GDOUBLE*)y.data(), (GDOUBLE*)Dg1.data(), &N1, (GDOUBLE*)u.data(), &N1, &NYZ, &szMatCache_);

    // tmp1 = I X D2 X I tmp:
    for ( auto k=0; k<N3; k++ ) {
      dmxm((GDOUBLE*)(tmp.data()+k*NXY), (GDOUBLE*)(y.data()+k*NXY), &N1, &N21, (GDOUBLE*)D2T.data().data(), &N21, &N22, &szMatCache_);
    }

    // y = Dg3 X I X I tmp1:
    dmxDm((GDOUBLE*)y.data(),  (GDOUBLE*)tmp.data(), &NXY, &N3, (GDOUBLE*)Dg3.data(), &N3, &szMatCache_);
  }
  else if ( std::is_same<T,GQUAD>::value ) {
    // tmp = I X I X Diag(D1) u:
    qDmxm((GQUAD*)y.data(), (GQUAD*)Dg1.data(), &N1, (GQUAD*)u.data(), &N1, &NYZ, &szMatCache_);

    // tmp1 = I X D2 X I tmp:
    for ( auto k=0; k<N3; k++ ) {
      qmxm((GQUAD*)(tmp.data()+k*NXY), (GQUAD*)(y.data()+k*NXY), &N1, &N21, (GQUAD*)D2T.data().data(), &N21, &N22, &szMatCache_);
    }

    // y = Dg3 X I X I tmp1:
    qmxDm((GQUAD*)y.data(),  (GQUAD*)tmp.data(), &NXY, &N3, (GQUAD*)Dg3.data(), &N3, &szMatCache_);
  }
  else {
    assert(FALSE);
  }

} // end, method Dg3_X_D2_X_Dg1


//**********************************************************************************
//**********************************************************************************
// METHOD : D3_X_Dg2_X_Dg1 
// DESC   : Apply tensor product operator to vector:
//            y = D3 X diag(D2) X diag(D1) u
// ARGS   : Dg1  : diag(D1)
//          Dg1  : diag(D2)
//          D3T  : GMatrix for 3-operator, transpose
//          u    : GVector argument
//          tmp  : tmp space, resized here if necessary
//          y    : GVector result
// RETURNS: none
//**********************************************************************************
template<typename T>
void D3_X_Dg2_X_Dg1(GTVector<T> &Dg1, GTVector<T> &Dg2, GTMatrix<T> &D3T,
                    GTVector<T> &u, GTVector<T> &tmp, GTVector<T> &y)
{
	GEOFLOW_TRACE();
  GSIZET N1, N2, N31, N32, NN, NXY, NYZ;

  N1  = Dg1.size();
  N2  = Dg2.size();
  N31 = D3T.size(1);
  N32 = D3T.size(2);
  NXY = N1*N2;
  NYZ = N2*N31;
  NN  = N1*N2*N31;
  ASSERT_MSG(!( u.size() < NN  || y.size() < NN ), "D3_X_Dg2_X_Dg1: incompatible vectors");

  tmp.resizem(NXY*N31);

  if      ( std::is_same<T,GFLOAT>::value ) {
    // tmp = I X I X Diag(D1) u:
    fDmxm((GFLOAT*)y.data(), (GFLOAT*)Dg1.data(), &N1, (GFLOAT*)u.data(), &N1, &NYZ, &szMatCache_);

    // tmp1 = I X D2 X I tmp:
    for ( auto k=0; k<N31; k++ ) {
      fmxDm((GFLOAT*)(tmp.data()+k*NXY), (GFLOAT*)(y.data()+k*NXY), &N1, &N2, (GFLOAT*)Dg2.data(), &N2, &szMatCache_);
    }

    // y = Dg3 X I X I tmp1:
    fmxm((GFLOAT*)y.data(), (GFLOAT*)tmp.data(), &NXY, &N31, (GFLOAT*)D3T.data().data(), &N31, &N32, &szMatCache_);
  }
  else if ( std::is_same<T,GDOUBLE>::value ) {
    // tmp = I X I X Diag(D1) u:
    dDmxm((GDOUBLE*)y.data(), Dg1.data(), &N1, (GDOUBLE*)u.data(), &N1, &NYZ, &szMatCache_);

    // tmp1 = I X D2 X I tmp:
    for ( auto k=0; k<N31; k++ ) {
      dmxDm((GDOUBLE*)(tmp.data()+k*NXY), (GDOUBLE*)(y.data()+k*NXY), &N1, &N2, (GDOUBLE*)Dg2.data(), &N2, &szMatCache_);
    }

    // y = Dg3 X I X I tmp1:
    dmxm((GDOUBLE*)y.data(), (GDOUBLE*)tmp.data(), &NXY, &N31, (GDOUBLE*)D3T.data().data(), &N31, &N32, &szMatCache_);
  }
  else if ( std::is_same<T,GQUAD>::value ) {
    // tmp = I X I X Diag(D1) u:
    qDmxm((GQUAD*)y.data(), (GQUAD*)Dg1.data(), &N1, (GQUAD*)u.data(), &N1, &NYZ, &szMatCache_);

    // tmp1 = I X D2 X I tmp:
    for ( auto k=0; k<N31; k++ ) {
      qmxDm((GQUAD*)tmp.data()+k*NXY, (GQUAD*)y.data()+k*NXY, &N1, &N2, (GQUAD*)Dg2.data(), &N2, &szMatCache_);
    }

    // y = Dg3 X I X I tmp1:
    qmxm((GQUAD*)y.data(), (GQUAD*)tmp.data(), &NXY, &N31, (GQUAD*)D3T.data().data(), &N31, &N32, &szMatCache_);
  }
  else {
    assert(FALSE);
  }

} // end, method D3_X_Dg2_X_Dg1


//**********************************************************************************
//**********************************************************************************
// METHOD : add 
// DESC   : point-by-point addition, returned in specified GTVector:
//            vret = a*va + b*vb, for GFLOAT types
// ARGS   : vret: GTVector<T> return vector
//          va  : first vector operatnd
//          vb  : second vector operatnd
//          a,b: factors multiplying va and vb, respectively
// RETURNS: GTVector & 
//**********************************************************************************
#pragma acc routine vector
template<typename T>
void add(GTVector<T> &vret, const GTVector<T> &va, const GTVector<T> &vb, T a, T b) 
{
	GEOFLOW_TRACE();
  ASSERT_MSG(!( va.size() < vret.size() || vb.size() < vret.size() ), "GTVector<T>::add: incompatible size");

#if defined(GEOFLOW_USE_CBLAS) || defined(GEOFLOW_USE_CUBLAS)

  for ( auto j=vret.getIndex().beg(); j<=vret.getIndex().end(); j+=vret.getIndex().stride() ) {
    vret[j] = a*va[j] + b*vb[j];
  }

#else

  GSIZET nn = vret.getIndex().end() - vret.getIndex().beg() + 1;
  if      ( std::is_same<T,GFLOAT>::value ) {
   fzaxpby(vret.data(), (GFLOAT*)(va.data()), &a, 
           (GFLOAT*)(vb.data()), &b, &nn, &szVecCache_);
  }
  else if ( std::is_same<T,GDOUBLE>::value ) {
   dzaxpby(vret.data(), (GFLOAT*)(va.data()), &a, 
           (GDOUBLE*)(vb.data()), &b, &nn, &szVecCache_);
  }
  else if ( std::is_same<T,GQUAD>::value ) {
   qzaxpby(vret.data(), (GFLOAT*)(va.data()), &a, 
           (GQUAD*)(vb.data()), &b, &nn, &szVecCache_);
  }
  else {
    assert(FALSE);
  }

#endif

} // end, add 


//**********************************************************************************
//**********************************************************************************
// METHOD : operator * mat-vec (float)
// DESC   : matrix-vector product, returns product,
//          without destroying *this data, for GFLOAT type
// ARGS   : GTVector &
// RETURNS: none
//**********************************************************************************
template<typename T>
void matvec_prod(GTVector<T> &vret, const GTMatrix<T> &A, const GTVector<T> &b) 
{
	GEOFLOW_TRACE();
  ASSERT_MSG(!( b.size() < A.size(2) ), "GMTK::matvec_prod: incompatible size");

#if defined(GEOFLOW_USE_CBLAS) || defined(GEOFLOW_USE_CUBLAS)

  GSIZET n1 = A.size(1);
  GSIZET n2 = A.size(2);
  for ( auto i=0; i<n1; i++ ) {
     vret[i] = 0;
     for ( auto j=0; j<n2; j++ ) {
       vret[i] += A(i,j) * b(j);
     }
   }

#else

  GSIZET n1 = A.size(1);
  GSIZET n2 = A.size(2);
  if      ( std::is_same<T,GFLOAT>::value ) {
    fmxv(vret.data(), (GFLOAT*)(A.data().data()), 
                      (GFLOAT*)(b.data())       , 
                      &n1, &n2, &szMatCache_);
  }
  else if ( std::is_same<T,GDOUBLE>::value ) {
    dmxv(vret.data(), (GDOUBLE*)(A.data().data()), 
                      (GDOUBLE*)(b.data())       , 
                      &n1, &n2, &szMatCache_);
  }
  else if ( std::is_same<T,GQUAD>::value ) {
    qmxv(vret.data(), (GQUAD*)(A.data().data()), 
                      (GQUAD*)(b.data())       , 
                      &n1, &n2, &szMatCache_);
  }
  else {
    assert(FALSE);
  }

#endif


} // end of operator * mat-vec


//**********************************************************************************
//**********************************************************************************
// METHOD : mat-mat prod 
// DESC   : Multiplies C = A * B and returns matrix prod C
//          result, for GFLOAT types
// ARGS   : GTMatrix m factor
// RETURNS: none
//**********************************************************************************
template<typename T>
void matmat_prod(GTMatrix<T> &C, const GTMatrix<T> &A, const GTMatrix<T> &B) 
{
	GEOFLOW_TRACE();
  ASSERT_MSG(!( A.size(2) != B.size(1) ), "GMTK::matmat_prod:incompatible matrix");

#if defined(GEOFLOW_USE_CBLAS) || defined(GEOFLOW_USE_CUBLAS)

  for ( auto i=0; i<C.size(1); i++ ) {
    for ( auto j=0; j<C.size(2); j++ ) {
      C(i,j) = 0.0;
      for ( auto k=0; k<A.size(2); k++ ) {
        C(i,j) += A(i,k) * B(k,j);
      }
    }
  }

#else

  GSIZET a1=A.size(1), a2 = A.size(2);
  GSIZET b1=B.size(1), b2 = B.size(2);
  if      ( std::is_same<T,GFLOAT>::value ) {
    fmxm(C.data().data(),
         (GFLOAT*)(A.data().data()),
         &a1,&a2,
         (GFLOAT*)(B.data().data()),
         &b1, &b2, &szMatCache_);
  }
  else if ( std::is_same<T,GDOUBLE>::value ) {
    dmxm(C.data().data(),
         (GDOUBLE*)(A.data().data()),
         &a1,&a2,
         (GDOUBLE*)(B.data().data()),
         &b1, &b2, &szMatCache_);
  }
  else if ( std::is_same<T,GQUAD>::value ) {
    qmxm(C.data().data(),
         (GQUAD*)(A.data().data()),
         &a1,&a2,
         (GQUAD*)(B.data().data()),
         &b1, &b2, &szMatCache_);
  }
  else {
    assert(FALSE);
  }

#endif

} // end of operator * mat-mat


//**********************************************************************************
//**********************************************************************************
// METHOD : dot
// DESC   : Compute dot prod of input vector fields,
//          with the understanding that either vector 
//          could have NULL or constant components:
//                r = f . v,
//          where return vector, r, is non-NULL and of
//          full length. A constant vector has a size of 1
// ARGS   : 
//          x  : first vector field; may be NULL or constant
//          y  : second vector  field; may be NULL or constant
//          tmp: tmp vector, of full size
//          r  : scalar field containing injection rate, of full size
// RETURNS: none.
//**********************************************************************************
template<typename T>
void dot(GTVector<GTVector<T>*> &x, GTVector<GTVector<T>*> &y, GTVector<T> &tmp, GTVector<T> &r)
{
	GEOFLOW_TRACE();
   assert(x.size() == y.size());

   r = 0.0;
   for ( auto j=0; j<x.size(); j++ ) {

     if ( x[j] == NULLPTR || y[j] == NULLPTR ) continue;
     if      ( x[j]->size() == 1 && y[j]->size() == 1 ) { // x , y constant
       tmp =  (*x[j])[0] * (*y[j])[0];
     }
     else if ( x[j]->size() >  1 && y[j]->size() == 1 ) { // y constant
       GMTK::saxpby<T>(tmp, 1.0, *x[j], (*y[j])[0]);  
       
     }
     else if ( x[j]->size() == 1 && y[j]->size() >  1 ) { // x constant
       GMTK::saxpby<T>(tmp, 1.0, *y[j], (*x[j])[0]);  
     }
     else {                                               // x, y of full length
       x[j]->pointProd(*y[j], tmp);
     }
     r += tmp;

   }

} // end of method dot


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
template<typename Grid, typename T>
void curl(Grid &grid, const GTVector<GTVector<T>*> &u, const GINT idir, 
          GTVector<GTVector<T>*> &tmp, GTVector<T> &curlc)
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
//                to refer to the Cartesian 'deriv' method in Grid
//          
// ARGS   : grid : grid
//          u    : input (scalar) field. 
//          idir : gradient component to compute. Must be appropriate for 
//                 problem dimension.
//          tmp  : tmp vector; must be of at least length 1.
//          gradc: result
// RETURNS: none.
//**********************************************************************************
template<typename Grid, typename T>
void grad(Grid &grid, GTVector<T> &u, const GINT idir, 
          GTVector<GTVector<T>*> &tmp, GTVector<T> &gradc)
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
template<typename Grid, typename T>
void constrain2sphere(Grid &grid, const GTVector<GTVector<T>*> &v, GTVector<GTVector<T>*> &Pv)
{
  GEOFLOW_TRACE();
  if ( grid.gtype() != GE_2DEMBEDDED || v.size() != 3 ) return;

  assert( v.size() >= 3 && "Incompatible dimensionality");

  GSIZET nxy = grid.ndof();
  T      r2, x, y, z;
  GTVector<GTVector<T>> *xnodes = &grid.xNodes();

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
template<typename Grid, typename T>
void constrain2sphere(Grid &grid, GTVector<GTVector<T>*> &v)
{
  GEOFLOW_TRACE();
  if ( grid.gtype() != GE_2DEMBEDDED || v.size() < 3 ) return;

  GSIZET nxy = grid.ndof();
  T      ri, r2, x, y, z;
  GTVector<GTVector<T>> *xnodes = &grid.xNodes();

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
template<typename Grid, typename T>
void vsphere2cart(Grid &grid, const GTVector<GTVector<T>*> &vsph, GVectorType vtype, GTVector<GTVector<T>*> &vcart)
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


GSIZET      nxy = grid.ndof();
T           x, y, z, tiny;
T           phi, r, theta;
T           gpp, gtt;
T           vthcontra, vphicontra;
GTVector<GTVector<T>> *xnodes = &grid.xNodes();

tiny = numeric_limits<T>::epsilon();

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
template<typename Grid, typename T>
void vcart2sphere(Grid &grid, const GTVector<GTVector<T>*> &vcart, GVectorType vtype, GTVector<GTVector<T>*> &vsph)
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

GSIZET      nxy = grid.ndof();
T           x, y, z, tiny;
T           phi, r, theta;
T           gpp, gtt;
T           vthcontra, vphicontra;
GTVector<GTVector<T>> *xnodes = &grid.xNodes();

tiny = numeric_limits<T>::epsilon();

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
template<typename T>
void cart2latlon(const GTVector<GTVector<T>*> &cart, GTVector<GTVector<T>*> &latlon)
{
  GEOFLOW_TRACE();
assert( cart.size() >= 3 && latlon.size() >= 2 && "Must have correct array sizes");


GSIZET      nxy = cart[0]->size();
T           x, y, z;
T           phi, r, theta;

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
template<typename T>
void rcart2sphere(const GTVector<GTVector<T>*> &cart, GTVector<GTVector<T>*> &rlatlon)
{
  GEOFLOW_TRACE();
assert( cart.size() >= 3 && rlatlon.size() >= 3 && "Must have correct array sizes");


GSIZET           nxy = cart[0]->size();
T           x, y, z;
T           phi, r, theta;

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
// RETURNS: T energy
//**********************************************************************************
template<typename Grid, typename T>
T energy(Grid &grid, const GTVector<GTVector<T>*> & u, GTVector<GTVector<T>*> &tmp, GBOOL isglobal, GBOOL ismax)
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

return static_cast<T>(ener);

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
// RETURNS: T enstrophy
//**********************************************************************************
template<typename Grid, typename T>
T enstrophy(Grid &grid, const GTVector<GTVector<T>*> & u, GTVector<GTVector<T>*> &tmp, GBOOL isglobal, GBOOL ismax)
{
    GEOFLOW_TRACE();
assert(tmp.size() >= 4 && "Insufficient temp space");


GINT                   ibeg, iend;
GDOUBLE                enst, local;
GC_COMM                comm = grid.get_comm();
GTVector<T>           *cc;
GTVector<GTVector<T>*> utmp(3);

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

return static_cast<T>(enst);

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
// RETURNS: T helicity
//**********************************************************************************
template<typename Grid, typename T>
T helicity(Grid &grid, const GTVector<GTVector<T>*> & u, GTVector<GTVector<T>*> &tmp, GBOOL isglobal, GBOOL ismax)
{
    GEOFLOW_TRACE();
assert(tmp.size() >= 4 && "Insufficient temp space");


GDOUBLE                hel, local;
GC_COMM                comm = grid.get_comm();
GTVector<T>           *cc;
GTVector<GTVector<T>*> utmp(3);

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

return static_cast<T>(hel);

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
// RETURNS: T helicity
//**********************************************************************************
template<typename Grid, typename T>
T relhelicity(Grid &grid, const GTVector<GTVector<T>*> & u, GTVector<GTVector<T>*> &tmp, GBOOL isglobal, GBOOL ismax)
{
    GEOFLOW_TRACE();
assert(tmp.size() >= 5 && "Insufficient temp space");


GDOUBLE                local, rhel;
GC_COMM                comm = grid.get_comm();
GTVector<T>           *cc;
GTVector<GTVector<T>*> utmp(3);

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
T tiny = 100.0*numeric_limits<T>::epsilon();

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

return static_cast<T>(rhel);

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
// RETURNS: T energy injection rate
//**********************************************************************************
template<typename Grid, typename T>
T energyinj(Grid &grid, const GTVector<GTVector<T>*> &u,  const GTVector<GTVector<T>*> &uf, GTVector<GTVector<T>*> &tmp, GBOOL isglobal, GBOOL ismax)
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

return static_cast<T>(einj);

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
template<typename Grid, typename T>
void domathop(Grid &grid, const GTVector<GTVector<T>*> &uin,  const GString sop, GTVector<GTVector<T>*> &utmp, GTVector<GTVector<T>*> &uout, std::vector<GINT> &iuout)
{
  GEOFLOW_TRACE();
  GINT                        ib, nxy;
  GTVector<GTVector<T>*>      tmp(3);

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
