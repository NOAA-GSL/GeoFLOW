//==================================================================================
// Module       : gmtk.ipp
// Date         : 1/31/18 (DLR)
// Description  : Math TooKit: namespace of C routines for various
//                math functions
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

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
  T fact;
  T pmm, pmmp1, pll, somx2;
  T colat, phi, r;
  T x, y, z;

  assert( m >= 0 && m <= l );

  for ( auto j=0; j<xnodes[0].size(); j++ ) {
    x     = xnodes[0][j]; y = xnodes[1][j]; z = xnodes[2][j];
    r     = sqrt(x*x +y*y + z*z);

    colat = acos(z/r); // Note: this is co-latitude
    phi   = atan2(y,x);
    x     = cos(colat);
    assert( fabs(x) <= 1.0  );
    pmm=1.0; // compute Pm_m 
    if ( m > 0 ) {
      somx2 = sqrt((1.0-x)*(1.0+x)); 
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
      pmmp1=x*(2*m+1)*pmm; 
      if (l == (m+1)) {
        plm[j] = pmmp1; 
        continue;
      }
      else { // compute P^m_l; l > m+1
        for ( auto ll=m+2; ll<=l; ll++) {
          pll   = (x*(2*ll-1)*pmmp1 - (ll+m-1)*pmm)/(ll-m); 
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
//              Ylm  =  sqrt{ (2l+1)(l-m)! / [4pi (l+m)!] } X
//                      P_lm exp(i m phi)
//          If iri = 0; else returns dYlm^*. If |m| > l, set dYlm=0
//          
//   
// ARGS   : l,m,  : orbital ang mom. quantum number, and azimuthal quantum number
//          xnodes: Cartesian coordinate arrays
//          iri   : if 0, return dYlm, else, return complex conjugate, Ylm^*
//          ylm_r: real comp. of dY^l_m for each grid point
//          ylm_i: imag. comp. of dY^l_m for each grid point
// RETURNS: none
//**********************************************************************************
template<typename T>
void Ylm_cart(GINT l, GINT m, GTVector<GTVector<T>> &xnodes, GINT iri, GTVector<T> &ylm_r, GTVector<T> &ylm_i)
{
  T rfact, phi, r;
  T x, y, z;
  T xl, xm;

  assert( iri >= 0 );

  if ( abs(m) > l ) {
    ylm_r = 0.0;
    ylm_i = 0.0;
    return;
  }
  
  xl = l; xm = m;
  rfact = sqrt( (2.0*xl+1.0)/(4.0*PI) 
        *       GMTK::fact<T>(xl-xm)/GMTK::fact<T>(xl+m) );

  GMTK::Plm_cart<T>(l, abs(m), xnodes, ylm_r);
  for ( auto j=0; j<xnodes[0].size(); j++ ) {
    x     = xnodes[0][j]; y = xnodes[1][j]; z = xnodes[2][j];
    r     = sqrt(x*x + y*y + z*z);
    phi   = atan2(y,x);
    ylm_r[j] *= rfact*cos(phi);
    ylm_i[j] *= rfact*sin(phi);
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
//                            sqrt[(l-m)(l+m_1)] e^(-i phi) Yl(m+1)
//              dYlm/dphi   = i m Ylm 
//          If iri = 0, return dYlm; else return dYlm^*. 
//          
//   
// ARGS   : l,m,  : orbital ang mom. quantum number, and azimuthal quantum number
//          xnodes: Cartesian coordinate arrays
//          idir  : If 1, take d/d\theta; if 2, take d/d\phi
//          iri   : if 0, return dYlm, else, return complex conjugate, dYlm^*
//          tmp   : tmp arrays. 2 are required.
//          dylm_r: real comp. of dY^l_m for each grid point
//          dylm_i: imag. comp. of dY^l_m for each grid point
// RETURNS: none
//**********************************************************************************
template<typename T>
void dYlm_cart(GINT l, GINT m, GTVector<GTVector<T>> &xnodes, GINT idir,  GINT iri, GTVector<GTVector<T>*> &tmp, GTVector<T> &dylm_r, GTVector<T> &dylm_i)
{
  T colat, phi, r;
  T x, y, z;
  T xl, xm;
  T cs, cotc, eps, epsi, sn, rfact;

  assert( iri >= 0 );
  assert( idir == 1 || idir == 2 );
  assert( tmp.size() >= 2 );

  eps   = 100.0*std::numeric_limits<T>::epsilon();
  epsi  = 1.0/eps;
  xl    = l; xm = m;
  rfact = sqrt( (xl-xm)*(xl+xm+1.0) );

  if ( idir == 1 ) {

    Ylm_cart<T>(l, m  , xnodes, 0, dylm_r  , dylm_i);   // Ylm
    Ylm_cart<T>(l, m+1, xnodes, 0, *tmp[0], *tmp[1]); // Y^l_m+1
    for ( auto j=0; j<xnodes[0].size(); j++ ) {
      x     = xnodes[0][j]; y = xnodes[1][j]; z = xnodes[2][j];
      r     = sqrt(x*x + y*y + z*z);
      cs    = cos(colat); sn = sin(colat);
      cotc  = fabs(sn) < eps ? cs*epsi : cs/sn;
      colat = acos(z/r);
      phi   = atan2(y,x);
      dylm_r[j] *= xm*cotc*dylm_r[j] 
                 + rfact*( cos(phi)(*tmp[0])[j] + sin(phi)*(*tmp[1])[j] );
      dylm_i[j] *= xm*cotc*dylm_i[j]
                 + rfact*( cos(phi)(*tmp[1])[j] - sin(phi)*(*tmp[0])[j] );
    } // end, coord loop

  }
  else if ( idir == 2 ) {

    Ylm_cart<T>(l, m+1, xnodes, 0, *tmp[0], *tmp[1]);   // Ylm
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
//          Cmputed as:
//              dYlm/dtheta = m cot theta Ylm + 
//                            sqrt[(l-m)(l+m_1)] e^(-i phi) Yl(m+1)
//              dYlm/dphi   = i m Ylm 
//          If iri = 0, return dYlm; else return dYlm^*. 
//          
//   
// ARGS   : l,m,  : orbital ang mom. quantum number, and azimuthal quantum number
//          xnodes: Cartesian coordinate arrays
//          idir  : If 1, take d^2/d\theta^2; if 2, take d^2/d\phi^2
//          iri   : if 0, return dYlm, else, return complex conjugate, dYlm^*
//          tmp   : tmp arrays. 2 are required.
//          dylm_r: real comp. of dY^l_m for each grid point
//          dylm_i: imag. comp. of dY^l_m for each grid point
// RETURNS: none
//**********************************************************************************
template<typename T>
void ddYlm_cart(GINT l, GINT m, GTVector<GTVector<T>> &xnodes, GINT idir, GINT iri, GTVector<GTVector<T>*> &tmp, GTVector<T> &dylm_r, GTVector<T> &dylm_i)
{
  T colat, phi, r;
  T x, y, z;
  T xl, xm;
  T cotc, cs, csc, eps, epsi, sn, rfact;

  assert( iri >= 0 );
  assert( idir == 1 || idir == 2 );
  assert( tmp.size() >= 2 );

  eps  = 100.0*std::numeric_limits<T>::epsilon();
  epsi = 1.0/eps;
  xl   = l; xm = m;

  if ( idir == 1 ) {

    Ylm_cart<T>(l, m  , xnodes, 0, dylm_r  , dylm_i);   // Ylm
    Ylm_cart<T>(l, m+1, xnodes, 0, *tmp[0], *tmp[1]); // Y^l_m+1
    rfact = sqrt( (xl-xm)*(xl+xm-1.0) ) * (2.0*xm+1.0);
    for ( auto j=0; j<xnodes[0].size(); j++ ) {
      x     = xnodes[0][j]; y = xnodes[1][j]; z = xnodes[2][j];
      r     = sqrt(x*x + y*y + z*z);
      colat = acos(z/r);
      cs    = cos(colat); sn = sin(colat);
      cotc  = fabs(sn) < eps ? cs*epsi: cs/sn;
      csc   = fabs(sn) < eps ? epsi : 1.0/sn;

      phi   = atan2(y,x);
      dylm_r[j] = xm*(xm*cotc*cotc - csc*csc)*dylm_r[j] 
                + rfact*cotc*( cos(phi)(*tmp[0])[j] + sin(phi)*(*tmp[1])[j] );
      dylm_i[j] = xm*(xm*cotc*cotc - csc*csc)*dylm_i[j]
                + rfact*cotc*( cos(phi)(*tmp[1])[j] - sin(phi)*(*tmp[0])[j] );
    } // end, coord loop

    Ylm_cart<T>(l, m+2, xnodes, 0, *tmp[0], *tmp[1]); // Y^l_m+2
    rfact = sqrt( (xl-xm)*(xl-xm-1.0)*(xm+xl+2.0)*(xm+xl+1.0) );
    for ( auto j=0; j<xnodes[0].size(); j++ ) {
      x     = xnodes[0][j]; y = xnodes[1][j]; z = xnodes[2][j];
      phi   = atan2(y,x);
      dylm_r[j] += rfact*( cos(2.0*phi)(*tmp[0])[j] + sin(2.0*phi)*(*tmp[1])[j] );
      dylm_i[j] += rfact*( cos(2.0*phi)(*tmp[1])[j] - sin(2.0*phi)*(*tmp[0])[j] );
    } // end, coord loop

  }
  else if ( idir == 2 ) {

    Ylm_cart<T>(l, m+1, xnodes, 0, *tmp[0], *tmp[1]);   // Ylm
    for ( auto j=0; j<xnodes[0].size(); j++ ) {
      dylm_r[j] = -xm*xm * (*tmp[0])[j];
      dylm_i[j] =  xm*xm * (*tmp[1])[j];
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
//              rYlm  =  (-1)^m/  sqrt(2) (Ylm + Ylm^*); m > 0
//                                         Yl0         ; m = 0
//                       (-1)^m/(isqrt(2) (Ylm - Ylm^*); m < 0
//          
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
  T rfact;


  if ( abs(m) > l ) {
    rylm = 0.0;
    return;
  }

  if ( m > 0 ) { 
    rfact = pow(-1.0,m)*sqrt(2.0);
    GMTK::Ylm_cart<T>(l, m, xnodes, rylm, tmp);
    rylm *= rfact;
  }
  else if ( m < 0 ) {
    rfact = pow(-1.0,m)*sqrt(2.0);
    GMTK::Ylm_cart<T>(l, m, xnodes, tmp, rylm);
    rylm *= rfact;
  }
  else {
    GMTK::Ylm_cart<T>(l, m, xnodes, rylm, tmp);
  }

} // end, method rYlm_cart


//**********************************************************************************
//**********************************************************************************
// METHOD : drYlm_cart
// DESC   : Compute idir--derivative real spherical harmonics at each grid point
//          specified by xnodes in Carteswian coordinates.
//          Computed as:
//              rYlm  =  (-1)^m/  sqrt(2) (d_idir Ylm + d_idir Ylm^*); m > 0
//                                         d_idir Yl0                ; m = 0
//                       (-1)^m/(isqrt(2) (d_idir Ylm - d_idir Ylm^*); m < 0
//          
//   
// ARGS   : l,m,  : orbital ang mom. quantum number, and azimuthal quantum number
//          xnodes: Cartesian coordinate arrays
//          idir  : If 1, computes d/d\theta; if 2, computes d/d\phi
//          tmp   : tmp arrays. 3 are required.
//          dylm  : deriv of real comp. of Y^l_m for each grid point
// RETURNS: none
//**********************************************************************************
template<typename T>
void drYlm_cart(GINT l, GINT m, GTVector<GTVector<T>> &xnodes, GINT idir, GTVector<GTVector<T>*>  &tmp, GTVector<T> &drylm)
{
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
    rfact = pow(-1.0,m)*sqrt(2.0);
    GMTK::dYlm_cart<T>(l, m, xnodes, idir, tmp, drylm, *di);
    drylm *= rfact;
  }
  else if ( m < 0 ) {
    rfact = pow(-1.0,m)*sqrt(2.0);
   *dr = tmp[2];
    GMTK::dYlm_cart<T>(l, m, xnodes, idir, *dr, drylm);
    drylm *= rfact;
  }
  else {
   *di = tmp[2];
    GMTK::dYlm_cart<T>(l, m, xnodes, idir, drylm, *di);
  }

} // end, method drYlm_cart


//**********************************************************************************
//**********************************************************************************
// METHOD : ddrYlm_cart
// DESC   : Compute idir--2nd derivative of real spherical harmonics at 
//          each grid point specified by xnodes in Carteswian coordinates.
//          Computed as:
//              rYlm  =  (-1)^m/  sqrt(2) (d^2_idir Ylm + d^2_idir Ylm^*); m > 0
//                                         d^2_idir Yl0                  ; m = 0
//                       (-1)^m/(isqrt(2) (d^2_idir Ylm - d^2_idir Ylm^*); m < 0
//          
//   
// ARGS   : l,m,  : orbital ang mom. quantum number, and azimuthal quantum number
//          xnodes: Cartesian coordinate arrays
//          idir  : If 1, computes d/d\theta; if 2, computes d/d\phi
//          tmp   : tmp arrays. 3 are required.
//          dylm  : deriv. of real comp. of Y^l_m for each grid point
// RETURNS: none
//**********************************************************************************
template<typename T>
void ddrYlm_cart(GINT l, GINT m, GTVector<GTVector<T>> &xnodes, GINT idir, GTVector<GTVector<T>*>  &tmp, GTVector<T> &drylm)
{
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
    rfact = pow(-1.0,m)*sqrt(2.0);
    GMTK::ddYlm_cart<T>(l, m, xnodes, idir, tmp, drylm, *di);
    drylm *= rfact;
  }
  else if ( m < 0 ) {
    rfact = pow(-1.0,m)*sqrt(2.0);
   *dr = tmp[2];
    GMTK::ddYlm_cart<T>(l, m, xnodes, idir, *dr, drylm);
    drylm *= rfact;
  }
  else {
   *di = tmp[2];
    GMTK::ddYlm_cart<T>(l, m, xnodes, idir, drylm, *di);
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
  assert( A.size() >= 2 && C.size() >= 2 &&  "Incompatible dimensionality");

  GSIZET n;
  T      fact = isgn < 0 ? -1.0 : 1.0;
  if ( iind != NULLPTR ) {
    for ( GSIZET k=0; k<nind; k++ ) { // cycle over all coord pairs
      n = iind[k];
      C[0][k] =  (*A[1])[n]*fact;
      C[1][k] = -(*A[0])[n]*fact;
      if ( C.size() > 2 ) C[2][n] = 0.0;
    }
  }
  else {
    for ( GSIZET n=0; n<A[0]->size(); n++ ) { // cycle over all coord pairs
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

  GSIZET n;
  T      fact = isgn < 0 ? -1.0 : 1.0;
  if ( iind != NULLPTR ) {
    for ( GSIZET k=0; k<nind; k++ ) { // cycle over all coord pairs
      n = iind[k];
      Cx[k] =  Ay[n]*fact;
      Cy[k] = -Ax[n]*fact;
    }
  }
  else {
    for ( GSIZET n=0; n<Ax.size(); n++ ) { // cycle over all coord pairs
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
  assert( A.size() >= 3 && B.size() && C.size() >= 3 && "Incompatible dimensionality");

  GSIZET n;
  T      x1, y1, z1;
  T      x2, y2, z2;

  if ( iind != NULLPTR ) {
    for ( GSIZET k=0; k<nind; k++ ) { // cycle over all coord pairs
      n = iind[k];
      x1 = (*A[0])[n]; y1 = (*A[1])[n]; z1 = (*A[2])[n];
      x2 = (*B[0])[n]; y2 = (*B[1])[n]; z2 = (*B[2])[n];
      C[0][k] = y1*z2 - z1*y2; 
      C[1][k] = z1*x2 - z2*x1; 
      C[2][k] = x1*y2 - x2*y1;
    }
  }
  else {
    for ( GSIZET n=0; n<A[0]->size(); n++ ) { // cycle over all coord pairs
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

  GSIZET n;
  T      x1, y1, z1;
  T      x2, y2, z2;

  if ( iind != NULLPTR ) {
    for ( GSIZET k=0; k<nind; k++ ) { // cycle over all coord pairs
      n = iind[k];
      x1 = Ax[n]; y1 = Ay[n]; z1 = Az[n];
      x2 = Bx[n]; y2 = By[n]; z2 = Bz[n];
      Cx[k] = y1*z2 - z1*y2; 
      Cy[k] = z1*x2 - z2*x1; 
      Cz[k] = x1*y2 - x2*y1;
    }
  }
  else {
    for ( GSIZET n=0; n<Ax.size(); n++ ) { // cycle over all coord pairs
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
  GSIZET n;
  T      xn;

  // NOTE: The better way to do this, especially for long lists is use
  //       outer loop over coord dimensions, and store progressive sum
  //       over each for each tuple. But this requires vector storage for
  //       each coordinate for all tuples.
  xn = 0.0;
  if ( iind != NULLPTR ) {
    for ( GSIZET k=0; k<nind; k++ ) { // cycle over all n-tuples
      n = iind[k];
      for ( GSIZET l=0, xn=0.0; l<x.size(); l++ ) xn += (*x[l])[n]*(*x[l])[n];
      xn = x0/pow(xn,0.5);
      for ( GSIZET l=0; l<x.size(); l++ ) (*x[l])[n] *= xn;
    }
  }
  else {
    for ( GSIZET n=0; n<x[0]->size(); n++ ) { // cycle over all n-tuples
      for ( GSIZET l=0, xn=0.0; l<x.size(); l++ ) xn += (*x[l])[n]*(*x[l])[n];
      xn = x0/pow(xn,0.5);
      for ( GSIZET l=0; l<x.size(); l++ ) (*x[l])[n] *= xn;
    }
  }

} // end of method normalize_euclidean


//**********************************************************************************
//**********************************************************************************
// METHOD :    saxpby
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
  assert(x.size() == y.size() && "Incompatible array sizes");
  for ( GSIZET j=0; j<x.size(); j++ ) { 
    x[j] = a*x[j] + b*y[j];
  }
} // end of method saxpby


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
template<typename T>
void normalizeL2(GGrid &grid, GTVector<GTVector<T>*> &u, GTVector<GTVector<T>*> &tmp, T E0)
{
  GSIZET  n;
  GDOUBLE xn, xint;

  // xinit gives 0.5  Int _u_^2 dV / Int dV:
  xint = static_cast<GDOUBLE>(GMTK::energy(grid, u, tmp, TRUE, FALSE));
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
  T tiny=std::numeric_limits<T>::epsilon();

  for ( auto i=0; i<u.size(); i++ ) {
    for ( auto j=0; j<u[i].size(); j++ ) {
      if ( fabs(u[i][j]) < tiny ) u[i][j] = 0.0;
    }
  }


} // end of method zero (2)



} // end, namespace GMTK
