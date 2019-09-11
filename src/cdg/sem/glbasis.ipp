//==================================================================================
// Module       : glbasis.ipp
// Date         : 1/18/18 (DLR)
// Description  : Encapsulates the methods and data associated with
//                a spectral element nodal Gauss-Legendre basis
//                Note: This basis is the Gauss-Legendre
//                basis, defined as
//                  h(xi)_j = L_N(xi) / ( dL_N(xi_j)/dxi * (xi-xi_j) )
//                where N is the order of the expansion, and L_N is the Legendre
//                polynomial of order N, and the xi_j are the nodal points. L_N
//                is obtained from the generalized Jacobi polynomial computed in
//                computeJacobi.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : GLLBasis
//==================================================================================
#include "glbasis.hpp"
#include <cstdlib>
#include <memory>
#include <cmath>
#include <cstdio>


//************************************************************************************
//************************************************************************************
// METHOD : Contructor (1)
// DESC   : Instantiate with polynomial order, and maxOrder
//          
// ARGS   :
// RETURNS: TRUE on success, else FALSE
//************************************************************************************
template<typename T, typename TE>
GLBasis<T,TE>::GLBasis(GINT  inOrder, GINT  MaxOrder)
: GLLBasis<T,TE>(inOrder,MaxOrder)
{
  resize(Np_);

} // end of constructor method


//************************************************************************************
//************************************************************************************
// METHOD : Contructor (2)
// DESC   : Default constructor
//          
// ARGS   :
// RETURNS:
//************************************************************************************
template<typename T, typename TE>
GLBasis<T,TE>::GLBasis()
: GLLBasis<T,TE>()
{
  resize(Np_);
} // end of constructor method

//************************************************************************************
//************************************************************************************
// METHOD : Contructor (3)
// DESC   : Instantiate with poly order (same for maxOrder)
//          
// ARGS   :
// RETURNS:
//************************************************************************************
template<typename T, typename TE>
GLBasis<T,TE>::GLBasis(GINT  Order)
: GLLBasis<T,TE>(Order)
{
  resize(Np_);
} // end of constructor method


//************************************************************************************
//************************************************************************************
// METHOD : Copy consructor
// DESC   : 
//          
// ARGS   :
// RETURNS: 
//************************************************************************************
template<typename T, typename TE>
GLBasis<T,TE>::GLBasis(const GLBasis &b)
{
   // copy data:
    alpha_    = b.alpha_;
    beta_     = b.beta_;
    ximin_    = b.ximin_;
    ximax_    = b.ximax_;
    Np_       = b.Np_;
    bInit_    = b.bInit_;
    bNeedDerivMatrix_ = b.bNeedDerivMatrix_;
    bNeedLegMatrix_  = b.bNeedLegMatrix_;

    //  matrices, and node data:
    xiNodes_     = b.xiNodes_;
    weights_     = b.weights_;
    Pn_          = b.Pn_;
    dPn_         = b.dPn_;
    Phi_         = b.Phi_;
    dPhi_        = b.dPhi_;
    dPhiT_       = b.dPhiT_;
    stiffMatrix_ = b.stiffMatrix_;
    LegMatrix_   = b.LegMatrix_;

}

//************************************************************************************
//************************************************************************************
// METHOD : Destructor
// DESC   : 
// ARGS   :
// RETURNS: 
//************************************************************************************
template<typename T, typename TE>
GLBasis<T,TE>::~GLBasis()
{
}

//************************************************************************************
//************************************************************************************
// METHOD : Assignment operator
// DESC   : 
// ARGS   :
// RETURNS: 
//************************************************************************************
template<typename T, typename TE>
void GLBasis<T,TE>::operator=(const GLBasis &b)
{
  if ( &b != this ) 
  {
   // copy data:
    ximin_    = b.ximin_;
    ximax_    = b.ximax_;
    Np_       = b.Np_;
    bInit_    = b.bInit_;
    bNeedDerivMatrix_ = b.bNeedDerivMatrix_;
    bNeedLegMatrix_  = b.bNeedLegMatrix_;

    //  matrices, and node data:
    xiNodes_     = b.xiNodes_;
    weights_     = b.weights_;
    Pn_          = b.Pn_;
    dPn_         = b.dPn_;
    Phi_         = b.Phi_;
    dPhi_        = b.dPhi_;
    dPhiT_       = b.dPhiT_;
    stiffMatrix_ = b.stiffMatrix_;
    LegMatrix_   = b.LegMatrix_;
  }

}


//************************************************************************************
//************************************************************************************
// METHOD : resize
// DESC   : resizes dynamically allocated quantities
//          if required
// ARGS   : GINT newOrder
// RETURNS: TRUE on success, else FALSE
//************************************************************************************
template<typename T, typename TE>
GBOOL GLBasis<T,TE>::resize(GINT  newOrder)
{

  // No resizing necessary if already less than
  // previously allocated quantities:

  Np_ = newOrder;

  GINT  NN = Np_+1;

  //  resize xiNodes_:
  xiNodes_.resize(NN);

  //  resize weights_:
  weights_.resize(NN);

  //  resize Pn_, derivatives::
  Pn_.resize(NN);
  dPn_.resize(NN);


  //  resize basis Phi_:
  Phi_.resize(NN,NN);

  //  resize basis dPhi_:
  dPhi_.resize(NN,NN);
  dPhiT_.resize(NN,NN);

  //  resize stiffMatrix_:
  stiffMatrix_.resize(NN,NN);

  if ( !init() ) return FALSE;
  
  return TRUE;

} // end of method resize


//************************************************************************************
//************************************************************************************
// METHOD : computeNodes 
// DESC   : computes nodes and weights based on a Gauss-Legendre integration scheme
//          Note: nodes are computed from smalles to largest xi
//          Note: this method was taken largely from Canuto et al, 1987, Appendix C.
//          Note: Legendre polynomials of order Np_, Np_-1 and Np_-2 and their
//                derivatives are computed here for use here and in 
//                computing integrals. These polynomials are used in
//                forming basis functions, Phi_, and their derivatives,
//                dPhi_, evaluated at each node.
// ARGS   : none
// RETURNS: TRUE on success, else FALSE
//************************************************************************************
template<typename T, typename TE>
GBOOL GLBasis<T,TE>::computeNodes()
{
  if ( Np_ < 1 ) return FALSE;
 
  GINT  i, j, k, np1, nh;
  T     det, pnp, pnp1p, pnp1m, rp, rm, dth, cd, sd, cs, ss, 
        pn , pnp1, pdnp1, pdn, pnm1, pdnm1, poly, pder,
        pdnp1p, pnm1p,  pdnp1m, pnm, pdnm, pnm1m, pdnp,
        recsum, x, delx, cssave, error, a, b;
 
  alpha_ = 0.0;
  beta_  = 0.0;  //alpha_=beta_=0 ==> Legendre fcn basis 

//rv    = 1.0 + alpha_;
  np1   = Np_ + 1;
  nh    = np1 / 2;

  computeJacobi(np1,alpha_, beta_, pnp1p, pdnp1p, pnp, pdnp, pnm1p, pdnm1, ximax_);
  computeJacobi(np1,alpha_, beta_, pnp1m, pdnp1m, pnm, pdnm, pnm1m, pdnm1, ximin_);
  det  = pnp*pnm1m - pnm*pnm1p;
  rp   = -pnp1p;
  rm   = -pnp1m;
  a    = ( rp*pnm1m - rm*pnm1p ) / det;
  b    = ( rm*pnp   - rp*pnm   ) / det;

  // set up recursion relation for the initial guesses for roots:
  dth = 3.1415926535898/( 2.0*Np_ + 1.0 );
  cd  = cos(2.0*dth);
  sd  = sin(2.0*dth);
  cs  = cos    (dth);
  ss  = sin    (dth);

  // compute the roots by polynomial deflation:
  // together with Newton-Raphson....
  for ( j=0; j<nh; j++ ) {
    x = cs;
    error = 1.0;
    for ( k=0; k<kstop_ && error>eps_; k++ ) {
      computeJacobi(np1, alpha_, beta_, pnp1, pdnp1, pn, pdn, pnm1, pdnm1, x);
//    poly = pnp1 + a*pn + b*pnm1;
//    pder = pdnp1 + a*pdn + b*pdnm1;
      poly = pnp1 ;
      pder = pdnp1;
      recsum = 0.0;
      for ( i=0; i<j; i++ ) {
        recsum += (  1.0/(x - xiNodes_[i]) );
      }
      delx = -poly/(pder - recsum*poly);
      x += delx;
      error = fabs(delx)/ fabs(x);    // MIN(fabs(ximin_),fabs(ximax_));
    }
    xiNodes_[j] = x;
    cssave     = cs*cd - ss*sd;
    ss         = cs*sd + ss*cd;
    cs         = cssave;
  }


#if 1
  // Symmetrize the nodes:  NOTE: this is correct only for the cases
  // where (xmin,xmax) = (-y, y), or something like this.
  for ( i=0; i<nh; i++ ) {
    xiNodes_(np1-i-1) = -xiNodes_[i];
  }


  if ( Np_ == ( 2*(Np_/2) ) ) {
    xiNodes_(nh) = 0.0;
  }

  // re-order from smallest to largest:
  GQUAD txi;
  for ( i=0; i<nh; i++ ) {
    txi           = xiNodes_(Np_-i);
    xiNodes_(Np_-i) = xiNodes_[i];
    xiNodes_   [i] = txi;
  }
#endif

  return TRUE;

} // end of method computeNodes


//************************************************************************************
//************************************************************************************
// METHOD : computeWeights
// DESC   : 
//          NOTE: method not really intended to be called publicly; 
//          it should be called only when computeNodes is called. 
//          For other bases, this method will change, however.
// ARGS   : none
// RETURNS: TRUE on success, else FALSE
//************************************************************************************
template<typename T, typename TE>
GBOOL GLBasis<T,TE>::computeWeights()
{
  if ( Np_ < 1 ) return FALSE;
 
  GINT  i, np1=Np_+1;
  GQUAD fact;
  GQUAD ppn, pder, pm1, pdm1, pm2, pdm2;

  for ( i=0; i<np1; i++ ) {
    computeJacobi(np1, alpha_, beta_, ppn, pder,pm1, pdm1, pm2, pdm2, xiNodes_[i]);
    Pn_ [i] = ppn;    
    dPn_[i] = pder;    
    fact   = 2.0/( 1.0 - xiNodes_[i]*xiNodes_[i] );
    weights_[i] = (pder==0.0)?0.0:fact/(pder*pder); 
  }

  return TRUE;

} // end of method computeWeights


#if 0
//************************************************************************************
//************************************************************************************
// METHOD : computeBasisAtNodes
// DESC   :
// ARGS   : none
// RETURNS: TRUE on success, else FALSE
//************************************************************************************
template<typename T, typename TE>
GBOOL GLBasis<T,TE>::computeBasisAtNodes()
{
  if ( !bInit_ && !init() ) return NULLPTR;


  GINT  i, j, np1=Np_+1;
  GQUAD ppn_i, pder_i, pm1, pdm1, pm2, pdm2, ppn_j, pder_j;
  GQUAD fact=1.0/(Np_*(Np_.0)), gfact;

  Phi_ = 0.0;

  for ( i=0; i<np1; i++ ) {
    Phi_(i,i) = 1.0;
    for ( j=0; j< i; j++ ) {
      gfact = 1.0/(xiNodes_[j]-xiNodes_[i]);
      Phi_(i,j) = 1.0/dPn_[i] * gfact * Pn_[j];
    }
    for ( j=i+1; j<np1; j++ ) {
      gfact = 1.0/(xiNodes_[j]-xiNodes_[i]);
      Phi_(i,j) = 1.0/dPn_[i] * gfact * Pn_[j];
    }
  }

} // end of method computeBasisAtNodes
#endif

//************************************************************************************
//************************************************************************************
// METHOD : computeDerivMatrix
// DESC   :
//          NOTE: computeWeights_ (computeNodes) must have been called
//                prior to entering this method, s.t. the Pn__i have been
//                calculated. 
// ARGS   : none
// RETURNS: TRUE on success, else FALSE
//************************************************************************************
template<typename T, typename TE>
GBOOL GLBasis<T,TE>::computeDerivMatrix()
{
  if ( !bInit_ && !init() ) return FALSE;

  GINT  l, j, np1=Np_+1;
  GQUAD delxi, xi, xisq, d2Pn_;

  dPhi_ = 0.0;
  // Note: index j sweeps over basis number; l sweeps over node number

  for ( j=0; j<np1; j++ ) {
    for ( l=0; l<j; l++ ) {
      delxi      = xiNodes_(l) - xiNodes_[j];
      dPhi_ (l,j) =  dPn_(l)/(dPn_[j]*delxi);
    }
    for ( l=j+1; l<np1; l++ ) {
      delxi       = xiNodes_(l) - xiNodes_[j];
      dPhi_  (l,j) =  dPn_(l)/(dPn_[j]*delxi);
    }
  }
  for ( l=0; l<np1; l++ ) {
    xi          = xiNodes_(l);
    xisq        = xi * xi;
    d2Pn_        = xi*dPn_(l)/(1.0-xisq);
    dPhi_  (l,l) =  d2Pn_  / dPn_(l);
  }

  dPhi_.transpose(dPhiT_);

  return TRUE;

} // end of method computeDerivMatrix


//************************************************************************************
//************************************************************************************
// METHOD : computeLegendreMatrix
// DESC   : computes matrix M_ij = P_i (xi_j),
//          where P_i is the Legendre polynomial of order i, and
//          xi_j and W_j are the j-th nodal point, and weight,
//          respectively.
// ARGS   : none
// RETURNS: TRUE on success; else FALSE
//************************************************************************************
template<typename T, typename TE>
GBOOL GLBasis<T,TE>::computeLegendreMatrix()
{
  if ( !bNeedLegMatrix_ ) return TRUE;

  if ( !bInit_ && !init() ) return FALSE;
  
  GINT  i, j;
  GQUAD ppn_i, pder_i, pm1, pdm1, pm2, pdm2;

  for ( i=0; i<Np_+1; i++ ) {
    for ( j=0; j<Np_+1; j++ ) {
      computeJacobi(i, 0.0, 0.0, ppn_i , pder_i ,pm1, pdm1, pm2, pdm2, xiNodes_[j]);
//    LegMatrix_(i,j) =  ppn_i * weights_[j];
      LegMatrix_(i,j) =  ppn_i;
    }
  }
  bNeedLegMatrix_ = FALSE;
  return TRUE;

} // end of method computeLegendreMatrix


