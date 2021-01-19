//=================================================================================
// Module       : gllbasis.ipp
// Date         : 1/19/18 (DLR)
// Description  : Encapsulates the methods and data associated with
//                a spectral element nodal basis
//                Note: This basis is the Gauss-Lobatto-Legendre
//                basis, defined as 
//                  h(xi)_j = -1/(N(N+1)*L_N(xi_j))  * (1-xi^2)*dL_N(xi)/dxi / (xi-xi_j)
//                where N is the order of the expansion, and L_N is the Legendre
//                polynomial of order N, and the xi_j are the nodal points. L_N
//                is obtained from the generalized Jacobi polynomial computed in
//                computeJacobi.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : GNBasis
//==================================================================================
#include <cstdlib>
#include <memory>
#include <cmath>
#include <cstdio>
#include "gllbasis.hpp"


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method
// DESC   : Instantitate with explicit order and max order
// ARGS   : GINT inOrder
//          GINT maxORdern
// RETURNS: none
//**********************************************************************************
template<typename T, typename TE>
GLLBasis<T,TE>::GLLBasis(GINT inOrder, GINT  MaxOrder)
:
Np_                    (MIN(inOrder,MaxOrder)),
kstop_                 (128),
bInit_                 (FALSE),
bNeedDerivMatrix_      (TRUE),
bNeedLegMatrix_        (TRUE),
alpha_                 (0.0),
beta_                  (0.0),
ximin_                 (-1.0),
ximax_                 (1.0),
eps_                   (1.0e-8),
ttiny_                 (100.0*std::numeric_limits<T>::epsilon()),
tetiny_                (100.0*std::numeric_limits<TE>::epsilon())
{
  if ( Np_ < 1 ) {
    std::cout << "GLLBasis<T,TE>::GLLBasis: invalid expansion order Np_=" << Np_ << std::endl;
    exit(1);
  }
  resize(Np_);

} // end of constructor method


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method
// DESC   : Default constructor
// ARGS   : GINT inOrder
//          GINT maxORdern
// RETURNS: none
//**********************************************************************************
template<typename T, typename TE>
GLLBasis<T,TE>::GLLBasis()
:
Np_                    (0),
kstop_                 (128),
bInit_                 (FALSE),
bNeedDerivMatrix_      (TRUE),
bNeedLegMatrix_        (TRUE),
alpha_                 (0.0),
beta_                  (0.0),
ximin_                 (-1.0),
ximax_                 (1.0),
eps_                   (1.0e-8),
ttiny_                 (100.0*std::numeric_limits<T>::epsilon()),
tetiny_                (100.0*std::numeric_limits<TE>::epsilon())
{
} // end of constructor method


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method
// DESC   : Instantiate with explicit order, taken to be max order as well
// ARGS   : GINT inOrder
// RETURNS: none
//**********************************************************************************
template<typename T, typename TE>
GLLBasis<T,TE>::GLLBasis(GINT Order)
:
Np_                    (Order),
kstop_                 (128),
bInit_                 (FALSE),
bNeedDerivMatrix_      (TRUE),
bNeedLegMatrix_        (TRUE),
alpha_                 (0.0),
beta_                  (0.0),
ximin_                 (-1.0),
ximax_                 (1.0),
eps_                   (1.0e-8),
ttiny_                 (100.0*std::numeric_limits<T>::epsilon()),
tetiny_                (100.0*std::numeric_limits<TE>::epsilon())
{
  if ( Np_ < 1 ) {
    std::cout << "GLLBasis<T,TE>::GLLBasis: invalid expansion order Np_=" << Np_ << std::endl;
    exit(1);
  }
  resize(Np_);
  init();

} // end of constructor method


//**********************************************************************************
//**********************************************************************************
// METHOD : Copy Constructor method
// DESC   : 
// ARGS   : GLLBasis &
// RETURNS: none
//**********************************************************************************
template<typename T, typename TE>
GLLBasis<T,TE>::GLLBasis(const GLLBasis &b)
{
   // copy data:
    alpha_    = b.alpha_;
    beta_     = b.beta_;
    ximin_    = b.ximin_;
    ximax_    = b.ximax_;
    eps_      = b.eps_;
    ttiny_    = b.ttiny_;
    tetiny_   = b.tetiny_;
    Np_       = b.Np_;
    kstop_    = b.kstop_;
    bInit_    = b.bInit_;
    bNeedDerivMatrix_ = b.bNeedDerivMatrix_;
    bNeedLegMatrix_   = b.bNeedLegMatrix_;


    //  matrices, and node data:
    xiNodes_       = b.xiNodes_;
    weights_       = b.weights_;
    Pn_            = b.Pn_;
    dPn_           = b.dPn_;
    Phi_           = b.Phi_;
    dPhi_          = b.dPhi_;
    dPhiT_         = b.dPhiT_;
    stiffMatrix_   = b.stiffMatrix_;
    LegMatrix_     = b.LegMatrix_;
    LegTransform_  = b.LegTransform_;
    iLegTransform_ = b.iLegTransform_;
    LegFilterMat_  = b.LegFilterMat_;
    LegFilterMatT_ = b.LegFilterMatT_;

}

//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename T, typename TE>
GLLBasis<T,TE>::~GLLBasis()
{
}

//**********************************************************************************
//**********************************************************************************
// METHOD : Assignment method
// DESC   : 
// ARGS   : GLLBasis &
// RETURNS: none
//**********************************************************************************
template<typename T, typename TE>
void GLLBasis<T,TE>::operator=(const GLLBasis &b)
{
  if ( &b != this ) 
  {
   // copy data:
    ximin_        = b.ximin_;
    ximax_        = b.ximax_;
    eps_          = b.eps_;
    Np_           = b.Np_;
    kstop_        = b.kstop_;
    bInit_        = b.bInit_;
    bNeedDerivMatrix_ = b.bNeedDerivMatrix_;
    bNeedLegMatrix_  = b.bNeedLegMatrix_;


    //  matrices, and node data:
    xiNodes_      = b.xiNodes_;
    weights_      = b.weights_;
    Pn_           = b.Pn_;
    dPn_          = b.dPn_;
    dPhi_         = b.dPhi_;
    dPhiT_        = b.dPhiT_;
    stiffMatrix_  = b.stiffMatrix_;
    LegMatrix_    = b.LegMatrix_;
    LegTransform_ = b.LegTransform_;
    iLegTransform_= b.iLegTransform_;
    LegFilterMat_ = b.LegFilterMat_;
    LegFilterMatT_= b.LegFilterMatT_;
  }

}

//************************************************************************************
//************************************************************************************
// METHOD : getXimin
// DESC   : Get minimim of reference interval 
// ARGS   : none
// RETURNS: element left boundary
//************************************************************************************
template<typename T, typename TE>
T GLLBasis<T,TE>::getXimin()
{
  return (T) ximin_;
} // end of method getXimin


//************************************************************************************
//************************************************************************************
// METHOD : getXimax
// DESC   : Get maximum of reference interval 
// ARG    : none
// RETURNS: element right boundary
//************************************************************************************
template<typename T, typename TE>
T GLLBasis<T,TE>::getXimax()
{
  return (T) ximax_;
} // end of method getXimax


//************************************************************************************
//************************************************************************************
// METHOD : getOrder
// DESC   : Get Lag. interp. polynomial order  
// ARG    : none
// RETURNS: GINT  element expansion order
//************************************************************************************
template<typename T, typename TE>
GINT  GLLBasis<T,TE>::getOrder()
{
  return Np_;
} // end of method getOrder


#if 0

//************************************************************************************
//************************************************************************************
// METHOD : getBasisAtNodes
// DESC   : Computes Gauss-Lobatto Legendre basis
//              as a function of  nodes in the parent domain. 
//              The returned quantity is a pointer to a matrix
//              Phi__i(xi_j)==Phi_(i,j),
//              where xi_j is the jth node element, and
//              Phi__i is given by (see Canuto, et. al.):
//
//              Phi__j (x) = (N*(N+1)Pn(xi_j))^-1  * (1-x^2)*dPn/dx(x)/(x - xi_j)
//              
// ARGS   : ret: array to T      
// RETURNS: pointer to ret data on success; else ; else NULLPTR 
//************************************************************************************
template<typename T, typename TE>
GBasisMatrix *GLLBasis<T,TE>::getBasisAtNodes(GBasisMatrix &ret)
{
  if ( !bInit_ && !init() ) return NULLPTR;
 
  if (!computeBasisAtNodes() ) return NULLPTR;


  if ( ret.size(1) < Phi_.size(1) || ret.size(2) < Phi_.size(2) ) return NULLPTR;

  GINT  i;
  T    *fptr=ret.data();
  T    *qptr=Phi_.data();
   
  for ( i=0; i<Phi_.size(1) * Phi_.size(2); i++ )
    *(fptr+i) = (T      )(*(qprt+i));

  return &ret;
} // end of method getBasisAtNodes


//************************************************************************************
//************************************************************************************
// METHOD : setXiDomain
// DESC   : Set reference domain 
// ARGS   : min, max 1d domain extent
// RETURNS: none
//************************************************************************************
template<typename T, typename TE>
void GLLBasis<T,TE>::setXiDomain(T min, T max)
{
  ximin_ = MIN(min,max);
  ximax_ = MAX(min,max);

} // end of method setXiDomain
#endif


//************************************************************************************
//************************************************************************************
// METHOD : SetOrder
// DESC   : 
// ARGS   : order: interp. polynomial order
// RETURNS: none
//************************************************************************************
template<typename T, typename TE>
void GLLBasis<T,TE>::setOrder(GINT order)
{
  if ( order != Np_ )
  {
     resize(order);
     if ( init() ) {
       std::cout << "GLLBasis<T,TE>::setOrder: Initialization failure" << std::endl;
       exit(1);
     }
     bNeedDerivMatrix_ = TRUE;
     bNeedLegMatrix_ = TRUE;
  }
  Np_ = order;
} // end of method setOrder


//************************************************************************************
//************************************************************************************
// METHOD : resize
// DESC   : resizes dynamically allocated quantities
//          if required
// ARGS   : newOrder: new interp. polynomial order
// RETURNS: TRUE on success, else FALSE
//************************************************************************************
template<typename T, typename TE>
GBOOL GLLBasis<T,TE>::resize(GINT  newOrder)
{

  // No resizing necessary if already less than
  // previously allocated quantities:

  Np_ = newOrder;

  //  resize xiNodes_:
  xiNodes_ .resize(Np_+1);
  xiNodesEv_ .resize(Np_+1);

  //  resize weights_:
  weights_.resize(Np_+1);
  weightsEv_.resize(Np_+1);
  iweightsEv_.resize(Np_+1);

  //  resize Pn_:
  Pn_.resize(Np_+1);

  //  resize dPn_:
  dPn_.resize(Np_+1);

  //  resize basis Phi_:
  Phi_.resize(Np_+1,Np_+1);

  //  resize basis dPhi_:
  dPhi_  .resize(Np_+1,Np_+1);
  dPhiT_ .resize(Np_+1,Np_+1);
  dPhiEv_ .resize(Np_+1,Np_+1);
  dPhiTEv_.resize(Np_+1,Np_+1);
  dPhiWEv_ .resize(Np_+1,Np_+1);
  dPhiWTEv_ .resize(Np_+1,Np_+1);
  dPhiiWEv_ .resize(Np_+1,Np_+1);
  dPhiiWTEv_ .resize(Np_+1,Np_+1);

  //  resize stiffMatrix_:
  stiffMatrix_  .resize(Np_+1,Np_+1);
  stiffMatrixEv_.resize(Np_+1,Np_+1);

  //  resize LegMatrix_:
  LegMatrix_    .resize(Np_+1,Np_+1);
  LegTransform_ .resize(Np_+1,Np_+1);
  iLegTransform_.resize(Np_+1,Np_+1);
  LegFilterMat_ .resize(Np_+1,Np_+1);
  LegFilterMatT_.resize(Np_+1,Np_+1);

  if ( !init() ) return FALSE; 

  return TRUE;

} // end of method resize


//************************************************************************************
//************************************************************************************
// METHOD : computeNodes 
// DESC   : Computes nodes and weights based on a Gauss-Lobatto Legendre integration scheme
//              Note: nodes are computed from smalles to largest xi
//              Note: this method was taken largely from Canuto et al, 1987, Appendix C.
//              Note: this method will be generalized to accept this as a base class,
//                    and derive other basis types from it....
//              Note: Legendre polynomials of order Np_, Np_-1 and Np_-2 and their
//                    derivatives are computed here for use here and in 
//                    computing integrals. These polynomials are used in
//                    forming basis functions, Phi_, and their derivatives,
//                    dPhi_, evaluated at each node.
// ARGS   : none
// RETURNS: TRUE on success, else FALSE
//************************************************************************************
template<typename T, typename TE>
GBOOL GLLBasis<T,TE>::computeNodes()
{
  GString serr = "GLLBasis<T,TE>::computeNodes: ";
  GBOOL bret;
  GINT  i, j, k, nh, np1;
  
  T     det, pnp, pnp1p, pnp1m, rp, rm, dth, cd, sd, cs, ss, 
        pn , pnp1, pdnp1, pdn, pnm1, pdnm1, a, b, poly, pder,
        pdnp1p, pnm1p,  pdnp1m, pnm, pdnm, pnm1m, pdnp,
        recsum, x, delx, cssave, error;


  if ( Np_ < 1 ) return FALSE;
  if ( Np_ == 1 ) {
    xiNodes_[0] = ximin_;
    xiNodes_[1] = ximax_;
    if ( computeWeights() ) bInit_ = TRUE;
    else                    return FALSE;
  }
 
  alpha_ = 0.0;
  beta_  = 0.0;  //alpha_=beta_=0 ==> Legendre fcn basis 

//rv    = 1.0 + alpha_;
  np1   = Np_ + 1;

  computeJacobi(np1,alpha_, beta_, pnp1p, pdnp1p, pnp, pdnp, pnm1p, pdnm1, ximax_);
  computeJacobi(np1,alpha_, beta_, pnp1m, pdnp1m, pnm, pdnm, pnm1m, pdnm1, ximin_);
  det  = pnp*pnm1m - pnm*pnm1p;
  rp   = -pnp1p;
  rm   = -pnp1m;
  a    = ( rp*pnm1m - rm*pnm1p ) / det;
  b    = ( rm*pnp   - rp*pnm   ) / det;

  // order nodes from largest to smallest:
  xiNodes_ [0] = ximax_;
  xiNodes_[Np_] = ximin_;
  nh = ( Np_ + 1 ) / 2;

  // set up recursion relation for the initial guesses for roots:
  dth = 3.1415926535898/( 2.0*Np_ + 1.0 );
  cd  = cos(2.0*dth);
  sd  = sin(2.0*dth);
  cs  = cos    (dth);
  ss  = sin    (dth);

  // compute the first half of the roots by polynomial deflation:
  for ( j=1; j<nh+1; j++ ) {
    x = cs;
    for ( k=0; k<kstop_; k++ ) {
      computeJacobi(np1, alpha_, beta_, pnp1, pdnp1, pn, pdn, pnm1, pdnm1, x);
      poly = pnp1 + a*pn + b*pnm1;
      pder = pdnp1 + a*pdn + b*pdnm1;
      recsum = 0.0;
      for ( i=0; i<j; i++ ) {
        recsum += (  1.0/(x - xiNodes_[i]) );
      }
      delx = -poly/(pder - recsum*poly);
      x += delx;
      error = fabs(delx)/ fabs(x);    // MIN(fabs(ximin_),fabs(ximax_));
      if ( error < eps_ ) break;
    }
    xiNodes_[j] = x;
    cssave     = cs*cd - ss*sd;
    ss         = cs*sd + ss*cd;
    cs         = cssave;
  }

  
  // Symmetrize the nodes:  NOTE: this is correct only for the cases
  // where (xmin,xmax) = (-y, y); ie, is symmetric interval.
  for ( i=0; i<nh; i++ ) {
    xiNodes_[np1-i-1] = -xiNodes_[i];
  }

  if ( Np_ == ( 2*(Np_/2) ) ) {
    xiNodes_[nh] = 0.0;
  }


  // re-order from smallest to largest:
  T txi;
  for ( i=0; i<nh; i++ ) {
    txi           = xiNodes_[Np_-i];
    xiNodes_[Np_-i] = xiNodes_[i];
    xiNodes_   [i] = txi;
  }

  return TRUE;

} // end of method computeNodes


//************************************************************************************
//************************************************************************************
// METHOD : computeWeights
// DESC   : 
//            NOTE: method not really intended to be called publicly; 
//            it should be called only when computeNodes is called. 
//            For other bases, this method will change, however.
// ARGS   : none
// RETURNS: TRUE on success, else FALSE
//************************************************************************************
template<typename T, typename TE>
GBOOL GLLBasis<T,TE>::computeWeights()
{
  if ( Np_ < 1 ) return FALSE;
 
  GINT  i;
  T     fact = 2.0/(Np_*(Np_ + 1.0));
  T     ppn, pder, pm1, pdm1, pm2, pdm2;

  for ( i=0; i<Np_+1; i++ ) {
    computeJacobi(Np_, alpha_, beta_, ppn, pder,pm1, pdm1, pm2, pdm2, xiNodes_[i]);
    Pn_ [i] = ppn;    
    dPn_[i] = pder;    
    weights_[i] = (Pn_[i]==0.0)?0.0:fact/(Pn_[i]*Pn_[i]); // Note: Pn_ computed in ComputNodes
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
GBOOL GLLBasis<T,TE>::computeBasisAtNodes()
{
  if ( !bInit_ && !init() ) return FALSE;


  GINT  i, j;
  T     ppn_i, pder_i, pm1, pdm1, pm2, pdm2, ppn_j, pder_j;
  T     fact=1.0/(Np_*(Np_+1.0)), gfact;

  Phi_ = 0.0;

  for ( i=0; i< Np_+1; i++ ) {
    computeJacobi(Np_, alpha_, beta_, ppn_i, pder_i,pm1, pdm1, pm2, pdm2, xiNodes_[i]);
    for ( j=0; j< Np_+1; j++ ) {
      Phi_(i,j) = 1.0;
      if ( i != j ) {
        computeJacobi(Np_, alpha_, beta_, ppn_j, pder_j,pm1, pdm1, pm2, pdm2, xiNodes_[j]);
        gfact = (1-xiNodes_[j]*xiNodes_[j])/(xiNodes_[j]-xiNodes_[i]);
        Phi_(i,j) = fact/ppn_i * gfact * pder_j;
      }
    }
  }

} // end of method computeBasisAtNodes
#endif

//************************************************************************************
//************************************************************************************
// METHOD : computeDerivMatrix
// DESC   :
//            NOTE: computeWeights (computeNodes) must have been called
//                  prior to entering this method, s.t. the Pn_i have been
//                  calculated. 
// ARGS   : none
// RETURNS: TRUE on success, else FALSE
//************************************************************************************
template<typename T, typename TE>
GBOOL GLLBasis<T,TE>::computeDerivMatrix()
{

  GINT  l, j;
  T     delxi;

  dPhi_ = 0.0;
  // Note: index j sweeps over basis number; l sweeps over node number

  for ( j=0; j<Np_+1; j++ ) {
    for ( l=0; l<j; l++ ) {
      delxi      = xiNodes_[l] - xiNodes_[j];
      dPhi_ (l,j) =  Pn_[l]/(Pn_[j]*delxi+ttiny_);
    }
    for ( l=j+1; l<Np_+1; l++ ) {
      delxi       = xiNodes_[l] - xiNodes_[j];
      dPhi_  (l,j) =  Pn_[l]/(Pn_[j]*delxi+ttiny_);
    }
  }
  dPhi_ (0 ,0 ) = -0.25*Np_*(Np_ + 1.0);
  dPhi_ (Np_,Np_) = -dPhi_(0,0);

  dPhi_.transpose(dPhiT_);
  bNeedDerivMatrix_ = FALSE;

  return TRUE;

} // end of method computeDerivMatrix


//************************************************************************************
//************************************************************************************
// METHOD : computeJacobi
// DESC   : Compute Jacobi polynomial nodes and derivatives for polynomial
//          type specified by alpha_, beta_. Taken from: 
// ARGS   : 
// RETURNS: none  
//************************************************************************************
template<typename T, typename TE>
void GLLBasis<T,TE>::computeJacobi(GINT  &N, T alpha_ , T beta_  , 
                                       T &poly  , T &pder  , T &polym1, 
                                       T &pderm1, T &polym2, T &pderm2, T &x)
{

  GINT  k; 
  T     apb  , polylist, pderlist, rv                ;
  T     polyn, pdern   , psave   , pdsave  , fk      ;
  T     a1   , a2      , a3      , a4      , b3      ;


  apb = alpha_ + beta_;
  rv = 1.0 + alpha_;

  poly = 1.0;
  pder = 0.0;

  if ( N == 0 ) return;

  polylist = poly;
  pderlist = pder;
  poly     = rv * x;
  pder     = rv;
  
  if ( N == 1 ) return;

  for ( k=2; k<=N; k++ ) {
    fk = static_cast<T>(k);
    a1 = 2.0*fk*(fk+apb)*(2.0*fk+apb-2.0);
    a2 = (2.0*fk+apb-1.0)*(alpha_*alpha_ -beta_*beta_);
    b3 = 2.0*fk+apb-2.0;
    a3 = b3*(b3+1.0)*(b3+2.0);
    a4 = 2.0*(fk+alpha_-1.0)*(fk+beta_-1.0)*(2.0*fk+apb);
    polyn    = ((a2+a3*x)*poly - a4*polylist)/a1;
    pdern    = ((a2+a3*x)*pder - a4*pderlist+a3*poly)/a1;
    psave    = polylist;
    pdsave   = pderlist;
    polylist = poly;
    poly     = polyn;
    pderlist = pder;
    pder     = pdern;
  }

  polym1 = polylist;
  pderm1 = pderlist;
  polym2 = psave;
  pderm2 = pdsave; 

} // end of method computeJacobi


//************************************************************************************
//************************************************************************************
// METHOD : computeStiffMatrix
// DESC   : computes stiffness matrix for this basis.
//          Method uses Gaussian integration, so the
//          nodes and the weights for the current Np_ order
//          must be calculated.
//
//          The stiffness matrix is defined as:
//              
//          A_(i,j )= Integral(ximin_,ximax_) { dPhi__i(xi)/dxi dPhi__j(xi)/dxi dxi }
//                  = Sum { w_k * dPhi__i(xi_k)/dxi dPhi__j(xi_k)/dxi }   (Gaussian quadrature)
//          where w_k are the weights for the basis associated with each node, xi_k.
// ARGS   : none
// RETURNS: TRUE on success, else FALSE
//************************************************************************************
template<typename T, typename TE>
GBOOL GLLBasis<T,TE>::computeStiffMatrix()
{

  GINT  i, j, k;

  stiffMatrix_ = 0.0;

  // Could use matrix algebra here....
  for ( i=0; i<Np_+1; i++ ) {
    for ( j=0; j<Np_+1; j++ ) {
       for ( k=0; k<Np_+1; k++ ) {
         stiffMatrix_(i,j) += weights_[k]*dPhi_(k,i)*dPhi_(k,j) ;  
       }
    } 
  } 
  return TRUE;
} // end of method computeStiffMatrix


//************************************************************************************
//************************************************************************************
// METHOD : init
// DESC   : Computes all quantities, nodes, weights, mass and
//          stiffness matrices, and derivative matrices Phi_, and dPhi_.
// ARGS   : none
// RETURNS: TRUE on success, else FALSE
//************************************************************************************
template<typename T, typename TE>
GBOOL GLLBasis<T,TE>::init()
{
  if ( !computeNodes() ) return FALSE;

  if ( !computeWeights() ) return FALSE;

  if ( !computeDerivMatrix() ) return FALSE;

  if ( !computeStiffMatrix() ) return FALSE; // stiffness matrix computed; dPhi_ also computed.

  // Copy computated data to the 'evaluated' structures:
  getXiNodes(xiNodesEv_);
  getWeights(weightsEv_);
  getiWeights(iweightsEv_); // computes inverse and cast
  getStiffMatrix(stiffMatrixEv_);
  getDerivMatrix(dPhiEv_,FALSE);     //  D
  getDerivMatrix(dPhiTEv_,TRUE);     //  D^T
  getDerivMatrixW(dPhiWEv_,FALSE);   //  Diag(W)*D
  getDerivMatrixW(dPhiWTEv_,TRUE);   // (Diag(W)*D)^T
  getDerivMatrixiW(dPhiiWEv_,FALSE); //  Diag(W^-1)*D
  getDerivMatrixiW(dPhiiWTEv_,TRUE); // (Diag(W^-1)*D)^T

  bInit_ = TRUE;

#if 0
  // Note: following call must have bInit_ = TRUE
  if ( !computeLegTransform(2) ) return FALSE; // Legendre transform matrix
#endif

  return TRUE;

} // end of method init


//************************************************************************************
//************************************************************************************
// METHOD : computeLegendreMatrix
// DESC   : Computes matrix M_ij = P_i (xi_j),
//          where P_i is the Legendre polynomial of order i, and
//          xi_j is the j-th nodal point
// ARGS   : none.
// RETURNS: TRUE on success; else FALSE 
//************************************************************************************
template<typename T, typename TE>
GBOOL GLLBasis<T,TE>::computeLegendreMatrix()
{

  GINT  i, j, p;
  T     ppn_i, pder_i, pm1, pdm1, pm2, pdm2;

  for ( i=0; i<Np_+1; i++ ) {
    for ( j=0; j<Np_+1; j++ ) {
      p = i;
      computeJacobi(p, 0.0, 0.0, ppn_i , pder_i ,pm1, pdm1, pm2, pdm2, xiNodes_[j]);
//    LegMatrix_(i,j) = ppn_i * weights_[j];
      LegMatrix_(i,j) = ppn_i ;
    }
  }

  return TRUE;

} // end of method computeLegendreMatrix


//************************************************************************************
//************************************************************************************
// METHOD : computeLegTransform
// DESC   : Computes Legendre transform matrix that
//          enables conversion to modal space. 
// ARGS   : ifilter: reference mode number 
// RETURNS: TRUE on success; else FALSE 
//************************************************************************************
template<typename T, typename TE>
GBOOL GLLBasis<T,TE>::computeLegTransform(GINT ifilter)
{

  GINT  i, j;
  
  assert(ifilter >= 0 && ifilter < Np_+1);

  for ( i=0; i<Np_+1; i++ ) {
    for ( j=0; j<Np_+1; j++ ) {
      if ( j < ifilter ) {
        LegTransform_(i,j) = evalBasis(j,xiNodes_[i]);
      }
      else {
        LegTransform_(i,j) = evalBasis(j,xiNodes_[i]) - evalBasis(j-ifilter,xiNodes_[i]);
      }
    }
  }
  LegTransform_.inverse(iLegTransform_);  


  return TRUE;

} // end of method computeLegTransform


//************************************************************************************
//************************************************************************************
// METHOD : getXiNodesComp
// DESC   : Get list (T vector) of reference interval _computational_ nodes  
// ARGS   : none
// RETURNS: GTVector*
//************************************************************************************
template<typename T, typename TE>
GTVector<T> *GLLBasis<T,TE>::getXiNodesComp()
{
  return &xiNodes_;
} // end of method getXiNodesComp


//************************************************************************************
//************************************************************************************
// METHOD : getXiNodes
// DESC   : Get list (TE vector) of reference interval _evaluation_ nodes
// ARGS   : none
// RETURNS: GTVector*
//************************************************************************************
template<typename T, typename TE>
GTVector<TE> *GLLBasis<T,TE>::getXiNodes()
{
  return &xiNodesEv_;
} // end of method getXiNodes


//************************************************************************************
//************************************************************************************
// METHOD : getXiNodes (2)
// DESC   : Get deep copy (TE array *) of reference interval nodes
// ARGS   : ret : array of TE type
//          num : num elements in ret array
// RETURNS: TE *
//************************************************************************************
template<typename T, typename TE>
 TE *GLLBasis<T,TE>::getXiNodes(TE *ret, GINT  num)
{

  if ( ret == NULLPTR || num < xiNodes_.size() ) return NULLPTR;

  GINT  i;
  T    *qptr=xiNodes_.data();

  for ( i=0; i<xiNodes_.size(); i++ )
    ret[i] = static_cast<TE>(*(qptr+i));

  return ret;
} // end of method getXiNodes (2)


//************************************************************************************
//************************************************************************************
// METHOD : getXiNodes (3)
// DESC   : Get deep copy of reference nodes, returning in specified array GTVector
// ARGS   : GTVector &ret
// RETURNS: pointer to input vector on success; else NULLPTR
//************************************************************************************
template<typename T, typename TE>
 void GLLBasis<T,TE>::getXiNodes(GTVector<TE> &ret)
{
  getXiNodes(ret.data(), ret.size());
} // end of method getXiNodes (3)


//************************************************************************************
//************************************************************************************
// METHOD : getWeightsComp
// DESC   : Get GTVector<T> member data vector weighters 
// ARGS   : none
// RETURNS: pointer to GTVector member data
//************************************************************************************
template<typename T, typename TE>
GTVector<T> *GLLBasis<T,TE>::getWeightsComp()
{
  return &weights_;
} // end of method getWeightsComp


//************************************************************************************
//************************************************************************************
// METHOD : getWeights
// DESC   : Get GTVector<TE> member data vector _evaluation_ weights
// ARGS   : none
// RETURNS: pointer to GTVector member data
//************************************************************************************
template<typename T, typename TE>
GTVector<TE> *GLLBasis<T,TE>::getWeights()
{
  return &weightsEv_;
} // end of method getWeights 


//************************************************************************************
//************************************************************************************
// METHOD : getWeights (2)
// DESC   : Get deep copy of weights in specified return array
// ARGS   : ret : TE type array
//          num : num elements in ret array
// RETURNS: TE *, if successful; else NULLPTR
//************************************************************************************
template<typename T, typename TE>
TE *GLLBasis<T,TE>::getWeights(TE *ret, GINT  num)
{
  GString serr = "GLLBasis<T,TE>::getWeights(2): ";
  if ( num < weights_.size() ) {
    std::cout << serr << "Incompatible vector length" << std::endl;
    exit(1);
  }


  for ( GINT  i=0; i<weights_.size(); i++ )
    ret[i] = static_cast<TE>(weights_[i]);

  return ret;
} // end of method getWeights (2)


//************************************************************************************
//************************************************************************************
// METHOD : getWeights (3)
// DESC   : Get deep copy of weights in specified return object
// ARGS   : ret : GTVector<TE> array
// RETURNS: pointer to ret GTVector
//************************************************************************************
template<typename T, typename TE>
void GLLBasis<T,TE>::getWeights(GTVector<TE> &ret)
{
  getWeights(ret.data(), ret.size());

} // end of method getWeights (3)


//************************************************************************************
//************************************************************************************
// METHOD : getiWeights (1)
// DESC   : Get GTVector<TE> member data pointer to  vector _evaluation_ 
//          inverse weights 
// ARGS   : none
// RETURNS: pointer to GTVector member data
//************************************************************************************
template<typename T, typename TE>
GTVector<TE> *GLLBasis<T,TE>::getiWeights()
{
  return &iweightsEv_;
} // end of method getiWeights (1) 


//************************************************************************************
//************************************************************************************
// METHOD : getiWeights (2)
// DESC   : Get deep copy of inverse weights in specified return object
// ARGS   : ret : GTVector<TE> array
// RETURNS: pointer to ret GTVector
//************************************************************************************
template<typename T, typename TE>
void GLLBasis<T,TE>::getiWeights(GTVector<TE> &ret)
{

  for ( GINT  i=0; i<weights_.size(); i++ )
    ret[i] = 1.0 / static_cast<TE>(weights_[i]);

} // end of method getiWeights (2)


//************************************************************************************
//************************************************************************************
// METHOD : getStiffMatrixComp
// DESC   : Get pointer to member data to comptation 1d stiffness GTMatrix<T> 
// ARGS   : none
// RETURNS: GTMatrix *
//************************************************************************************
template<typename T, typename TE>
GTMatrix<T> *GLLBasis<T,TE>::getStiffMatrixComp()
{
  return &stiffMatrix_;
} // end of method getStiffMatrixComp 


//************************************************************************************
//************************************************************************************
// METHOD : getStiffMatrix (1)
// DESC   : Get pointer to member data 1d stiffness GTMatrix<T      > 
// ARGS   : none
// RETURNS: GTMatrix *
//************************************************************************************
template<typename T, typename TE>
GTMatrix<TE> *GLLBasis<T,TE>::getStiffMatrix()
{
  return &stiffMatrixEv_;
} // end of method getStiffMatrix (1)


//************************************************************************************
//************************************************************************************
// METHOD : getStiffMatrix (2)
// DESC   : Deep copy 1d stiffness matrix to supplied object
// ARGS   : ret: GTMatrix<TE> 
// RETURNS: none.
//************************************************************************************
template<typename T, typename TE>
void GLLBasis<T,TE>::getStiffMatrix(GTMatrix<TE> &ret)
{
  GString serr = "GLLBasis<T,TE>::getStiffMatrix(2): ";
  if ( ret.size(1) < stiffMatrix_.size(1) || ret.size(2) < stiffMatrix_.size(2) ) {
    std::cout << serr << "Incompatible matrix sizes" << std::endl;
    exit(1);
  }

  for ( GINT i=0; i<stiffMatrix_.size(1); i++ ) 
    for ( GINT  j=0; j<stiffMatrix_.size(2); j++ ) 
      ret(i,j) = static_cast<TE>(stiffMatrix_(i,j));

} // end of method getStiffMatrix (2)


//************************************************************************************
//************************************************************************************
// METHOD : getLegMatrix 
// DESC   : Deep copy of Legendre poly. matrix: Lp = Lp(nLegOrder,iNodeIndex)
// ARGS   : ret: GTMatrix<TE> 
// RETURNS: none.
//************************************************************************************
template<typename T, typename TE>
void GLLBasis<T,TE>::getLegMatrix(GTMatrix<TE> &ret)
{
  if ( ret.size(1) < LegMatrix_.size(1) || ret.size(2) < LegMatrix_.size(2) ) return;
  if ( bNeedLegMatrix_ && !computeLegendreMatrix() ) return;

  for ( GINT  i=0; i<LegMatrix_.size(1); i++ )
    for ( GINT  j=0; j<LegMatrix_.size(2); j++ )
      ret(i,j) = static_cast<TE>(LegMatrix_(i,j));

} // end of method getLegMatrix


//************************************************************************************
//************************************************************************************
// METHOD : getFilterMat
// DESC   : Get pointer to Legendre transformation filter matrix. This
//          isn't filled here, but may be set from caller, and stored
//          here
// ARGS   : btranspose: return transpose?
// RETURNS: GTMatrix *
//************************************************************************************
template<typename T, typename TE>
GTMatrix<TE> *GLLBasis<T,TE>::getFilterMat(GBOOL btranspose)
{
  if ( btranspose ) return &LegFilterMatT_;
  else              return &LegFilterMat_;

} // end of method getFilterMat


//************************************************************************************
//************************************************************************************
// METHOD : getLegTransform
// DESC   : Get pointer to Legendre transformation matrix
// ARGS   : none
// RETURNS: GTMatrix *
//************************************************************************************
template<typename T, typename TE>
GTMatrix<TE> *GLLBasis<T,TE>::getLegTransform()
{
  return &LegTransform_;

} // end of method getLegTransform


//************************************************************************************
//************************************************************************************
// METHOD : getiLegTransform
// DESC   : Get pointer to inverse of Legendre transformation matrix
// ARGS   : none
// RETURNS: GTMatrix *
//************************************************************************************
template<typename T, typename TE>
GTMatrix<TE> *GLLBasis<T,TE>::getiLegTransform()
{
  return &iLegTransform_;

} // end of method getiLegTransform


//************************************************************************************
//************************************************************************************
// METHOD : getDerivMatrixComp
// DESC   : Get derivative matrix _compujtational_ member data
// ARGS   : bTranspose: TRUE==>return transpose; else don't
// RETURNS: pointer to member data GTMatrix
//************************************************************************************
template<typename T, typename TE>
GTMatrix<T> *GLLBasis<T,TE>::getDerivMatrixComp(GBOOL bTranspose)
{
  if ( bTranspose ) return &dPhiT_;
  else              return &dPhi_;
} // end of method getDerivMatrixComp


//************************************************************************************
//************************************************************************************
// METHOD : getDerivMatrix (1)
// DESC   : Get derivative matrix member data
// ARGS   : bTranspose: TRUE==>return transpose; else don't
// RETURNS: pointer to member data GTMatrix
//************************************************************************************
template<typename T, typename TE>
GTMatrix<TE> *GLLBasis<T,TE>::getDerivMatrix(GBOOL bTranspose)
{
  if ( bTranspose ) return &dPhiTEv_;
  else              return &dPhiEv_;
} // end of method getDerivMatrix


//************************************************************************************
//************************************************************************************
// METHOD : getDerivMatrix (2)
// DESC   : Get deep copy of deriv matrix in specified return object
// ARGS   : ret: GTMatrix to return data
//          bTranspose: flag to get transpose (TRUE); else, not
// RETURNS: pointer to ret GTMatrix on success; else NULLPTR
//************************************************************************************
template<typename T, typename TE>
void GLLBasis<T,TE>::getDerivMatrix(GTMatrix<TE> &ret, GBOOL bTranspose)
{
  if ( ret.size(1) < dPhi_.size(1) || ret.size(2) < dPhi_.size(2) ) return ;

  if ( bTranspose ) {
    for ( GINT i=0; i<dPhi_.size(1); i++ ) 
      for ( GINT j=0; j<dPhi_.size(2); j++ ) 
        ret(i,j) = static_cast<TE>(dPhiT_(i,j));
  } 
  else {
    for ( GINT i=0; i<dPhi_.size(1); i++ ) 
      for ( GINT j=0; j<dPhi_.size(2); j++ ) 
        ret(i,j) = static_cast<TE>(dPhi_(i,j));
  }
  
} // end of method getDerivMatrix (2)


//************************************************************************************
//************************************************************************************
// METHOD : getDerivMatrixW (1)
// DESC   : Get W X derivative matrix member data
// ARGS   : bTranspose: TRUE==>return transpose; else don't
// RETURNS: pointer to member data GTMatrix
//************************************************************************************
template<typename T, typename TE>
GTMatrix<TE> *GLLBasis<T,TE>::getDerivMatrixW(GBOOL bTranspose)
{
  if ( bTranspose ) return &dPhiWTEv_;
  else              return &dPhiWEv_;
} // end of method getDerivMatrixW (1)



//************************************************************************************
//************************************************************************************
// METHOD : getDerivMatrixW (2)
// DESC   : Get deep copy of W X deriv matrix in specified return object
// ARGS   : ret: GTMatrix to return data
//          bTranspose: flag to get transpose (TRUE); else, not
// RETURNS: pointer to ret GTMatrix on success; else NULLPTR
//************************************************************************************
template<typename T, typename TE>
void GLLBasis<T,TE>::getDerivMatrixW(GTMatrix<TE> &ret, GBOOL bTranspose)
{
  if ( ret.size(1) < dPhi_.size(1) || ret.size(2) < dPhi_.size(2) ) return ;

  if ( bTranspose ) {
    for ( GINT i=0; i<dPhi_.size(1); i++ ) 
      for ( GINT j=0; j<dPhi_.size(2); j++ ) 
        ret(i,j) = static_cast<TE>(dPhiT_(i,j)*weights_[j]);
  } 
  else {
    for ( GINT i=0; i<dPhi_.size(1); i++ ) 
      for ( GINT j=0; j<dPhi_.size(2); j++ ) 
        ret(i,j) = static_cast<TE>(dPhi_(i,j)*weights_[i]);
  }
  
} // end of method getDerivMatrixW (2)


//************************************************************************************
//************************************************************************************
// METHOD : getDerivMatrixiW (1)
// DESC   : Get W^-1 X derivative matrix member data
// ARGS   : bTranspose: TRUE==>return transpose; else don't
// RETURNS: pointer to member data GTMatrix
//************************************************************************************
template<typename T, typename TE>
GTMatrix<TE> *GLLBasis<T,TE>::getDerivMatrixiW(GBOOL bTranspose)
{
  if ( bTranspose ) return &dPhiiWTEv_;
  else              return &dPhiiWEv_;
} // end of method getDerivMatrixiW (1)


//************************************************************************************
//************************************************************************************
// METHOD : getDerivMatrixiW (2)
// DESC   : Get deep copy of W^-1  X deriv matrix in specified return object
// ARGS   : ret: GTMatrix to return data
//          bTranspose: flag to get transpose (TRUE); else, not
// RETURNS: pointer to ret GTMatrix on success; else NULLPTR
//************************************************************************************
template<typename T, typename TE>
void GLLBasis<T,TE>::getDerivMatrixiW(GTMatrix<TE> &ret, GBOOL bTranspose)
{
  if ( ret.size(1) < dPhi_.size(1) || ret.size(2) < dPhi_.size(2) ) return ;

  if ( bTranspose ) {
    for ( GINT i=0; i<dPhi_.size(1); i++ ) 
      for ( GINT j=0; j<dPhi_.size(2); j++ ) 
        ret(i,j) = static_cast<TE>(dPhiT_(i,j)/weights_[j]);
  } 
  else {
    for ( GINT i=0; i<dPhi_.size(1); i++ ) 
      for ( GINT j=0; j<dPhi_.size(2); j++ ) 
        ret(i,j) = static_cast<TE>(dPhi_(i,j)/weights_[i]);
  }
  
} // end of method getDerivMatrixiW (2)


//************************************************************************************
//************************************************************************************
// METHOD : evalBasis (1)
// DESC   : 
//                h(xi)_i = 1/(N(N+1)*L_N(xi_i))  * (1-xi^2)*dL_N(xi)/dxi / (xi-xi_i)
// ARGS   : i  : which polynomial to evaluat (0,... Np_)
//          eta: reference interval value at which to evaluate 
// RETURNS: scalar result of evaluation   
//************************************************************************************
template<typename T, typename TE>
TE GLLBasis<T,TE>::evalBasis (GINT  i, TE eta)
{
  GString serr = "GLLBasis::evalBasis(1): ";
  T  ppn_i, pm1, pdm1, pm2, pdm2, ppn_xi, pder_xi, del;
  T  fact=-1.0/(Np_*(Np_+1.0)), gfact, xi;
  TE fRet;

  if ( !bInit_  && !init() ) {
    std::cout << serr << "GLLBasis::basis data incomplete" << std::endl;
    exit(1);
  }
  xi = static_cast<T>(eta);

  fRet = 1.0;
  del  = xi-xiNodes_[i];
  if ( xi < ximin_ || xi > ximax_ ) fRet = 0.0;
  else if ( fabs(del) > tetiny_ ) {
    ppn_i = Pn_[i];
    computeJacobi(Np_, alpha_, beta_, ppn_xi, pder_xi,pm1, pdm1, pm2, pdm2, xi);
    gfact = (1.0-xi*xi)/del; 
    fRet  = static_cast<TE>( fact * gfact * pder_xi/ppn_i );
  }
  return fRet;

} // end of method evalBasis (1)


//************************************************************************************
//************************************************************************************
// METHOD     : evalBasis (2)
// DESC   : 
//                h(xi)_i = 1/(N(N+1)*L_N(xi_j))  * (1-xi^2)*dL_N(xi)/dxi / (xi-xi_j)
// ARGS   : i   : which polynomial to evaluat (0,... Np_)
//          eta : deref GTVector to hold input ref interval points
//          vret: GTVector to hold results of evaluations at eta
// RETURNS: pointer to vret on success; else NULLPTR  
//************************************************************************************
template<typename T, typename TE>
GTVector<TE> *GLLBasis<T,TE>::evalBasis (GINT  i, GTVector<TE> &eta, GTVector<TE> &vret)
{
  GString serr = "GLLBasis::evalBasis(2): ";
  GINT  j;
  T     ppn_i, pm1, pdm1, pm2, pdm2, ppn_xi, pder_xi;
  T     fact=-1.0/(Np_*(Np_+1.0)), gfact, xi;
  TE    fRet;

  if ( !bInit_ &&  !init() ) {
    std::cout << serr << "basis data incomplete" << std::endl;
    exit(1);
  }

  for ( j=0; j<eta.size(); j++ ) {
    xi = static_cast<T>(eta[j]);
    fRet = 1.0;
    if ( xi < ximin_ || xi > ximax_ ) fRet = 0.0;
    else if ( fabs(xi-xiNodes_[i]) > tetiny_ ) {
      ppn_i = Pn_[i];
      computeJacobi(Np_, alpha_, beta_, ppn_xi, pder_xi,pm1, pdm1, pm2, pdm2, xi);
      gfact = (1.0-xi*xi)/(xi-xiNodes_[i]+ttiny_); 
      fRet  = static_cast<TE>( fact * gfact * pder_xi/(ppn_i+ttiny_) );
    }
    vret[j] = fRet;
  }
  return &vret;

} // end of method evalBasis (2)


//************************************************************************************
//************************************************************************************
// METHOD : evalBasis (3)
// DESC   : Evaluates basis at input parent domain points , eta_i
//          For Gauss-Lobatto, the basis is:
//          h_j(eta) = -1/(Np_*(Np_+1)) * (1-eta**2) dL_Np_ (eta)dxi / (L_Np_(xi_j) (eta-xi_j))
// ARGS   : eta  : ref interval points at which to eval all polynomials
//          mret : GTMatrix return for evaluated polynomials
// RETURNS: pointer to GTMatrix, M_ij = h_j(eta_i) on success; else NULLPTR
//************************************************************************************
template<typename T, typename TE>
GTMatrix<TE> *GLLBasis<T,TE>::evalBasis (GTVector<TE> &eta, GTMatrix<TE> &mret)
{
  GString serr = "GLLBasis::evalBasis(3): ";
  GSIZET i, j;
  T      ppn_j, pm1, pdm1, pm2, pdm2, ppn_xi, pder_xi;
  T      fact=-1.0/(Np_*(Np_+1.0)), gfact, xi;
  TE     fRet;

  if ( !bInit_ && !init() ) {
    std::cout << serr << "evalBasis (3): basis data incomplete" << std::endl;
    exit(1);
  }

  for ( i=0; i<eta.size(); i++ ) {
    for ( j=0; j<Np_+1; j++ )  {
      fRet = 1.0;
      xi    = static_cast<T>(eta[i]);
      if ( xi < ximin_ || xi > ximax_ ) fRet = 0.0;
      else if ( fabs(xi-xiNodes_[j]) > tetiny_ ) {
        ppn_j = Pn_[j];
        computeJacobi(Np_, alpha_, beta_, ppn_xi, pder_xi,pm1, pdm1, pm2, pdm2, xi);
        gfact = (1.0-xi*xi)/(xi-xiNodes_[j]); 
        fRet  = static_cast<TE>( fact * gfact * pder_xi/ppn_j );
      }
      mret(i,j) = fRet;
    }
  }
  return &mret;

} // end of method evalBasis (3)

//************************************************************************************
//************************************************************************************
// METHOD : evalBasis (4)
// DESC   : Evaluates basis at input parent domain points , eta_i
//          For Gauss-Lobatto, the basis is:
//          h_j(eta) = -1/(Np_*(Np_+1)) * (1-eta**2) dL_Np_ (eta)dxi / (L_Np_(xi_j) (eta-xi_j))
// ARGS   : eta  : ref interval points at which to eval all polynomials
//          vret : GTVector return for evaluated polynomials; same  
//                 storage format as in evalBasis(3)
// RETURNS: pointer to GTMatrix, M_ij = h_j(eta_i) on success; else NULLPTR
//************************************************************************************
template<typename T, typename TE>
GTVector<TE> *GLLBasis<T,TE>::evalBasis (GTVector<TE> &eta, GTVector<TE> &vret)
{
  GString serr = "GLLBasis::evalBasis(4): ";
  GINT  i, j;
  T     ppn_j, pm1, pdm1, pm2, pdm2, ppn_xi, pder_xi;
  T     fact=-1.0/(Np_*(Np_+1.0)), gfact, xi;
  TE    fRet;

  if ( !bInit_ && !init() ) {
    std::cout << serr << "evalBasis (4): basis data incomplete" << std::endl;
    exit(1);
  }

  for ( i=0; i<eta.size(); i++ ) {
    for ( j=0; j<Np_+1; j++ )  {
      fRet = 1.0;
      xi    = static_cast<T>(eta[i]);
      if ( xi < ximin_ || xi > ximax_ ) fRet = 0.0;
      else if ( fabs(xi-xiNodes_[j]) > tetiny_ ) {
        ppn_j = Pn_[j];
        computeJacobi(Np_, alpha_, beta_, ppn_xi, pder_xi,pm1, pdm1, pm2, pdm2, xi);
        gfact = (1.0-xi*xi)/(xi-xiNodes_[j]); 
        fRet  = static_cast<TE>( fact * gfact * pder_xi/ppn_j );
      }
      vret[i+j*(Np_+1)] = fRet;
    }
  }
  return &vret;

} // end of method evalBasis (4)


//************************************************************************************
//************************************************************************************
// METHOD : evalBasis (5)
// DESC   : Evaluates basis at input parent domain points , eta_i
//          For Gauss-Lobatto, the basis is:
//          h_j(eta) = -1/(Np_*(Np_+1)) * (1-eta**2) dL_Np_ (eta)dxi / (L_Np_(xi_j) (eta-xi_j))
// ARGS   : eta  : array of ref interval points, xi_j
//          neta : num elements in eta array
//          mret : GTMatrix return for evaluated polynomials
// RETURNS: pointer to return matrix, mret on success; else NULPTR 
//************************************************************************************
template<typename T, typename TE>
GTMatrix<TE> *GLLBasis<T,TE>::evalBasis (TE eta[], GINT neta, GTMatrix<TE> &mret)
{
  GString serr = "GLLBasis<T,TE>::evalBasis: ";
  GINT  i, j;
  T     ppn_j, pm1, pdm1, pm2, pdm2, ppn_xi, pder_xi;
  T     fact=-1.0/(Np_*(Np_+1.0)), gfact, xi;
  TE    fRet;

  if ( !bInit_ && !init() ) {
    std::cout << serr << "evalBasis (5): basis data incomplete" << std::endl;
    exit(1);
  }

  for ( i=0; i<neta; i++ ) {
    for ( j=0; j<Np_+1; j++ )  {
      fRet = 1.0;
      xi    = static_cast<T>(eta[i]);
      if ( xi < ximin_ || xi > ximax_ ) fRet = 0.0;
      else if ( fabs(xi-xiNodes_[j]) > tetiny_ ) {
        ppn_j = Pn_[j];
        computeJacobi(Np_, alpha_, beta_, ppn_xi, pder_xi,pm1, pdm1, pm2, pdm2, xi);
        gfact = (1.0-xi*xi)/(xi-xiNodes_[j]); 
        fRet  = static_cast<TE>( fact * gfact * pder_xi/ppn_j );
      }
      mret(i,j) = fRet;
    }
  }

  return &mret;

} // end of method evalBasis (5)


//************************************************************************************
//************************************************************************************
// METHOD : evalDBasis (1)
// DESC   : Evaluates basis j, derivative at input parent domain point , eta
//              Deriv. is derived from :
//              dh_j(eta)/dxi =  -1/(Np_*(Np_-1)) * (1-eta**2) dL_Np_ (eta)dxi / (L_Np_(xi_j) (eta-xi_j))
// ARGS   : 
//           j   : basis function to evaluate derivative of
//           eta : refereince interval at which evaluate
// RETURNS: derivative value 
//************************************************************************************
template<typename T, typename TE>
TE GLLBasis<T,TE>::evalDBasis (GINT j, TE eta)
{
  GString serr = "GLLBasis::evalDBasis(1): ";
  T     ppn_j, pm1, pdm1, pm2, pdm2, ppn_xi, pder_xi, pdd;
  T     fact=-1.0/(Np_*(Np_+1.0)), gfact, g1, g1i, g2, xi ;
  TE    fRet;
  
  if ( !bInit_ && !init() ) {
    std::cout << serr << "basis data incomplete" << std::endl;
    exit(1);
  }

  fRet = 0.0;
  xi     = eta;
  g1     = xi - xiNodes_[j];
  if      ( xi == ximin_ && j == 0 ) {
    fRet = -0.25*Np_*(Np_+1.0); 
  }
  else if ( xi == ximax_ && j == Np_ ) {
    fRet = 0.25*Np_*(Np_+1.0); 
  }
  else if ( fabs(g1) > tetiny_) {
    ppn_j = Pn_[j];
    computeJacobi(Np_, alpha_, beta_, ppn_xi, pder_xi,pm1, pdm1, pm2, pdm2, xi);
    gfact = fact / ppn_j;
    g1i   = 1.0/g1;
    g2    = (1.0 - xi*xi)*g1i;
    pdd   = 2.0*xi*pder_xi - Np_*(Np_ + 1.)*ppn_xi; 
    fRet  = static_cast<TE>(gfact * g1i * ( pdd - (2.0*xi + g2 )*pder_xi) );
  }
  return fRet;

} // end of method evalDBasis (1)


//************************************************************************************
//************************************************************************************
// METHOD : evalDBasis (2)
// DESC   : Evaluates basis derivative at input parent domain points , eta_i
//              Deriv. is derived from :
//              dh_j(eta)/dxi =  -1/(Np_*(Np_-1)) * (1-eta**2) dL_Np_ (eta)dxi / (L_Np_(xi_j) (eta-xi_j))
// ARGS   : eta  : GTVector of ref interval points at which to evaluate derivative basis
//          mret : GTMatrix return object to hold evaluation results
// RETURNS: pointer to mret on success; else NULLPTR
//************************************************************************************
template<typename T, typename TE>
GTMatrix<TE> *GLLBasis<T,TE>::evalDBasis (GTVector<TE> &eta, GTMatrix<TE> &mret)
{
  GString serr = "GLLBasis::evalDBasis(2): ";
  GINT  i, j, mm, nn;
  T     ppn_j, pm1, pdm1, pm2, pdm2, ppn_xi, pder_xi, pdd;
  T     fact=-1.0/(Np_*(Np_+1.0)), gfact, g1, g1i, g2, xi ;
  TE    fRet;
  
  if ( !bInit_ && !init() ) {
    std::cout << serr << "basis data incomplete" << std::endl;
    exit(1);
  }

  nn = MIN(eta.size(),mret.size(1));
  mm = MIN(Np_+1,mret.size(2)); 
  for ( i=0; i<nn; i++) {
    xi     = static_cast<T>(eta[i]);
    for ( j=0; j<mm;  j++) {
      fRet = 0.0;
      g1     = xi - xiNodes_[j];
      if      ( xi == ximin_ && j == 0 ) {
        fRet = -0.25*Np_*(Np_+1.0); 
      }
      else if ( xi == ximax_ && j == Np_ ) {
        fRet = 0.25*Np_*(Np_+1.0); 
      }
      else if ( fabs(g1) > tetiny_) {
        ppn_j = Pn_[j];
        computeJacobi(Np_, alpha_, beta_, ppn_xi, pder_xi,pm1, pdm1, pm2, pdm2, xi);
        gfact = fact / ppn_j;
        g1i   = 1.0/g1;
        g2    = (1.0 - xi*xi)*g1i;
        pdd   = 2.0*xi*pder_xi - Np_*(Np_ + 1.)*ppn_xi; 
        fRet  = static_cast<TE>(gfact * g1i * ( pdd - (2.0*xi + g2 )*pder_xi) );
      } 
      mret(i,j) = fRet;
    }
  }
  return &mret;

} // end of method evalDBasis (2)


//************************************************************************************
//************************************************************************************
// METHOD : evalDBasis (3)
// DESC   : Evaluates basis derivative at input parent domain points , eta_i
//          Deriv. is derived from :
//          h_j(eta) =  -1/(Np_*(Np_-1)) * (1-eta**2) dL_Np_ (eta)dxi / (L_Np_(xi_j) (eta-xi_j))
// ARGS   : eta : array of reference interval points at which to evaluate
//          n   : size of eta array
//          mret: Matrix of basis evaluations, M_ij = dh_j(eta_i)/dxi
// RETURNS: Pointer to mret on success; else NULLPTR
//************************************************************************************
template<typename T, typename TE>
GTMatrix<TE> *GLLBasis<T,TE>::evalDBasis (TE eta[], GINT n, GTMatrix<TE> &mret)
{
  GString serr = "GLLBasis::evalDBasis(3): ";
  GINT  i, j, mm, nn;
  T     ppn_j, pm1, pdm1, pm2, pdm2, ppn_xi, pder_xi, pdd;
  T     fact=-1.0/(Np_*(Np_+1.0)), gfact, g1, g1i, g2, xi ;
  TE    fRet;
  
  if ( !bInit_ && !init() ) {
    std::cout << serr << "basis data incomplete" << std::endl;
    exit(1);
  }

  nn = MIN(n,mret.size(1));
  mm = MIN(Np_+1,mret.size(2)); 
  for ( i=0; i<nn; i++) { // loop over nodes
    xi     = static_cast<T>(eta[i]);
    for ( j=0; j<mm;  j++) { // loop over modes
      fRet = 0.0;
      g1     = xi - xiNodes_[j];
      if      ( xi == ximin_ && j == 0 ) {
        fRet = -0.25*Np_*(Np_+1.0); 
      }
      else if ( xi == ximax_ && j == Np_ ) {
        fRet = 0.25*Np_*(Np_+1.0); 
      }
      else if ( fabs(g1) > tetiny_ ) {
        ppn_j = Pn_[j];
        computeJacobi(Np_, alpha_, beta_, ppn_xi, pder_xi,pm1, pdm1, pm2, pdm2, xi);
        gfact = fact / ppn_j;
        g1i   = 1.0/g1;
        g2    = (1.0 - xi*xi)*g1i;
        pdd   = 2.0*xi*pder_xi - Np_*(Np_ + 1.)*ppn_xi; 
        fRet  = static_cast<TE>( gfact * g1i * ( pdd - (2.0*xi + g2 )*pder_xi) );
      } 
      mret(i,j) = fRet;
    }
  }
  return &mret;
} // end of method evalDBasis (3)


//************************************************************************************
//************************************************************************************
// METHOD : evalDBasis (4)
// DESC   : Evaluates j-th basis derivative at input parent domain points , eta_i
//          Deriv. is derived from :
//          h_i(eta) =  -1/(Np_*(Np_-1)) * (1-eta**2) dL_Np_ (eta)dxi / (L_Np_(xi_j) (eta-xi_j))
// ARGS   : j   : which polynomial to evaluat (0,... Np_)
//          eta : deref GTVector to hold input ref interval points
//          vret: GTVector to hold results of evaluations at eta
// RETURNS: Pointer to vret on success; else NULLPTR
//************************************************************************************
template<typename T, typename TE>
GTVector<TE> *GLLBasis<T,TE>::evalDBasis (GINT j, GTVector<TE> &eta, GTVector<TE> &vret)
{
  GString serr = "GLLBasis::evalDBasis(4): ";
  GINT  i;
  T     ppn_j, pm1, pdm1, pm2, pdm2, ppn_xi, pder_xi, pdd;
  T     fact=-1.0/(Np_*(Np_+1.0)), gfact, g1, g1i, g2, xi ;
  TE    fRet;
  
  if ( !bInit_ && !init() ) {
    std::cout << serr << "basis data incomplete" << std::endl;
    exit(1);
  }


  for ( i=0; i<eta.size(); i++) {
    fRet = 0.0;
    xi     = static_cast<T>(eta[i]);
    g1     = xi - xiNodes_[j];
    if      ( xi == ximin_ && j == 0 ) {
      fRet = -0.25*Np_*(Np_+1.0); 
    }
    else if ( xi == ximax_ && j == Np_ ) {
      fRet = 0.25*Np_*(Np_+1.0); 
    }
    else if ( fabs(g1) > tetiny_) {
      ppn_j = Pn_[j];
      computeJacobi(Np_, alpha_, beta_, ppn_xi, pder_xi,pm1, pdm1, pm2, pdm2, xi);
      gfact = fact / ppn_j;
      g1i   = 1.0/g1;
      g2    = (1.0 - xi*xi)*g1i;
      pdd   = 2.0*xi*pder_xi - Np_*(Np_ + 1.)*ppn_xi; 
      fRet  = static_cast<TE>(gfact * g1i * ( pdd - (2.0*xi + g2 )*pder_xi) );
    } 
    vret[i] = fRet;
  }

  return &vret;
} // end of method evalDBasis (4)


