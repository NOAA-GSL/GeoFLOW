//************************************************************************************//
// Module       : gregelem.hpp
// Date         : 9/14/01 (DLR)
// Copyright    : 2001-2006 Copyright University Corporation for Atmospheric
//                Research
// Description  : Encapsulates the methods and data associated with
//                a single regular rectangular 2D spectral element
//                The following is the ordering for the vertices (Vi)
//                and segment midpoints (Mi):
//
//          V7 ________________ V6
//            /|     M6       /|
//       M7  / |             / |
//          /  |        M5  /  |
//         /M11|  M4       /   |
//     V4 /____|__________/    | M10
//        |    |         | V5  |
//        |    |     M2  |     |
//        | V3 |_________|_____| V2
//    M8  |    /         |    /
//        |   /       M9 |   /
//        |  / M3        |  / M1
//        | /            | /
//        |______________|/
//       V0       M0     V1
//
// Faces are labeled s.t. F0-F3 correspond to orientation of edges on bottom plane;
// F4 and F5 are bottom and top faces, respectively.
// [Note that in 2d, we use just the bottom plane.]

// Derived From : GElemB.
// Modifications:
//************************************************************************************//
#include "gregelem.hpp"
#include <cstdlib>
#include <memory>
#include <cmath>
#include <cstdio>
#include "gcutils.hpp"

//************************************************************************************
//************************************************************************************
// Constructor Method (1)
//************************************************************************************
GRegElem::GRegElem()
{
  elemtype_ = RECT_QUAD;

} // end of constructor method (1)

//************************************************************************************
//************************************************************************************
// Constructor Method (2)
//************************************************************************************
GRegElem::GRegElem(GNBasis *inBasis1, GNBasis *inBasis2, GNBasis *inBasis3)
: GElemB(inBasis1, inBasis2, inBasis3)
{
  elemtype_ = RECT_QUAD;
  CreateBaseElemDynamic();
} // end of constructor method (2)

//************************************************************************************
//************************************************************************************
// Constructor Method (3)
//************************************************************************************
GRegElem::GRegElem(GNBasis *b[], GINT nb)
: GElemB(b,nb)
{
  elemtype_ = RECT_QUAD;
  CreateBaseElemDynamic();
} // end of constructor method (2)


//************************************************************************************
//************************************************************************************
// Destructor
//************************************************************************************
GRegElem::~GRegElem()
{
   DeleteDynamic();
}

//************************************************************************************
//************************************************************************************
// Assignment operator method (1)
//************************************************************************************
void GRegElem::operator=(const GRegElem &elem)
{
  GINT  i;

  if ( &elem != this ) {
    // copy data:
    for ( i=0; i<GDIM; i++ ) {
      Np_       [i] = elem.Np_[i];
      gllbasis_ [i] = elem.gllbasis_[i];
      glbasis_  [i] = elem.glbasis_ [i];
      xiNodes_  [i] = gllbasis_[i]->GetXiNodes(); 
      Weights_  [i] = gllbasis_[i]->GetWeights(); 
      D1d_      [i] = gllbasis_[i]->GetDerivMatrix(); 
      D1dT_     [i] = gllbasis_[i]->GetDerivMatrix(TRUE); 
      StiffMat_ [i] = gllbasis_[i]->GetDerivMatrix(); 
      StiffMatT_[i] = gllbasis_[i]->GetDerivMatrix(TRUE); 
    }
    bSolved_  = elem.bSolved_;

  }

} // end of = operator


//************************************************************************************
//************************************************************************************
// METHOD     : ComputeSpNodes
// DESCRIPTION: Computes real space locations of parent space nodes.
//              Two matrices are created, spNodes1(2), of the form
//                x_i(0,0), x_i(1,0), ..., x_i(Np_[0],0),
//                x_i(0,1), ...
//                .
//                .
//                .
//                x_i(0,Np_[1]), x_i(1,Np_[0]),..., x_i(Np_[0],Np_[1])
//
//               where x_i(l,m) is x1(2) evaluated at parent node
//               xi1_l, xi2_m, and xi are increasing with index in
//               each direction. 
//            
// ARGUMENTS  :
// RETURNS    :  TRUE on success; else FALSE. Failure occurs if
//               nodes and weights cannot be computed, or if memory
//               cannot be allocated.
//************************************************************************************
GBOOL GRegElem::ComputeSpNodes()
{
  GDOUBLE xip[GDIM], xim[GDIM];
  GINT    i, j, k, r;
  char   *serr = "GRegElem::ComputeSpNodes: ";

  // NOTE: Geometry quantities assumed to be computed in SolveFE

  if ( !bSolved_ )
    if ( SolveFE() <= 0 ) return FALSE; 

#if defined(G_IS2D)
  for ( j=0; j<Np_[1]; j++ ) {
    xip[1] = 1.0 + (*xiNodes_[1])[j];
    xim[1] = 1.0 - (*xiNodes_[1])[j];
    for ( i=0; i<Np_[0]; i++ ) {
      xip[0] = 1.0 + (*xiNodes_[0])(i);
      xim[0] = 1.0 - (*xiNodes_[0])(i);
      r = i + j*(Np_[0]);
      (*spNodes_[0])(r) = 0.25* ( spVertices_[0].x1*xim[0]*xim[1] + spVertices_[1].x1*xip[0]*xim[1]
                                + spVertices_[3].x1*xim[0]*xip[1] + spVertices_[2].x1*xip[0]*xip[1] );
      (*spNodes_[1])(r) = 0.25* ( spVertices_[0].x2*xim[0]*xim[1] + spVertices_[1].x2*xip[0]*xim[1]
                                + spVertices_[3].x2*xim[0]*xip[1] + spVertices_[2].x2*xip[0]*xip[1] );
    } 
  }
#elif defined(G_IS3D)
  for ( k=0; k<Np_[2]; k++ ) {
    xip[2] = 1.0 + (*xiNodes_[2])[k];
    xim[2] = 1.0 - (*xiNodes_[2])[k];
    for ( j=0; j<Np_[1]; j++ ) {
      xip[1] = 1.0 + (*xiNodes_[1])[j];
      xim[1] = 1.0 - (*xiNodes_[1])[j];
      for ( i=0; i<Np_[0]; i++ ) {
        xip[0] = 1.0 + (*xiNodes_[0])(i);
        xim[0] = 1.0 - (*xiNodes_[0])(i);
        r = i + j*Np_[0]  + k*Np_[0]*Np_[1];
        (*spNodes_[0])(r) = 0.125* ( ( spVertices_[0].x1*xim[0]*xim[1] + spVertices_[1].x1*xip[0]*xim[1]
                                     + spVertices_[3].x1*xim[0]*xip[1] + spVertices_[2].x1*xip[0]*xip[1] )*xim[2] 
                                   + ( spVertices_[4].x1*xim[0]*xim[1] + spVertices_[5].x1*xip[0]*xim[1]
                                     + spVertices_[7].x1*xim[0]*xip[1] + spVertices_[6].x1*xip[0]*xip[1] )*xip[2] );
        (*spNodes_[1])(r) = 0.125* ( ( spVertices_[0].x2*xim[0]*xim[1] + spVertices_[1].x2*xip[0]*xim[1]
                                     + spVertices_[3].x2*xim[0]*xip[1] + spVertices_[2].x2*xip[0]*xip[1] )*xim[2] 
                                   + ( spVertices_[4].x2*xim[0]*xim[1] + spVertices_[5].x2*xip[0]*xim[1]
                                     + spVertices_[7].x2*xim[0]*xip[1] + spVertices_[6].x2*xip[0]*xip[1] )*xip[2] );
        (*spNodes_[2])(r) = 0.125* ( ( spVertices_[0].x3*xim[0]*xim[1] + spVertices_[1].x3*xip[0]*xim[1]
                                     + spVertices_[3].x3*xim[0]*xip[1] + spVertices_[2].x3*xip[0]*xip[1] )*xim[2] 
                                   + ( spVertices_[4].x3*xim[0]*xim[1] + spVertices_[5].x3*xip[0]*xim[1]
                                     + spVertices_[7].x3*xim[0]*xip[1] + spVertices_[6].x3*xip[0]*xip[1] )*xip[2] );
      } 
    }
  }
#endif
  return TRUE;
} // end of method ComputeSpNodes


#if 0
//************************************************************************************
//************************************************************************************
// METHOD     : GetBasisAtXi
// DESCRIPTION: Compute expansion basis at the parent domanin (xi) points.
//              This quantity really only has meaning for modal basis,
//              since for nodal bases, the basis is 1 at nodal points
//              (and 0 when not at nodal points). Provision is made for
//              this quantity, but only for specified 1- and 2-indices.
//              
// ARGUMENTS  : GINT  i, j, representing the 1- and 2-indices for labeling
//              the basis function, phi(i,j). The matrix ret is of type GMatrix,
//              and contains the basis fcn evaulated at the quadrature points, and
//              arranged s.t.
//
//              B(0,0) B(0,1)    ...   B(0,Np_[1])
//              B(1,0)           ...
//              .
//              .
//              B(Np_[0],0)         ...   B(Np_[0],Np_[1])
//
//              where the value of the quadrature 1- &  2- coord. increases with index.
// RETURNS    : GMatrix pointer of return matrix on success; else NULL.
//************************************************************************************
GMatrix *GRegElem::GetBasisAtXi(GINT  i, GINT  j, GMatrix *B)
{
  if ( B == NULL ) return NULL;
  if ( B->dim(1) != Np_[0]  || B->dim(2) != Np_[1] ) return NULL;

  if ( !bSolved_ )
    if ( !SolveFE() ) return NULL;

  GINT  l, m;
  GDOUBLE *b_data, *b1_data, *b2_data;
  GMatrix B1(Np_[0],Np_[0]), B2(Np_[1],Np_[1]);

  // Get 1- and 2- bases evaluated at the respective 1- and 2- quadrature points:
  if ( basis1->GetBasisAtNodes(&B1) == NULL ) return NULL;
  if ( basis2->GetBasisAtNodes(&B2) == NULL ) return NULL;
  
  // Note that each Ba from GetBasis is s.t.
  //   Ba(i,j) = Phi_i(xi_j)

  // When doing the multiplication here, it would be best not to make
  // the method calls with each (*,*) invocation, but instead work
  // direcly on the data.... This is done here in order to be explicit.
  b_data  = B->Data();
  b1_data = B1.Data();
  b2_data = B2.Data();
  for ( l=0; l<Np_[0]; l++ ) {
    for ( m=0; m<Np_[1]; m++ ) {
//    (*B)(l,m) = B1(i,l)*B2(j,m);
      *(b_data+i*(B->dim(2))+j) = (*(b1_data+i*Np_[0]+l)) * (*(b2_data+j*Np_[1]+m));
    }
  }
 
  return B;

} // end of method GetBasisAtXi
#endif


//************************************************************************************
//************************************************************************************
// METHOD     : SetVertices
// DESCRIPTION: Set elements real-space vertices. These should be
//              set in counter clockwise order, in 2d, as as in the above
//              figure in 3D.
// ARGUMENTS  : Array of GFPoint points representing vertex points.
//              There must be at least 4 of these, and only the
//              first 2*GDIM are used. Not very robust. 
// RETURNS    :  none
//************************************************************************************
GBOOL GRegElem::SetVertices(GFPoint p[], GINT  num)
{
  GINT  j, k, m, n;
  char *serr = "GRegElem::SetVertices: ";


  if ( num < nVertices_ ) {
    cout << serr << "Incorrect number of vertices" << endl;
    exit(1);
  }

  for ( j=0; j<nVertices_; j++ ) {
    spVertices_[j] = p[j]; 
  }
#if 0
  for ( n=0; n<GDIM*nVertices_; n+=GDIM ) {
    for ( j=0; j<GDIM; j++ ) spvVertices_[n] = spVertices[n/GDIM][j];
  }
#endif

  // Compute the edge segment midpoints:
  for ( n=0; n<4; n++ ) {   // bottom plane
    m = (n+1) % 4;
    for ( j=0; j<GDIM; j++) {
      spEMidpoints_[n][j]= 0.5*(spVertices_[m][j]+spVertices_[n][j]);
    }
  }

  // Compute element center:
  for ( j=0; j<GDIM; j++ )
    elemCenter_[j] = 0.25*( spVertices_[0][j]+spVertices_[1][j]
                   +        spVertices_[2][j]+spVertices_[3][j] );

#if defined(G_IS3D)
  for ( n=4; n<8; n++ ) {   // top plane, if in 3D
    m = 4 + (n+1) % 4;
    for ( j=0; j<GDIM; j++ ) spEMidpoints_[n][j]= 0.5*(spVertices_[m][j]+spVertices_[n][j]);
  }
  for ( n=8,k=0; n<12; n++,k++ ) { // vertical edges, if in 3D
    m = 0.5*(n+k);   
    for ( j=0; j<GDIM; j++ ) spEMidpoints_ [n][j]= 0.5*(spVertices_[m][j]+spVertices_[m-4][j]);
  }
  
  // Face midpoints: must be interesection of lines connecting 
  // bdy midpoints. For regular element, it can be written:
  for ( j=0; j<GDIM; j++ ) spFMidpoints_[0][j]= 0.25*(spVertices_[0][j]+spVertices_[5][j]
                                                    + spVertices_[1][j]+spVertices_[4][j]);
  for ( j=0; j<GDIM; j++ ) spFMidpoints_[1][j]= 0.25*(spVertices_[1][j]+spVertices_[6][j]
                                                    + spVertices_[2][j]+spVertices_[5][j]);
  for ( j=0; j<GDIM; j++ ) spFMidpoints_[2][j]= 0.25*(spVertices_[3][j]+spVertices_[6][j]
                                                    + spVertices_[2][j]+spVertices_[7][j]);
  for ( j=0; j<GDIM; j++ ) spFMidpoints_[3][j]= 0.25*(spVertices_[3][j]+spVertices_[4][j]
                                                    + spVertices_[0][j]+spVertices_[7][j]);
  for ( j=0; j<GDIM; j++ ) spFMidpoints_[4][j]= 0.25*(spVertices_[0][j]+spVertices_[2][j]
                                                    + spVertices_[1][j]+spVertices_[3][j]);
  for ( j=0; j<GDIM; j++ ) spFMidpoints_[5][j]= 0.25*(spVertices_[4][j]+spVertices_[6][j]
                                                    + spVertices_[5][j]+spVertices_[7][j]);

  for ( j=0; j<GDIM; j++ ) {
    elemCenter_[j] = 0.125*(
          spVertices_[0][j]+spVertices_[1][j]+spVertices_[2][j]+spVertices_[3][j] 
        + spVertices_[4][j]+spVertices_[5][j]+spVertices_[6][j]+spVertices_[7][j] 
                           );
  }
#endif
  

  return TRUE;

} // end of method SetVertices 


//************************************************************************************
//************************************************************************************
// METHOD     : DeleteDynamic
// DESCRIPTION: deletes dynamically allocated quantities
// ARGUMENTS  :
// RETURNS    :  none
//************************************************************************************
void GRegElem::DeleteDynamic()
{

  for ( GINT i=0; i<GDIM; i++ ) {
    if ( opInterp_ [i]   != NULL ) delete opInterp_ [i];
    if ( opInterpT_[i]   != NULL ) delete opInterpT_[i];
    if ( dInterp_  [i]   != NULL ) delete dInterp_  [i];
    if ( dInterpT_ [i]   != NULL ) delete dInterpT_ [i];
  }

} // end of method DeleteDynamic


//************************************************************************************
//************************************************************************************
// METHOD     : XToXi() 
// DESCRIPTION: converts the real space values into reference-domain values,
//              on [-1,1].
// ARGUMENTS  : Point *x = pointer to real space points;
//              Point *xi= pointer to parent domain points
//              GINT  num= number of points to invert
//
//              NOTE: there will be problems if the number of
//              points in each array is not the same.
// RETURNS    : TRUE on success; FALSE on failure if the 
//              solution is imaginary, or undefined, or if input
//              point lies outside of element.
//************************************************************************************
GBOOL GRegElem::XToXi(GFPoint x[], GFPoint xi[], const GINT num)
{
  GINT     i, j;
  GDOUBLE  eps, L[GDIM], iL[GDIM], x0[GDIM];
  char    *serr = "GRegElem::XToXi: ";

  // Check that point is within element boundaries:
  if ( !Point_in_elem(x, num) ) {
    cout << serr <<  "One or more points lies outside element!" << endl;
    return FALSE;
  }

  if ( !bSolved_ )
    if ( !SolveFE() ) return FALSE;
  

  L [0] = spVertices_[1].x1 - spVertices_[0].x1;
  L [1] = spVertices_[3].x2 - spVertices_[0].x2;
  x0[0] = spVertices_[0].x1;
  x0[1] = spVertices_[0].x2;
#if defined(G_IS3D)
  L [2] = spVertices_[4].x3 - spVertices_[0].x3;
  x0[2] = spVertices_[0].x3;
#endif
  for ( j=0; j<GDIM; j++ ) iL[j] = 1.0 / L[j];

  eps = 10.0 * GTINY;
  for ( i=0; i<num; i++ ) {
//  xi[i].Bracket(GPTINY);
    for ( j=0; j<GDIM; j++ ) {
      xi[i][j]   = 2.0*( x[i][j] - x0[j] ) * iL[j] - 1.0;
      if ( xi[i][j] <  1+eps  && xi[i][j] >  1.0-eps ) xi[i][j] = 1.0;
      if ( xi[i][j] < -1+eps  && xi[i][j] > -1.0-eps ) xi[i][j] = -1.0;
    }
//  xi[i].Truncate(); // remove all data outside xi' bracket
  }

  return TRUE;
} // end of method XToXi 



//************************************************************************************
//************************************************************************************
// METHOD     : SolveFE
// DESCRIPTION: computes quadrature points, weights, deriv. matrices 
//              for element. These are computed for the 1d component bases, 
//              and scatter arrays are set up that can determine the
//              V, E, and I parent nodal locations from the components.
// ARGUMENTS  :
// RETURNS    : integer status flag
//************************************************************************************
GINT GRegElem::SolveFE()
{
  GINT     i, ir, iret=1, I[3], j, k, r;
  GDOUBLE  Jac;
  char    *serr = "GRegElem::SolveFE: ";
 
  bSolved_ = TRUE;

 
  // Compute quadrature points and corresponding weights for each coordinate
  // basis:

  for ( i=0; i<GDIM; i++ ) {
    xiNodes_  [i] = gllbasis_[i]->GetXiNodes(); 
    Weights_  [i] = gllbasis_[i]->GetWeights(); 
    D1d_      [i] = gllbasis_[i]->GetDerivMatrix(); 
    D1dT_     [i] = gllbasis_[i]->GetDerivMatrix(TRUE); 
    StiffMat_ [i] = gllbasis_[i]->GetStiffMatrix(); 
    StiffMatT_[i] = gllbasis_[i]->GetStiffMatrix(TRUE); 
  }

  Jac     = 1.0 / pow(2.0,GDIM);
  volume_ = 1.0;
  for ( i=0; i<GDIM; i++ ) {
    L_[i]    = fabs(spVertices_[i+1][i] - spVertices_[i][i]);
    if ( i == 2 ) 
    L_[i]    = fabs(spVertices_[4][i] - spVertices_[0][i]);

    if ( L_[i] <= 0.0 ) {
      cout << serr << "Invalid element length" << endl;
      cout << serr << "Element=" << *this << endl;
      exit(1);
    }
    Jac     *= L_[i];
    volume_ *= L_[i];
  }

  for ( k=0; k<Np_[2]; k++ ) {
    I[2] = k;
    for ( j=0; j<Np_[1]; j++ ) {
      I[1] = j;
      for ( i=0; i<Np_[0]; i++ ) {
        r = i + j*Np_[0] + k*Np_[0]*Np_[1];
        I[0] = i;
        MassMatrix_[r] = Jac;
        for ( ir=0; ir<GDIM; ir++ ) MassMatrix_[r] *= (*Weights_[ir])[I[ir]];
//      MassMatrix_     (k) = (*Weights_[0])(i) * (*Weights_[1])(j)* Jac; 
      }
    }
  }
 
  
  gMassMatrix_ = MassMatrix_;  // make deep copy

  bSolved_ = ComputeSpNodes();
  InitMortars();

  // Compute dealiasing qantities:
  ComputeDealias();

  return iret;
  
} // end of method SolveFE


//************************************************************************************
//************************************************************************************
// METHOD     : SetInterpBasis
// DESCRIPTION: 
// ARGUMENTS  :
// RETURNS    : 
//************************************************************************************
void GRegElem::SetInterpBasis(GNBasis *b[], GINT nb)
{
  GINT      i, j, k, pNp[GDIM], r;
  GBOOL     bValid=TRUE;
  GVector  *pxi[GDIM], *gW[GDIM];
  char     *serr = "GRegElem::SetInterpBasis: ";


  if ( nb < GDIM ) {
    cout << serr << "Incorrect number of interpolation bases" << endl;
    exit(1);
  }
  for ( i=0; i<GDIM; i++ ) glbasis_[i] = b[i];

  for ( i=0; i<GDIM; i++ ) {
    if ( glbasis_[i] == NULL ) {
      cout << serr << "NULL p-basis" << endl;
      exit(1);
    }
    pNp[i] = glbasis_[i]->GetOrder() + 1;
  }

  if ( !bSolved_ )
    if ( !SolveFE() ) return;                    // ensure that all is ok with gll basis

  // If basis to which to interpolate required 
  // quantities is valid, allocate and compute the 
  // required quantities:

  
  // Instantiate or resize interp. basis quantities:
  for ( i=0; i<GDIM; i++ ) {
    if ( opInterp_ [i] == NULL  ) opInterp_ [i] = new GMatrix(pNp [i], Np_[i]); 
    else                          opInterp_ [i]->      Resize(pNp [i], Np_[i]);
    if ( opInterpT_[i] == NULL  ) opInterpT_[i] = new GMatrix( Np_[i],pNp [i]); 
    else                          opInterpT_[i]->      Resize( Np_[i],pNp [i]);
    if ( dInterp_  [i] == NULL  ) dInterp_  [i] = new GMatrix(pNp [i], Np_[i]); 
    else                          dInterp_  [i]->      Resize(pNp [i], Np_[i]);
    if ( dInterpT_ [i] == NULL  ) dInterpT_ [i] = new GMatrix( Np_[i],pNp [i]); 
    else                          dInterpT_ [i]->      Resize( Np_[i],pNp [i]);
    vptmp_[i].Resize( pNp[i]*Np_[i] );
    Ivp_  [i].Resize(  Np_[i],  pNp[i]);
    IvpT_ [i].Resize(  pNp[i],  Np_[i]);
    Ipv_  [i].Resize(  pNp[i],  Np_[i]);
    IpvT_ [i].Resize(  Np_[i],  pNp[i]);
    pxi   [i] = glbasis_[i]->GetXiNodes();
    gW    [i] = glbasis_[i]->GetWeights();
  }
  
  for ( i=0; i<GDIM; i++ ) {
    gllbasis_[i]->EvalBasis (pxi[i], opInterp_[i]);
    gllbasis_[i]->EvalDBasis(pxi[i], dInterp_ [i]);
    gllbasis_[i]->EvalBasis (pxi[i],&Ipv_     [i]);
    Ipv_[i].Transpose(IpvT_[i]); 
    glbasis_ [i]->EvalBasis (pxi[i],&Ivp_     [i]);
    Ivp_[i].Transpose(IvpT_[i]); 
  } 

  // Interp and derivative operators include the weights:
  for ( r=0; r<GDIM; r++ ) {
    for ( j=0; j<opInterp_[r]->dim(2); j++ ) {
      for ( i=0; i<opInterp_[r]->dim(1); i++ ) {
        (*opInterp_[r])(i,j) *= (*gW[r])(i);
        (*dInterp_ [r])(i,j) *= (*gW[r])(i);
      }
    }
    opInterp_[r]->Transpose(*opInterpT_[r]);
    dInterp_ [r]->Transpose(*dInterpT_ [r]);
  }
  
  
} // end of method SetInterpBasis


//************************************************************************************
//************************************************************************************
// METHOD     : Integrate
// DESCRIPTION: Computes integral of input vector over full geometric element.
//              s.t.
//                int  = Int   v  dx dy 
//              Note: this could be placed in the base class, if MassMatrix_ is
//                    also included.
// ARGUMENTS  : v            : integrand as a GVector
//              imultiplicity: inverse of node multiplicity; used if non-NULL. 
//                             Must be of the same length as v.
// RETURNS    :
//************************************************************************************
GDOUBLE GRegElem::Integrate(GVector *v, GDOUBLE *imultiplicity)
{

  GINT    j;
  GDOUBLE sum;
  char   *serr = "GRegElem::Integrate: ";

  if ( !bSolved_ &&  !SolveFE() ) {
    cout << serr << "Element solve failed" << endl;
    exit(1);
  }

  if ( v == NULL ) {
    cout << serr << "NULL integrand" << endl;
    exit(1);
  }

  if ( imultiplicity == NULL ) {
    sum = MTK::fvec_dot(*v, MassMatrix_);
  }
  else {
    for ( j=0, sum=0.0; j<v->dim(); j++ ) {
      sum += (*v)(j)*MassMatrix_(j)*imultiplicity[j];
    }
  }
  
  return sum;

} // end of method Integrate


//************************************************************************************
//************************************************************************************
// METHOD     : PIntegrate
// DESCRIPTION: Computes integral of input vector over 'parent' or reference
//              element subdomain
//              s.t.
//                int  = Int_{-1}^{1}   v  dxi_1 dxi_2 
//              Note: this could be placed in the base class, if Weights_ are
//                    also included.
// ARGUMENTS  : v            : integrand as a GVector
//              imultiplicity: node multiplicity inverse; used if non-NULL. Must be of the same
//                            length as v.
// RETURNS    :
//************************************************************************************
GDOUBLE GRegElem::PIntegrate(GVector *v, GDOUBLE *imultiplicity)
{

  GINT     i, j, k, r;
  GDOUBLE  sum;
  GVector *W[GDIM];
  char    *serr = "GRegElem::PIntegrate: ";

  if ( !bSolved_ &&  !SolveFE() ) {
    cout << serr << "Element solve failed" << endl;
    exit(1);
  }

  if ( v == NULL ) {
    cout << serr << "NULL integrand" << endl;
    exit(1);
  }

  for ( i=0; i<GDIM; i++ ) W[i] = gllbasis_[i]->GetWeights();

  if ( imultiplicity == NULL ) {
    for ( k=0, r=0, sum=0.0; k<Np_[2]; k++ ) 
      for ( j=0; j<Np_[2]; j++ ) 
        for ( i=0; i<Np_[2]; i++, r++ ) 
#if defined(G_IS2D)
          sum += (*v)[r] * (*W[0])[i] * (*W[1])[j];
#elif defined(G_IS3D)
          sum += (*v)[r] * (*W[0])[i] * (*W[1])[j] * (*W[2])[k];
#endif
  }
  else {
    for ( k=0, r=0, sum=0.0; k<Np_[2]; k++ ) 
      for ( j=0; j<Np_[2]; j++ ) 
        for ( i=0; i<Np_[2]; i++, r++ ) 
#if defined(G_IS2D)
          sum += (*v)[r] * (*W[0])[i] * (*W[1])[j];
#elif defined(G_IS3D)
          sum += (*v)[r] * (*W[0])[i] * (*W[1])[j] * (*W[2])[k] * imultiplicity[r];
#endif
  }
  
  return sum;

} // end of method PIntegrate


//************************************************************************************
//************************************************************************************
// METHOD     : Differentiate
// DESCRIPTION: Computes derivative of input vector in direction idir
// ARGUMENTS  : dv  : output arguement: derivative in idir direction of v
//              v   : input argument as a GVector
//              idir: coordinate direction
// RETURNS    : TRUE on success; else FALSE
//************************************************************************************
inline GBOOL GRegElem::Differentiate(GVector *dv, GVector *v, GINT  idir)
{
  GINT    i;
  GDOUBLE iL[GDIM];
  char   *serr = "GRegElem::Differentiate: ";

  if ( v == NULL || dv == NULL ) {
    cout << serr << "NULL vectors not allowed" << endl;
    exit(1);
  }

  for ( i=0; i<GDIM; i++ ) iL[i] = 2.0/L_[i];


  switch ( idir ) {
    case 1:
#if defined(G_IS2D)
      // Do I2 X D1 term: 
      MTK::I2_X_D1(*D1d_[0], *v, Np_[0], Np_[1], *dv);
#elif defined(G_IS3D)
      // Do I3 X I2 X D1 term: 
      MTK::I3_X_I2_X_D1(*D1d_[0], *v, Np_[0], Np_[1], Np_[2], *dv);
#endif
      break;
    case 2: 
#if defined(G_IS2D)
      // Do D2 X I1 term: 
      MTK::D2_X_I1(*D1dT_[1], *v, Np_[0], Np_[1], *dv);
#elif defined(G_IS3D)
      // Do I3 X D2 X I1 term: 
      MTK::I3_X_D2_X_I1(*D1dT_[1], *v, Np_[0], Np_[1], Np_[2], *dv);
#endif
      break;
#if defined(G_IS3D)
    case 3:
      // Do D3 X I2 X I1 term: 
      MTK::D3_X_I2_X_I1(*D1dT_[2], *v, Np_[0], Np_[1], Np_[2], *dv);
      break;
#endif
  }
  if ( iL[idir-1] != 1.0 ) MTK::fvec_const_prod_rep(*dv, iL[idir-1]);

  return TRUE;

} // end of method Differentiate


//************************************************************************************
//************************************************************************************
// METHOD     : DifferentiateD
// DESCRIPTION: Computes derivative of input vector in direction idir using
//              on basis 'dbasis'. dv, and v must both have size at least (dbasis->dim() +1)^d.
//              It is assumed that the basis is isotropic (the same in each coord
//              direction).
//
//              Note: It is best to call 'SetDxiDxD' method before entry, in case the
//                    element type requires these to be set.
// ARGUMENTS  : dv    : output arguement: derivative in idir direction of v
//              v     : input argument as a GVector
//              idir  : coordinate direction
//              dbasis: basis on which to compute derivative
// RETURNS    : TRUE on success; else FALSE
//************************************************************************************
inline GBOOL GRegElem::DifferentiateD(GVector *dv, GVector *v, GINT  idir, GNBasis *dbasis)
{
  GINT    i, Np[GDIM];
  GDOUBLE iL[GDIM];
  GMatrix *D1d[GDIM], *D1dT[GDIM];
  char   *serr = "GRegElem::DifferentiateD: ";

  for ( i=0; i<GDIM; i++ ) {
    Np  [i] = dbasis->GetOrder()+1; 
    D1d [i] = dbasis->GetDerivMatrix();
    D1dT[i] = dbasis->GetDerivMatrix(TRUE);
    iL  [i] = 2.0/L_[i];
  }

  switch ( idir ) {
    case 1:
#if defined(G_IS2D)
      // Do I2 X D1 term: 
      MTK::I2_X_D1(*D1d[0], *v, Np[0], Np[1], *dv);
#elif defined(G_IS3D)
      // Do I3 X I2 X D1 term: 
      MTK::I3_X_I2_X_D1(*D1d[0], *v, Np[0], Np[1], Np[2], *dv);
#endif
      break;
    case 2: 
#if defined(G_IS2D)
      // Do D2 X I1 term: 
      MTK::D2_X_I1(*D1dT[1], *v, Np[0], Np[1], *dv);
#elif defined(G_IS3D)
      // Do I3 X D2 X I1 term: 
      MTK::I3_X_D2_X_I1(*D1dT[1], *v, Np[0], Np[1], Np[2], *dv);
#endif
      break;
#if defined(G_IS3D)
    case 3:
      // Do D3 X I2 X I1 term: 
      MTK::D3_X_I2_X_I1(*D1d_[2], *v, Np[0], Np[1], Np[2], *dv);
      break;
#endif
  }

  if ( iL[idir-1] != 1.0 ) MTK::fvec_const_prod_rep(*dv, iL[idir-1]);


  return TRUE;

} // end of method Differentiate


#if 0
//************************************************************************************
//************************************************************************************
/**
 * METHOD     : DifferentiateWithMass
 * DESCRIPTION: Computes derivative of input vector in direction idir with mass matrix
 * ARGUMENTS  : dv  : output arguement: derivative in idir direction of v
 *              v   : input argument as a GVector
 *              idir: coordinate direction
 * RETURNS    : TRUE on success; else FALSE
 */
//************************************************************************************
GBOOL GRegElem::DifferentiateWithMass(GVector *dv, GVector *v, GVector* tmp, GINT  idir)
{
  GINT    NN, N1, N2;
  GDOUBLE   L1, L2, iL1, iL2;

  if ( v == NULL || dv == NULL ) return FALSE;

  L1   = PDISTANCE(spVertices_[1],spVertices_[0]);
  L2   = PDISTANCE(spVertices_[2],spVertices_[1]);
  if ( L1 == 0 || L2 == 0 ) return FALSE;

  iL1 = L2/2.0;
  iL2 = L1/2.0;

  N1 = Np_[0];
  N2 = Np_[1];
  NN = N1  * N2;
  if ( v->dim() != NN || dv->dim() != NN ) return FALSE;

  if ( idir == 1 ) {
    // Do Dg2 X D1 term: 
    MTK::Dg2_X_D1(*MD1,*Weights_[1], *v, *tmp, *dv);
    if ( iL1 != 1.0 ) MTK::fvec_const_prod_rep(*dv, iL1);
  }
  else if ( idir == 2 ) {
    // Do D2 X D1 term: 
    MTK::D2_X_Dg1(*D2TM,*Weights_[0], *v, *tmp, *dv);
    if ( iL2 != 1.0 ) MTK::fvec_const_prod_rep(*dv, iL2);
  }

  return TRUE;

} // end of method DifferentiateWithMass


//************************************************************************************
//************************************************************************************
/**
 * METHOD     : DifferentiateWeak
 * DESCRIPTION: Computes weak derivative of input vector in direction idir with mass matrix
 * ARGUMENTS  : dv  : output arguement: derivative in idir direction of v
 *              v   : input argument as a GVector
 *              idir: coordinate direction
 * RETURNS    : TRUE on success; else FALSE
 */
//************************************************************************************
GBOOL GRegElem::DifferentiateWeak(GVector *dv, GVector *v, GVector *tmp,GINT  idir)
{
  GINT    NN, N1, N2;
  GDOUBLE   L1, L2, iL1, iL2;

  if ( v == NULL || dv == NULL ) return FALSE;

  L1   = PDISTANCE(spVertices_[1],spVertices_[0]);
  L2   = PDISTANCE(spVertices_[2],spVertices_[1]);
  if ( L1 == 0 || L2 == 0 ) return FALSE;

  iL1 = L2/2.0;
  iL2 = L1/2.0;

  N1 = Np_[0];
  N2 = Np_[1];
  NN = N1  * N2;
  if ( v->dim() != NN || dv->dim() != NN ) return FALSE;

  if ( idir == 1 ) {
    // Do Dg2 X D1 term: 
    MTK::Dg2_X_D1(*D1TM,*Weights_[1], *v, *tmp, *dv);
    if ( iL1 != 1.0 ) MTK::fvec_const_prod_rep(*dv, -iL1);// Notice minus sign here
  }
  else if ( idir == 2 ) {
    // Do D2 X D1 term: 
    MTK::D2_X_Dg1(*MD2,*Weights_[0], *v, *tmp, *dv);
    if ( iL2 != 1.0 ) MTK::fvec_const_prod_rep(*dv, -iL2);//Notice minus sign here
  }
  return TRUE;

} // end of method DifferentiateWeak
#endif


//************************************************************************************
//************************************************************************************
// METHOD     : ComputeDealias
// DESCRIPTION: Computes quantities used in performing dealiasing, in particular, 
//              the dealiased mass matrices, for all classes (orders) of nonlinearity.
//              Method GElemB::SetDealiasBasis must have been called prior to entry.
// ARGUMENTS  : 
// RETURNS    : TRUE on success; else FALSE
//************************************************************************************
void GRegElem::ComputeDealias()
{
  GINT             i, im, j, k, DNp[3], m, n;
  GBOOL            bDealiased;
  GDOUBLE          Jac, tmp;
  GTVector<GQUAD> *dqweights[3];
  GVector          pxi, dpxi;
  char            *serr = "GRegElem::ComputeDealias: ";

  for ( j=0, bDealiased=TRUE; j<GDIM; j++ ) {
    bDealiased = bDealiased && dealias_basis_[j] != NULL;
  }

  if ( nnonlin_ == 0 ) return;  

  Jac     = 1.0 / pow(2.0,GDIM);
  for ( k=0; k<GDIM; k++ ) Jac     *= L_[k];

  for ( m=0; m<nnonlin_; m++ ) {

    DNp[2] = 1;
    for ( j=0; j<GDIM; j++ ) {
      DNp      [j] = dealias_basis_[m][j]->GetOrder() + 1;
      dqweights[j] = dealias_basis_[m][j]->GetQWeights();
    } 

    DMassMatrix_[m].Resize(DNp[0]*DNp[1]*DNp[2]);
    for ( k=0; k<DNp[2]; k++ ) {
      for ( j=0; j<DNp[1]; j++ ) {
        for ( i=0; i<DNp[0]; i++ ) {
          n   = i + j*DNp[0] + k*DNp[0]*DNp[1];
          tmp = (*dqweights[0])[i] * (*dqweights[1])[j] * Jac;
          DMassMatrix_[m][n] = GDIM == 3 ? tmp*(*dqweights[2])[k] : tmp; // include z-component, if necessary
        }
      }
    }

  }

} // end of method ComputeDealias


#if 0
//************************************************************************************
//************************************************************************************
// METHOD     : isCpoint
// DESCRIPTION: Determines whether specified vertex id identifies a 'C'-type point or
//              not. C-point is a real vertex id that does not line up with its mortar
//              endpoint. Method returns 1d mortar that doesn't line up, as as well as 
//              the mortar boundary (0,1) point that fails to match up.
// ARGUMENTS  :  jvert    : vertex index
//               imortfail: mortar index returned
//               imbdy    : which mortar bdy index fails to match up
// RETURNS    : TRUE on success; else FALSE
//************************************************************************************
inline GBOOL GRegElem::isCpoint(GINT jvert, GINT &imortfail, GINT &imbdy)
{
  GINT     j, nbdy, imort[2], ibdy[2];
  GBOOL    bCpoint;
  GFPoint *pb;
  char    *serr = "GRegElem::isCpoint: ";

  // Check all 1d mortars for this vertex:
  if        ( jvert == 0 ) {
    imort[0] = 0; imort[1] = 3;
    ibdy [0] = 0; ibdy [1] = 0 ;
  } else if ( jvert == 1 ) {
    imort[0] = 0; imort[1] = 1;
    ibdy [0] = 1; ibdy [1] = 0 ;
  } else if ( jvert == 2 ) {
    imort[0] = 1; imort[1] = 2;
    ibdy [0] = 1; ibdy [1] = 1 ;
  } else if ( jvert == 3 ) {
    imort[0] = 2; imort[1] = 3;
    ibdy [0] = 0; ibdy [1] = 1 ;
  }
  else {
    cout << serr << "invalid vertex id" << endl;
    exit(1);
  }
  
  for ( j=0, bCpoint=FALSE; j<2 && !bCpoint; j++ ) {
    pb = edge_mortars_[imort[j]].GetMortarMidpoint();
    if ( spVertices_[jvert] == *pb ) { bCpoint = TRUE; imortfail = imort[j]; imbdy = ibdy[j]; }
  }
#if 0
if ( bCpoint ) {
  cout << serr << "spVert[" << jvert << "]=" << spVertices_[jvert] << " mort_mid=" << *pb << " imortfail=" << imortfail
<< " imbdy=" << imbdy << endl;
}
#endif

  return bCpoint;
} // end of method isCpoint
#endif


//************************************************************************************
//************************************************************************************
// METHOD     : SetDxiDxD
// DESCRIPTION: Sets inverse Jaobian matrix elements (e.g., for derivatives).
//              Is intended to be used for non-'native' (member) bases, such
//              as dealiasing bases. Each time this method is called a new set
//              of matrix elements is computed and stored. In this way, sets of
//              Jacobian matrices for different bases (orders) may be used to integrate
//              equations.
//
//              To use these Jacobian matrices, one may either use the 'DifferentiateD' method
//              or get a pointer with a call to the 'GetDxiDxD' method.
// ARGUMENTS  :  dbasis   : basis to use in computing dXi_i/dX_j
// RETURNS    : none
//************************************************************************************
void GRegElem::SetDxiDxD(GNBasis *dbasis)
{
  GINT     j;
  char    *serr = "GRegElem::SetDxiDxD: ";

  if ( elemtype_ == RECT_QUAD ) {
    return;
  }

  cout << serr << "Element types other than RECT_QUAD/CUBE not currently supported" << endl;
  exit(1);


} // end of method SetDxiDxD


