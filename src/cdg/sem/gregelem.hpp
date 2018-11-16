//************************************************************************************//
// Module       : rectquad2d.hpp
// Date         : 9/14/01 (DLR)
// Copyright    : 2001-2006 Copyright University Corporation for Atmospheric
//                Research
// Description  : Encapsulates the methods and data associated with
//                a single regular rectangular 2D spectral 
//                element. The following is the ordering for the vertices (Vi),
//                segment midpoints (Mi), and faces (Fi):
// 
//
//          V7 ________________ V6
//            /|     M6       /|
//       M7  / |             / |
//          /  |    F5   M5 /  |
//         /M11|  M4   F2  /   |
//     V4 /____|__________/    | M10
//        |    |         | V5  |
//        |    |     M2  | F1  |
//        | V3 |_________|_____| V2
//    M8  | F3 / F0      |    /
//        |   /       M9 |   /
//        |  / M3   F4   |  / M1
//        | /            | /
//        |/_____________|/
//       V0       M0     V1
//
// Faces are labeled s.t. F0-F3 correspond to orientation of edges on bottom plane;
// F4 and F5 are bottom and top faces, respectively.
//
// [Note that in 2d, we use just the bottom plane.]
//                  
// Derived From : GElemB.
// Modifications:
//************************************************************************************//
#if !defined(GREGELEM_HPP)
#define GREGELEM_HPP 
#include "gelemb.hpp"
#include "gcutils.hpp"

class GRegElem : public GElemB
{
public:
                 GRegElem();
                 GRegElem(GNBasis *b1, GNBasis *b2, GNBasis *b3=NULL);
                 GRegElem(GNBasis *b[], GINT nb);
                ~GRegElem();
void             operator=(const GRegElem &);
                   
//GINT           GetOrder(const GINT  idir);
GDOUBLE          Integrate(GVector *v, GDOUBLE *multiplicity=NULL);
GDOUBLE          PIntegrate(GVector *v, GDOUBLE *multiplicity=NULL);
inline GBOOL     Differentiate(GVector *dv, GVector *v, GINT  idir);
inline GBOOL     DifferentiateD(GVector *dv, GVector *v, GINT  idir, GNBasis *dbasis);
void             SetDxiDxD(GNBasis *dbasis);

#if 0
GBOOL            DifferentiateWithMass(GVector *dv, GVector *v, GVector* tmp, GINT  idir);
GBOOL            DifferentiateWeak(GVector *dv, GVector *v, GVector* tmp, GINT  idir);
#endif

GVector         *GetJacobian();                             
GVector         *GetdXidX(GMatrix **, GINT  );             
GVector         *GetdXidX   (const GINT  i, const GINT  j );
GVector         *GetMetric  (const GINT  i, const GINT  j );
GVector         *GetWJMetric(const GINT  i, const GINT  j );
GMatrix         *GetInterpOp(GINT  idir, GBOOL Transpose=FALSE);
GMatrix         *GetInterpDeriv(GINT  idir, GBOOL Transpose=FALSE);


GBOOL            XToXi(GFPoint pX[], GFPoint pXi[], const GINT  n);
//GMatrix       *GetBasisAtXi(GINT  i, GINT  j, GMatrix *B); // native ordering


GBOOL            ComputeSpNodes();
GBOOL            ComputeLaplacian();
void             ComputeDealias();
GBOOL            SetVertices(GFPoint P[], GINT  num);
//void             SetBasis(GNBasis *b1, GNBasis *b2);
//void             SetBasis(GNBasis *b, GINT  idir);
//void             SetOrder(GINT  iorder1, GINT  iorder2);
void             SetInterpBasis(GNBasis *b[], GINT nb);
GBOOL            isCpoint(GINT ivertex, GINT &imort, GINT &ibdy);
GBOOL            Point_in_elem(GFPoint pt[], GINT np);
GBOOL            Point_in_elem(GDOUBLE x, GDOUBLE y, GDOUBLE z=0.0);
GBOOL            AnyPoint_in_elem(GFPoint p[], GINT n);

GINT             SolveFE();

protected:

private:

// Private methods:
GBOOL             Initialize(GINT  order1, GINT  order2);
void              DeleteDynamic();
GBOOL             ComputeGeo();

GDOUBLE           L_[GDIM];

// Private data:


};

//************************************************************************************
//************************************************************************************
// METHOD     : isCpoint
// DESCRIPTION: Determines whether specified vertex id identifies a 'C'-type point or
//              not. C-point is a real vertex id that does not line up with a mortar
//              node (a hanging node). Method returns mortar that doesn't line up, 
//              as as well as the mortar boundary (0,1) point that fails to match up.
// ARGUMENTS  :  jvert    : vertex index
//               imortfail: mortar index returned
//               imbdy    : which mortar bdy index fails to match up
// RETURNS    : TRUE on success; else FALSE
//************************************************************************************
inline GBOOL GRegElem::isCpoint(GINT jvert, GINT &imortfail, GINT &imbdy) {
  GINT     j, nbdy, imort[GDIM], ibdy[GDIM];
  GBOOL    bCpoint;
  GFPoint *pb;
  char  *serr = "GRegElem::isCpoint: ";

#if defined(G_IS2D)
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

  for ( j=0, bCpoint=FALSE; j<2*(GDIM-1) && !bCpoint; j++ ) {
    pb = edge_mortars_[imort[j]].GetMortarMidpoint();
    if ( spVertices_[jvert] == *pb ) { bCpoint = TRUE; imortfail = imort[j]; imbdy = ibdy[j]; }
  }
  return bCpoint;
#elif defined(G_IS3D)
  // Check all 2d mortars for this vertex:
  if        ( jvert == 0 ) {
    imort[0] = 0; imort[1] = 3; imort[2] = 4;
    ibdy [0] = 0; ibdy [1] = 0; ibdy [2] = 0; 
  } else if ( jvert == 1 ) {
    imort[0] = 0; imort[1] = 1; imort[2] = 4;
    ibdy [0] = 1; ibdy [1] = 0; ibdy [2] = 1;
  } else if ( jvert == 2 ) {
    imort[0] = 1; imort[1] = 2; imort[2] = 4;
    ibdy [0] = 1; ibdy [1] = 1; ibdy [2] = 2;
  } else if ( jvert == 3 ) {
    imort[0] = 2; imort[1] = 3; imort[2] = 4;
    ibdy [0] = 0; ibdy [1] = 1; ibdy [2] = 3;
  } else if ( jvert == 4 ) {
    imort[0] = 0; imort[1] = 3; imort[2] = 5;
    ibdy [0] = 0; ibdy [1] = 0; ibdy [2] = 0;
  } else if ( jvert == 5 ) {
    imort[0] = 0; imort[1] = 1; imort[2] = 5;
    ibdy [0] = 2; ibdy [1] = 2; ibdy [2] = 1;
  } else if ( jvert == 6 ) {
    imort[0] = 1; imort[1] = 2; imort[2] = 5;
    ibdy [0] = 2; ibdy [1] = 2; ibdy [2] = 2;
  } else if ( jvert == 7 ) {
    imort[0] = 2; imort[1] = 3; imort[2] = 5;
    ibdy [0] = 3; ibdy [1] = 2; ibdy [2] = 3;
  }
  else {
    cout << serr << "invalid vertex id" << endl;
    exit(1);
  }

  for ( j=0, bCpoint=FALSE; j<GDIM && !bCpoint; j++ ) {
    pb = face_mortars_[imort[j]].GetMortarMidpoint();
    if ( spVertices_[jvert] == *pb ) { 
      bCpoint = TRUE; imortfail = imort[j]; imbdy = ibdy[j]; 
//    cout << serr << "spVert=" << spVertices_[jvert] << " mmid=" << *pb << endl;
    }
  }
  return bCpoint;
#endif
} 



//************************************************************************************
//************************************************************************************
// METHOD     : Point_in_elem
// DESCRIPTION: Determines whether specified point is in the element (or on its bdy).
//              Note: This method valid only for regular (polygonal) hexagons/quads
// ARGUMENTS  : pt : point array
//              np : number points in array
// RETURNS    : TRUE on success (all points in element); else FALSE
//************************************************************************************
inline GBOOL GRegElem::Point_in_elem(GDOUBLE x, GDOUBLE y, GDOUBLE z)
{
  GFPoint  P;
  char    *serr = "GRegElem::Point_in_elem: ";

  P.x1 = x; P.x2 = y; P.x3 = z;

  return Point_in_elem(&P,1);

} // end of method Point_in_elem



//************************************************************************************
//************************************************************************************
// METHOD     : Point_in_elem
// DESCRIPTION: Determines whether specified point is in the element (or on its bdy).
//              Note: This method valid only for regular (polygonal) hexagons/quads
// ARGUMENTS  : pt : point array
//              np : number points in array
// RETURNS    : TRUE on success (all points in element); else FALSE
//************************************************************************************
inline GBOOL GRegElem::Point_in_elem(GFPoint pt[], GINT np) 
{
  GBOOL    bRet=TRUE;
  GINT     j, n; 
  GFPoint  N, pverts[4], R, V1, V2;
  GDOUBLE  fact, RdotN;
  char    *serr = "GRegElem::Point_in_elem: ";

#if defined(G_IS2D)
  bRet = GCUtils::Point_in_poly(spVertices_, nVertices_, pt, np); // bottom face
#elif defined(G_IS3D)
  // Check if in volume:
  for ( n=0; n<np && bRet; n++ ) {
    // For each face, find R = p - V0, and inner-pointing normal to plane(
    // compute N = V1 X V2 / ( |V1||V2| ) ). Take R.N; must be >=0 for each face:
    
    // bottom face:
    for ( j=0; j<GDIM; j++ ) R [j] = pt[n][j] - spVertices_[0][j];
    GCUtils::PlaneNormal(N, spVertices_[0], spVertices_[1], spVertices_[3]);
    for ( j=0, RdotN=0.0; j<GDIM; j++ ) RdotN += N[j] * R[j];
//  bRet = bRet && ( fabs(RdotN) < R.GetBracket() && RdotN > 0 );
    bRet = bRet && (RdotN >= -R.GetBracket());
   
    // top face:
    for ( j=0; j<GDIM; j++ ) R [j] = pt[n][j] - spVertices_[4][j];
    GCUtils::PlaneNormal(N, spVertices_[4], spVertices_[7], spVertices_[5]);
    for ( j=0, RdotN=0.0; j<GDIM; j++ ) RdotN += N[j] * R[j];
    bRet = bRet && (RdotN >= -spVertices_[0].GetBracket() );


    // front face:
    for ( j=0; j<GDIM; j++ ) R [j] = pt[n][j] - spVertices_[0][j];
    GCUtils::PlaneNormal(N, spVertices_[0], spVertices_[4], spVertices_[1]);
    for ( j=0, RdotN=0.0; j<GDIM; j++ ) RdotN += N[j] * R[j];
    bRet = bRet && (RdotN >= -R.GetBracket() );

    // right face:
    for ( j=0; j<GDIM; j++ ) R [j] = pt[n][j] - spVertices_[1][j];
    GCUtils::PlaneNormal(N, spVertices_[1], spVertices_[5], spVertices_[2]);
    for ( j=0, RdotN=0.0; j<GDIM; j++ ) RdotN += N[j] * R[j];
    bRet = bRet && (RdotN >= -R.GetBracket() );

    // back face:
    for ( j=0; j<GDIM; j++ ) R [j] = pt[n][j] - spVertices_[3][j];
    GCUtils::PlaneNormal(N, spVertices_[3], spVertices_[2], spVertices_[7]);
    for ( j=0, RdotN=0.0; j<GDIM; j++ ) RdotN += N[j] * R[j];
    bRet = bRet && (RdotN >= -R.GetBracket() );

    // left face:
    for ( j=0; j<GDIM; j++ ) R [j] = pt[n][j] - spVertices_[0][j];
    GCUtils::PlaneNormal(N, spVertices_[0], spVertices_[3], spVertices_[4]);
    for ( j=0, RdotN=0.0; j<GDIM; j++ ) RdotN += N[j] * R[j];
    bRet = bRet && (RdotN >= -R.GetBracket() );
  }
#endif

  return bRet;
} // end of method Point_in_elem


//************************************************************************************
//************************************************************************************
// METHOD     : AnyPoint_in_elem
// DESCRIPTION: Determines whether any of the specified points is in the element
//              (or on its bdy).  Note: This method valid only for regular (polygonal) 
//              hexagons/quads
// ARGUMENTS  : pt : point array
//              np : number points in array
// RETURNS    : TRUE on success (at least one point in element); else FALSE
//************************************************************************************
inline GBOOL GRegElem::AnyPoint_in_elem(GFPoint pt[], GINT np) 
{
  GBOOL    bRet, bFace;
  GINT     j, n; 
  GFPoint  N, pverts[4], R, V1, V2;
  GDOUBLE  fact, RdotN;
  char    *serr = "GRegElem::AnyPoint_in_elem: ";

  bRet = GCUtils::AnyPoint_in_poly(spVertices_  , 4, pt, np); // bottom face

#if defined(G_IS3D)
  // Check if in volume:
  for ( n=0, bRet=FALSE; n<np && !bRet; n++ ) {
    // For each face, find R = p - V0, and inner-pointing normal to plane(
    // compute N = V1 X V2 / ( |V1||V2| ) ). Take R.N; must be >=0 for each face:
    bFace = TRUE;
    // bottom face:
    for ( j=0; j<GDIM; j++ ) R [j] = pt[n][j] - spVertices_[0][j];
    GCUtils::PlaneNormal(N, spVertices_[0], spVertices_[1], spVertices_[3]);
    for ( j=0, RdotN=0.0; j<GDIM; j++ ) RdotN += N[j] * R[j];
//  bRet = bRet && ( fabs(RdotN) < R.GetBracket() && RdotN > 0 );
    bFace = bFace && (RdotN >= 0);
   
    // top face:
    for ( j=0; j<GDIM; j++ ) R [j] = pt[n][j] - spVertices_[5][j];
    GCUtils::PlaneNormal(N, spVertices_[5], spVertices_[4], spVertices_[6]);
    for ( j=0, RdotN=0.0; j<GDIM; j++ ) RdotN += N[j] * R[j];
    bFace = bFace && (RdotN >= -R.GetBracket());


    // front face:
    for ( j=0; j<GDIM; j++ ) R [j] = pt[n][j] - spVertices_[1][j];
    GCUtils::PlaneNormal(N, spVertices_[1], spVertices_[0], spVertices_[5]);
    for ( j=0, RdotN=0.0; j<GDIM; j++ ) RdotN += N[j] * R[j];
    bFace = bFace && (RdotN >= R.GetBracket());

    // right face:
    for ( j=0; j<GDIM; j++ ) R [j] = pt[n][j] - spVertices_[1][j];
    GCUtils::PlaneNormal(N, spVertices_[1], spVertices_[2], spVertices_[5]);
    for ( j=0, RdotN=0.0; j<GDIM; j++ ) RdotN += N[j] * R[j];
    bFace = bFace && (RdotN >= R.GetBracket());

    // back face:
    for ( j=0; j<GDIM; j++ ) R [j] = pt[n][j] - spVertices_[2][j];
    GCUtils::PlaneNormal(N, spVertices_[2], spVertices_[6], spVertices_[3]);
    for ( j=0, RdotN=0.0; j<GDIM; j++ ) RdotN += N[j] * R[j];
    bFace = bFace && (RdotN >= R.GetBracket());

    // left face:
    for ( j=0; j<GDIM; j++ ) R [j] = pt[n][j] - spVertices_[0][j];
    GCUtils::PlaneNormal(N, spVertices_[0], spVertices_[3], spVertices_[4]);
    for ( j=0, RdotN=0.0; j<GDIM; j++ ) RdotN += N[j] * R[j];
    bFace = bFace && (RdotN >= R.GetBracket());
    
    bRet = bRet || bFace;
  }
#endif

  return bRet;
} // end of method AnyPoint_in_elem

#endif
