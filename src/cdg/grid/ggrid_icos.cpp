//==================================================================================
// Description  : Object defining a (global) icosahedral grid, that
//                uses gnomonic projections to locate element vertices.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : GGrid.
//==================================================================================

#include <cstdlib>
#include <memory>
#include <cmath>
//#include "omp.h"  // Are we calling API functions ?
#include "geoflow.hpp"
#include "gspecbdy_factory.hpp"
#include "gelem_base.hpp"
#include "ggrid_icos.hpp"
#include "gtpoint.hpp"

using namespace std;

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Instantiate for 2d grid
// ARGS   : ptree: main property tree
//          b     : vector of basis pointers, of size at least ndim=2.Determies 
//                  dimensionality
//          comm  : communicator
// RETURNS: none
//**********************************************************************************
GGridIcos::GGridIcos(const geoflow::tbox::PropertyTree &ptree, GTVector<GNBasis<GCTYPE,GFTYPE>*> &b, GC_COMM &comm)
:          GGrid(ptree, b, comm),
ilevel_                      (0),
ndim_                     (GDIM),
radiusi_                   (0.0),
radiuso_                   (0.0), // don't need this one in 2d
gdd_                   (NULLPTR),
lshapefcn_             (NULLPTR)
{
  assert(b.size() == GDIM && "Basis has incorrect dimensionality");
  
  GString gname   = ptree.getValue<GString>("grid_type");
  assert(gname == "grid_icos" || gname == "grid_sphere");
  geoflow::tbox::PropertyTree gridptree = ptree.getPropertyTree(gname);

  gbasis_.resize(b.size());
  gbasis_ = b;
  lshapefcn_ = new GShapeFcn_linear();
  ilevel_  = gridptree.getValue<GINT>("ilevel");
  
  if ( ndim_ == 2 ) {
    assert(GDIM == 2 && "GDIM must be 2");
    radiusi_ = gridptree.getValue<GFTYPE>("radius");
    init2d();
  }
  else if ( ndim_ == 3 ) {
    assert(GDIM == 3 && "GDIM must be 3");
    std::vector<GINT> ne(3);
    radiusi_   = gridptree.getValue<GFTYPE>("radiusi");
    radiuso_   = gridptree.getValue<GFTYPE>("radiuso");
    nradelem_  = gridptree.getValue<GINT>("num_radial_elems");
    init3d();
  }
  else {
    assert(FALSE && "Invalid dimensionality");
  }

} // end of constructor method (1)



//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
GGridIcos::~GGridIcos()
{
  if ( lshapefcn_ != NULLPTR ) delete lshapefcn_;
  if ( gdd_       != NULLPTR ) delete gdd_;
} // end, destructor


//**********************************************************************************
//**********************************************************************************
// METHOD :  << operator method (1)
// DESC   : output stream operator
// ARGS   :
// RETURNS: ostream &
//**********************************************************************************
std::ostream &operator<<(std::ostream &str, GGridIcos &e)
{
  
  str << " radiusi: " << e.radiusi_;
  str << " radiusi: " << e.radiuso_;
  str << " level  : " << e.ilevel_;
  str << std::endl << " Centroids: " ;
  for ( GSIZET i=0; i<e.ftcentroids_.size(); i++ )
     str << (e.ftcentroids_[i]) << " " ;
  str << std::endl;

  return str;
} // end of operator <<


//**********************************************************************************
//**********************************************************************************
// METHOD : set_partitioner
// DESC   : Set domain decomposition object
// ARGS   : GDD_base pointer
// RETURNS: none
//**********************************************************************************
void GGridIcos::set_partitioner(GDD_base *gdd)
{

  gdd_ = gdd;

} // end of method set_partitioner


//**********************************************************************************
//**********************************************************************************
// METHOD : init2d
// DESC   : Initialize base state/base icosahedron
// ARGS   : none.
// RETURNS: none.
//**********************************************************************************
void GGridIcos::init2d()
{
  GString serr = "GridIcos::init2d: ";

  GFTYPE phi = (1.0+sqrt(5.0))/2.0;  // Golden ratio


  // Base vertex list:
  fv0_.resize(12,3);
#if 0
  // Following give orientation s.t. edge lies at the top:
  fv0_(0 ,0) = 0  ; fv0_(0 ,1) = 1  ; fv0_(0 ,2) = phi;
  fv0_(1 ,0) = 0  ; fv0_(1 ,1) =-1  ; fv0_(1 ,2) = phi;
  fv0_(2 ,0) = 0  ; fv0_(2 ,1) = 1  ; fv0_(2 ,2) =-phi;
  fv0_(3 ,0) = 0  ; fv0_(3 ,1) =-1  ; fv0_(3 ,2) =-phi;
  fv0_(4 ,0) = 1  ; fv0_(4 ,1) = phi; fv0_(4 ,2) = 0;
  fv0_(5 ,0) =-1  ; fv0_(5 ,1) = phi; fv0_(5 ,2) = 0;
  fv0_(6 ,0) = 1  ; fv0_(6 ,1) =-phi; fv0_(6 ,2) = 0;
  fv0_(7 ,0) =-1  ; fv0_(7 ,1) =-phi; fv0_(7 ,2) = 0;
  fv0_(8 ,0) = phi; fv0_(8 ,1) = 0  ; fv0_(8 ,2) = 1;
  fv0_(9 ,0) =-phi; fv0_(9 ,1) = 0  ; fv0_(9 ,2) = 1;
  fv0_(10,0) = phi; fv0_(10,1) = 0  ; fv0_(10,2) =-1;
  fv0_(11,0) =-phi; fv0_(11,1) = 0  ; fv0_(11,2) =-1;
#endif

// Following give orientation s.t. vertex is at top:
fv0_(0 ,0) =  0.000000000000000; fv0_(0 ,1) =  0.000000000000000; fv0_(0 ,2) =  1.000000000000000;
fv0_(1 ,0) =  0.894427190999916; fv0_(1 ,1) =  0.000000000000000; fv0_(1 ,2) =  0.447213595499958;
fv0_(2 ,0) = -0.894427190999916; fv0_(2 ,1) =  0.000000000000000; fv0_(2 ,2) = -0.447213595499958;
fv0_(3 ,0) =  0.000000000000000; fv0_(3 ,1) =  0.000000000000000; fv0_(3 ,2) = -1.000000000000000;
fv0_(4 ,0) = -0.723606797749979; fv0_(4 ,1) =  0.525731112119134; fv0_(4 ,2) =  0.447213595499958;
fv0_(5 ,0) = -0.723606797749979; fv0_(5 ,1) = -0.525731112119134; fv0_(5 ,2) =  0.447213595499958;
fv0_(6 ,0) =  0.723606797749979; fv0_(6 ,1) =  0.525731112119134; fv0_(6 ,2) = -0.447213595499958;
fv0_(7 ,0) =  0.723606797749979; fv0_(7 ,1) = -0.525731112119134; fv0_(7 ,2) = -0.447213595499958;
fv0_(8 ,0) =  0.276393202250021; fv0_(8 ,1) =  0.850650808352040; fv0_(8 ,2) =  0.447213595499958;
fv0_(9 ,0) =  0.276393202250021; fv0_(9 ,1) = -0.850650808352040; fv0_(9 ,2) =  0.447213595499958;
fv0_(10,0) = -0.276393202250021; fv0_(10,1) =  0.850650808352040; fv0_(10,2) = -0.447213595499958;
fv0_(11,0) = -0.276393202250021; fv0_(11,1) = -0.850650808352040; fv0_(11,2) = -0.447213595499958;

  
  // Normalize vertices to be at 1.0:
  GFTYPE fact = 1.0/sqrt(phi*phi + 1.0);
  fv0_ *= fact;

  // Points in verts array that make up 
  // each face of base icosahedron:
  ifv0_.resize(20,3);
  ifv0_(0 ,0) = 0 ; ifv0_(0 ,1) = 1 ; ifv0_(0 ,2) = 8 ;
  ifv0_(1 ,0) = 0 ; ifv0_(1 ,1) = 8 ; ifv0_(1 ,2) = 4 ;
  ifv0_(2 ,0) = 0 ; ifv0_(2 ,1) = 4 ; ifv0_(2 ,2) = 5 ;
  ifv0_(3 ,0) = 0 ; ifv0_(3 ,1) = 5 ; ifv0_(3 ,2) = 9 ;
  ifv0_(4 ,0) = 0 ; ifv0_(4 ,1) = 9 ; ifv0_(4 ,2) = 1 ;
  ifv0_(5 ,0) = 1 ; ifv0_(5 ,1) = 6 ; ifv0_(5 ,2) = 8 ;
  ifv0_(6 ,0) = 8 ; ifv0_(6 ,1) = 6 ; ifv0_(6 ,2) = 10;
  ifv0_(7 ,0) = 8 ; ifv0_(7 ,1) = 10; ifv0_(7 ,2) = 4 ;
  ifv0_(8 ,0) = 4 ; ifv0_(8 ,1) = 10; ifv0_(8 ,2) = 2 ;
  ifv0_(9 ,0) = 4 ; ifv0_(9 ,1) = 2 ; ifv0_(9 ,2) = 5 ;
  ifv0_(10,0) = 5 ; ifv0_(10,1) = 2 ; ifv0_(10,2) = 11;
  ifv0_(11,0) = 5 ; ifv0_(11,1) = 11; ifv0_(11,2) = 9 ;
  ifv0_(12,0) = 9 ; ifv0_(12,1) = 11; ifv0_(12,2) = 7 ;
  ifv0_(13,0) = 9 ; ifv0_(13,1) = 7 ; ifv0_(13,2) = 1 ;
  ifv0_(14,0) = 1 ; ifv0_(14,1) = 7 ; ifv0_(14,2) = 6 ;
  ifv0_(15,0) = 3 ; ifv0_(15,1) = 6 ; ifv0_(15,2) = 7 ;
  ifv0_(16,0) = 3 ; ifv0_(16,1) = 7 ; ifv0_(16,2) = 11;
  ifv0_(17,0) = 3 ; ifv0_(17,1) = 11; ifv0_(17,2) = 2 ;
  ifv0_(18,0) = 3 ; ifv0_(18,1) = 2 ; ifv0_(18,2) = 10;
  ifv0_(19,0) = 3 ; ifv0_(19,1) = 10; ifv0_(19,2) = 6 ;

  // Copy base data to new structure:
  GTPoint<GFTYPE> pt;
  tbase_.resize(ifv0_.size(1));
  for ( GSIZET j=0; j<tbase_.size(); j++ ) tbase_[j].resize(3);
  for ( GSIZET i=0; i<ifv0_.size(1); i++ ) { // for all triangles:
    for ( GSIZET j=0; j<ifv0_.size(2); j++ ) { // for each vertex:
      for ( GSIZET k=0; k<3; k++ ) pt[k] = fv0_(ifv0_(i,j),k);
      *tbase_[i].v[j] = pt;
    }
  }

  lagrefine();

} // end of method init2d


//**********************************************************************************
//**********************************************************************************
// METHOD : init3d
// DESC   : Initialize for 3d elements
// ARGS   : none.
// RETURNS: none.
//**********************************************************************************
void GGridIcos::init3d()
{
  GString serr = "GridIcos::init3d: ";

  init2d();

} // end, method init3d


//**********************************************************************************
//**********************************************************************************
// METHOD : lagrefine
// DESC   : Use 'Lagrangian polynomial'-like placement of 'nodes' in 
//          base face/triangle to refine, before doing projection of 
//          vertices. An alternative might be, say, a self-similar 
//          (recursive) refinement of every triangle into 4 triangles.
// ARGS   : none.
// RETURNS: none.
//**********************************************************************************
void GGridIcos::lagrefine()
{
  GString serr = "GridIcos::lagrefine: ";
   
  GLLONG ibeg, j, l, m, n, t;

  // Re-dimension mesh points to be 3d: Expect
  // 20 * (ileve+1)^2, and 20 * (ilevel+1)^2 * 3 quads
  tmesh_.resize(20*(ilevel_*(ilevel_+2)+1)); // refined triangular mesh
  for ( j=0; j<tmesh_.size(); j++ ) tmesh_[j].resize(3);

  GTPoint<GFTYPE> a(3), b(3), c(3); // base vertices

  n = 0;
  GTVector<GTPoint<GFTYPE>> R0(2*(ilevel_+1)-1), R1(2*(ilevel_+1));
  GTVector<GTPoint<GFTYPE>> Rz(2*(ilevel_+1)+1); // interleave R0, R1

#pragma omp parallel private (ibeg,l,m,t,a,b,c,R0,R1,Rz) reduction(+: n)
  R0.resize(2*(ilevel_+1)-1);
  R1.resize(2*(ilevel_+1));
  Rz.resize(2*(ilevel_+1)+1);

  // Do refinement of base mesh triandles:
#pragma omp for
  for ( t=0; t<tbase_.size(); t++ ) { // for each base triangle 
    a = tbase_[t].v1; b = tbase_[t].v2; c = tbase_[t].v3;
    for ( l=0; l<ilevel_+1; l++ ) { // for each triangle 'row'
      lagvert(a,b,c,l   ,R0); // get previous row of points
      lagvert(a,b,c,l+1 ,R1); // get current row of points
      interleave(R0, R1, l, Rz); // zig-zag from R1 to R0, to R1, etc
      for ( m=0; m<2*l+1; m++ ) { // loop over all tri's on this row
        ibeg = m;
        tmesh_[n].v1 = Rz[ibeg];
        tmesh_[n].v2 = Rz[ibeg+1];
        tmesh_[n].v3 = Rz[ibeg+2];
        n++;
      }
    } // end, l-loop
  } // end, t-loop

  
  // Project all vertices to unit sphere:
  project2sphere(tmesh_,1.0);

  // Order triangles (set iup_ flags):
  order_triangles(tmesh_);

  // Compute centroids of all triangles:
  ftcentroids_.clear();
  GFTYPE fact = 1.0/3.0;
  for ( j=0; j<tmesh_.size(); j++ ) { // for each triangle
    a =  tmesh_[j].v1 + tmesh_[j].v2;
    a += tmesh_[j].v3;
    a *= fact;
    ftcentroids_.push_back(a);
  }

} // end of method lagrefine


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
void
GGridIcos::interleave(GTVector<GTPoint<GFTYPE>> &R0, GTVector<GTPoint<GFTYPE>> &R1,
                   GINT I, GTVector<GTPoint<GFTYPE>> &Rz)
{
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
void
GGridIcos::lagvert(GTPoint<GFTYPE>&a, GTPoint<GFTYPE> &b, GTPoint<GFTYPE> &c,
                   GINT I, GTVector<GTPoint<GFTYPE>> &R)
{
  GString serr = "GridIcos::lagvert: ";

  GFTYPE xI, xJ;

  GTPoint<GFTYPE> rL(3);
  GTPoint<GFTYPE> rR(3);

  R.resizem(I+1);

  GFTYPE fact = 1.0/static_cast<GFTYPE>(ilevel_+1);

  // Build 'rail' points on L and R:
  xI = static_cast<GFTYPE>(I);
  rL = a + ( (b - a) * (xI * fact) );
  rR = a + ( (c - a) * (xI * fact) );

  // Compute R vertices based on refinement indices:
  fact = I > 0 ? 1.0/static_cast<GFTYPE>(I) : 1.0;
  for ( GSIZET j=0; j<I+1; j++ ) {
    xJ = static_cast<GFTYPE>(j);
    R[j] = rL + (rR - rL)*(xJ*fact);
  }

} // end of method lagvert


//**********************************************************************************
//**********************************************************************************
// METHOD : do_elems (1)
// DESC   : Public entry point for grid computation
// ARGS   : 
//          rank: MPI rank or partition id
// RETURNS: none.
//**********************************************************************************
void GGridIcos::do_elems()
{
  if ( ndim_ == 2 ) do_elems2d(irank_);
  if ( ndim_ == 3 ) do_elems3d(irank_);

} // end, method do_elems (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : do_elems (2)
// DESC   : Public entry point for grid element computation, for restart
// ARGS   : p     : matrix of size the number of elements X GDIM containing 
//                  the poly expansion order in each direction
//          xnodes: vector of GDIM vectors containing Cartesian coords of elements
//                  for each node point
// RETURNS: none.
//**********************************************************************************
void GGridIcos::do_elems(GTMatrix<GINT> &p,
                        GTVector<GTVector<GFTYPE>> &xnodes)
{
  if ( ndim_ == 2 ) do_elems2d(p, xnodes);
  if ( ndim_ == 3 ) do_elems3d(p, xnodes);

} // end, method do_elems (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : do_elems2d (1)
// DESC   : Build 2d elemental grid on base mesh
// ARGS   : 
//          irank: MPI rank or partition id
// RETURNS: none.
//**********************************************************************************
void GGridIcos::do_elems2d(GINT irank)
{
  GString           serr = "GridIcos::do_elems2d (1): ";
  GFTYPE            fact, xlatc, xlongc;
  GTVector<GFPoint> cverts(4), gverts(4), tverts(4);
  GTPoint<GFTYPE>   c(3), ct(3), v1(3), v2(3), v3(3); // 3d points
  GElem_base        *pelem;
  GTVector<GINT>    iind;
  GTVector<GINT>    I(1);
  GTVector<GFTYPE>  Ni;
  GTVector<GTVector<GFTYPE>>  *xNodes;
  GTVector<GTVector<GFTYPE>*> *xiNodes;
  GTVector<GTVector<GFTYPE>>   xgtmp(3);

  // Do eveything on unit sphere, then project to radiusi_
  // as a final step.

  assert(gbasis_.size()>0 && "Must set basis first");

  if ( gdd_ == NULLPTR ) gdd_ = new GDD_base(nprocs_);

  // Resize points to appropriate size:
  for ( GSIZET j=0; j<tmesh_.size(); j++ ) tmesh_[j].resize(3);
  for ( GSIZET j=0; j<4; j++ ) {
    cverts[j].resize(3); // is a 3d point
    gverts[j].resize(2); // is only a 2d point
    tverts[j].resize(2); // is only a 2d point
  }

  // Get indirection indices to create elements
  // only for this task:
  gdd_->doDD(ftcentroids_, irank, iind);

  GTVector<GSIZET> isort;


  // When setting elements, must first construct Cartesian
  // coordinates at interior node points. This is done in
  // following steps:
  //   (0) find element vertices from triangular mesh
  //   (1) find centroid of element
  //   (2) use centroid to find gnomonic (2d) vertices of element
  //   (3) construct interior node points from gnomonic vertices and basis
  //   (4) transform inter node coords back to 3D Cartesian space from gnomonic space
  //   (5) project the node coords to sphere in Cart. coords 
  fact = 1.0/3.0;
  GSIZET i;
  GSIZET nfnodes;   // no. face nodes
  GSIZET icurr = 0; // current global index
  GSIZET fcurr = 0; // current global face index
  // For each triangle in base mesh owned by this rank...
  for ( GSIZET n=0; n<iind.size(); n++ ) { 
    i = iind[n];
    v1 = *tmesh_[i].v[0];
    v2 = *tmesh_[i].v[1];
    v3 = *tmesh_[i].v[2];       
    ct = (v1 + (v2 + v3)) * fact;  // triangle centroid; don't overwrite later on
    // Compute element vertices:
    // NOTE: Is this ordering compatible with shape functions 
    // (placement of +/-1 with physical point)?
    for ( GSIZET j=0; j<3; j++ ) { // 3 new elements for each triangle
      pelem = new GElem_base(GE_2DEMBEDDED, gbasis_);
      cverts[0] = v1; cverts[1] = (v1+v2)*0.5; cverts[2] = ct; cverts[3] = (v1+v3)*0.5;
#if 0
      order_latlong2d(cverts);
#else
      if ( iup_[i] == 1 ) {
      switch (j) {
      case 0: 
      cverts[0] = v1; cverts[1] = (v1+v2)*0.5; cverts[2] = ct; cverts[3] = (v1+v3)*0.5; break;
      case 1: 
      cverts[0] = (v1+v2)*0.5; cverts[1] = v2; cverts[2] = (v2+v3)*0.5; cverts[3] = ct; break;
      case 2: 
      cverts[0] = ct; cverts[1] = (v2+v3)*0.5; cverts[2] = v3; cverts[3] = (v1+v3)*0.5; break;
      }
      }
      else {
      switch (j) {
      case 0: 
      cverts[0] = v1; cverts[1] = (v1+v2)*0.5; cverts[2] = ct; cverts[3] = (v1+v3)*0.5; break;
      case 1: 
      cverts[0] = ct; cverts[1] = (v1+v2)*0.5; cverts[2] = v2; cverts[3] = (v2+v3)*0.5; break;
      case 2: 
      cverts[0] = (v1+v3)*0.5; cverts[1] = ct; cverts[2] = (v2+v3)*0.5; cverts[3] = v3; break;
      }
      }
#endif
      
      xNodes  = &pelem->xNodes();  // node spatial data
      xiNodes = &pelem->xiNodes(); // node ref interval data
      Ni.resize(pelem->nnodes());

      project2sphere(cverts, 1.0); // project verts to unit sphere     
      c = (cverts[0] + cverts[1] + cverts[2] + cverts[3])*0.25; // elem centroid
      xlatc  = asin(c.x3)         ; // reference lat/long
      xlongc = atan2(c.x2,c.x1);
      xlongc = xlongc < 0.0 ? 2*PI+xlongc : xlongc;

      cart2gnomonic(cverts, 1.0, xlatc, xlongc, tverts); // gnomonic vertices of quads
      reorderverts2d(tverts, isort, gverts); // reorder vertices consistenet with shape fcn
      for ( GSIZET l=0; l<2; l++ ) { // loop over available gnomonic coords
        xgtmp[l].resizem(pelem->nnodes());
        xgtmp[l] = 0.0;
        (*xNodes)[l] = 0.0;
        for ( GSIZET m=0; m<4; m++ ) { // loop over gnomonic vertices
          I[0] = m;
          lshapefcn_->Ni(I, *xiNodes, Ni);
          xgtmp[l] += (Ni * (gverts[m][l]*0.25)); // gnomonic node coordinate
        }
      }
      gnomonic2cart(xgtmp, 1.0, xlatc, xlongc, *xNodes); //
      project2sphere(*xNodes, radiusi_);
      pelem->init(*xNodes);
      gelems_.push_back(pelem);

      nfnodes = 0;
      for ( GSIZET j=0; j<gelems_[n]->nfaces(); j++ )  // get # face nodes
        nfnodes += gelems_[n]->face_indices(j).size();
      pelem->igbeg() = icurr;      // beginning global index
      pelem->igend() = icurr+pelem->nnodes()-1; // end global index
      pelem->ifbeg() = fcurr;
      pelem->ifend() = fcurr+nfnodes-1; // end global face index
      icurr += pelem->nnodes();
      fcurr += nfnodes;

    } // end of element loop for this triangle
  } // end of triangle base mesh loop

} // end of method do_elems2d (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : do_elems3d (1)
// DESC   : Build 3d elemental grid. It's assumed that init3d has been called
//          prior to entry, and that the icos base grid has been set.
// ARGS   : 
//          irank: MPI rank or partition id
// RETURNS: none.
//**********************************************************************************
void GGridIcos::do_elems3d(GINT irank)
{
  GString           serr = "GridIcos::do_elems3d (1): ";
  GSIZET            nxy;
  GFTYPE            fact, r0, rdelta, xlatc, xlongc;
  GTVector<GFPoint> cverts(4), gverts(4), tverts(4);
  GTPoint<GFTYPE>   c(3), ct(3), v1(3), v2(3), v3(3); // 3d points
  GElem_base        *pelem;
  GElem_base        *pelem2d;
  GTVector<GINT>    iind;
  GTVector<GINT>    I(1);
  GTVector<GFTYPE>  Ni;
  GTVector<GTVector<GFTYPE>>  *xNodes;
  GTVector<GTVector<GFTYPE>>   xNodes2d(2);
  GTVector<GTVector<GFTYPE>*>  xiNodes2d(2);
  GTVector<GFTYPE>            *xiNodesr;
  GTVector<GTVector<GFTYPE>>   xgtmp(3);

  // Do eveything on unit sphere, then project to radiusi_
  // as a final step.

  assert(gbasis_.size()>0 && "Must set basis first");

  if ( gdd_ == NULLPTR ) gdd_ = new GDD_base(nprocs_);

  // Resize points to appropriate size:
  for ( GSIZET j=0; j<tmesh_.size(); j++ ) tmesh_[j].resize(3);
  for ( GSIZET j=0; j<4; j++ ) {
    cverts[j].resize(3); // is a 3d point
    gverts[j].resize(2); // is only a 2d point
    tverts[j].resize(2); // is only a 2d point
  }

  // Get indirection indices to create elements
  // only for this task:
  gdd_->doDD(ftcentroids_, irank, iind);

  GTVector<GSIZET> isort;

  // Make radial dimension of elements the same:
  rdelta = (radiuso_ - radiusi_)/nradelem_;


  // When setting elements, must first construct Cartesian
  // coordinates at interior node points. This is done in
  // following steps:
  //   (0) find element vertices from triangular mesh
  //   (1) find centroid of element
  //   (2) use centroid to find gnomonic (2d) vertices of element
  //   (3) construct interior node points from gnomonic vertices and basis
  //   (4) transform inter node coords back to 3D Cartesian space from gnomonic space
  //   (5) project the node coords to sphere in sph. coords 
  //   (6) extend 'patch' in radial direction for all nodes in each
  //       radial element, transforming each element to Cartesian coords
  fact = 1.0/3.0;
  GSIZET i;
  GSIZET nfnodes;   // no. face nodes
  GSIZET icurr = 0; // current global index
  GSIZET fcurr = 0; // current global face index
  // For each triangle in base mesh owned by this rank...
  for ( GSIZET n=0; n<iind.size(); n++ ) { 
    i = iind[n];
    v1 = *tmesh_[i].v[0];
    v2 = *tmesh_[i].v[1];
    v3 = *tmesh_[i].v[2];       
    ct = (v1 + (v2 + v3)) * fact;  // triangle centroid; don't overwrite later on
    // Compute element vertices:
    // NOTE: Is this ordering compatible with shape functions 
    // (placement of +/-1 with physical point)?
    for ( GSIZET j=0; j<3; j++ ) { // 3 new elements for each triangle
      pelem2d = new GElem_base(GE_2DEMBEDDED, gbasis_);
      cverts[0] = v1; cverts[1] = (v1+v2)*0.5; cverts[2] = ct; cverts[3] = (v1+v3)*0.5;
      if ( iup_[i] == 1 ) {
      switch (j) {
      case 0: 
      cverts[0] = v1; cverts[1] = (v1+v2)*0.5; cverts[2] = ct; cverts[3] = (v1+v3)*0.5; break;
      case 1: 
      cverts[0] = (v1+v2)*0.5; cverts[1] = v2; cverts[2] = (v2+v3)*0.5; cverts[3] = ct; break;
      case 2: 
      cverts[0] = ct; cverts[1] = (v2+v3)*0.5; cverts[2] = v3; cverts[3] = (v1+v3)*0.5; break;
      }
      }
      else {
      switch (j) {
      case 0: 
      cverts[0] = v1; cverts[1] = (v1+v2)*0.5; cverts[2] = ct; cverts[3] = (v1+v3)*0.5; break;
      case 1: 
      cverts[0] = ct; cverts[1] = (v1+v2)*0.5; cverts[2] = v2; cverts[3] = (v2+v3)*0.5; break;
      case 2: 
      cverts[0] = (v1+v3)*0.5; cverts[1] = ct; cverts[2] = (v2+v3)*0.5; cverts[3] = v3; break;
      }
      }
      
      nxy     = pelem2d->nnodes();
      Ni.resize(nxy);
      for ( GINT m=0; m<2; m++ ) {
        xiNodes2d[m] = &pelem2d->xiNodes(m);
        xNodes2d [m].resizem(nxy);
      }

      project2sphere(cverts, 1.0); // project verts to unit sphere     
      c = (cverts[0] + cverts[1] + cverts[2] + cverts[3])*0.25; // elem centroid
      xlatc  = asin(c.x3)         ; // reference lat/long
      xlongc = atan2(c.x2,c.x1);
      xlongc = xlongc < 0.0 ? 2*PI+xlongc : xlongc;

      cart2gnomonic(cverts, 1.0, xlatc, xlongc, tverts); // gnomonic vertices of quads
      reorderverts2d(tverts, isort, gverts); // reorder vertices consistenet with shape fcn
      for ( GINT l=0; l<2; l++ ) { // loop over available gnomonic coords
        xgtmp[l].resizem(nxy);
        xgtmp[l] = 0.0;
        xNodes2d[l] = 0.0;
        for ( GSIZET m=0; m<4; m++ ) { // loop over gnomonic vertices
          I[0] = m;
          lshapefcn_->Ni(I, xiNodes2d, Ni);
          xgtmp[l] += (Ni * (gverts[m][l]*0.25)); // gnomonic node coordinate
        }
      }

      // Convert 2d plane back to Cart coords and project to 
      // surface of inner sphere:
      gnomonic2cart(xgtmp, 1.0, xlatc, xlongc, xNodes2d); 
      project2sphere(xNodes2d, radiusi_);
      xyz2spherical(xNodes2d);

      // Loop over radial elements and build all elements 
      // based on this patch:
      for ( GSIZET e=0; e<nradelem_; e++ ) {
        pelem = new GElem_base(GE_DEFORMED, gbasis_);
        xiNodesr = &pelem->xiNodes(2); // get radial reference nodes
        xNodes   = &pelem->xNodes();  // node spatial data
        r0       = radiusi_ + e*rdelta;
        for ( GSIZET m=0; m<pelem->size(2); m++ ) { // for each radial node
          for ( GSIZET n=0; n<nxy; n++ ) { // find sph coords for each horizontal node
            (*xNodes)[0][n+m*nxy] =  0.5*rdelta*((*xiNodesr)[m] + 1.0);
            (*xNodes)[1][n+m*nxy] =  xNodes2d[1][n];
            (*xNodes)[2][n+m*nxy] =  xNodes2d[2][n];
          }
        }
        spherical2xyz(*xNodes);

        pelem->init(*xNodes);
        gelems_.push_back(pelem);

        nfnodes = 0;
        for ( GSIZET j=0; j<gelems_[e]->nfaces(); j++ )  // get # face nodes
          nfnodes += gelems_[e]->face_indices(j).size();
        pelem->igbeg() = icurr;      // beginning global index
        pelem->igend() = icurr+pelem->nnodes()-1; // end global index
        pelem->ifbeg() = fcurr;
        pelem->ifend() = fcurr+nfnodes-1; // end global face index
        icurr += pelem->nnodes();
        fcurr += nfnodes;
      }
      delete pelem2d;

    } // end of element loop for this triangle
  } // end of triangle base mesh loop

} // end of method do_elems3d (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : do_elems2d (2)
// DESC   : Build 2d elemental grid on base mesh
// ARGS   : p      : matrix of size the number of elements X GDIM containing 
//                   the poly expansion order in each direction
//          gxnodes: vector of GDIM vectors containing Cartesian coords of elements
// RETURNS: none.
//**********************************************************************************
void GGridIcos::do_elems2d(GTMatrix<GINT> &p,
                           GTVector<GTVector<GFTYPE>> &gxnodes)
{
  GString                     serr = "GridIcos::do_elems2d (2): ";
  GElem_base                  *pelem;
  GTVector<GTVector<GFTYPE>>  *xNodes;
  GTVector<GTVector<GFTYPE>*> *xiNodes;
  GTVector<GNBasis<GCTYPE,GFTYPE>*>
                               gb(GDIM);
  GTVector<GINT>               ppool(gbasis_.size());

  // Now, treat the gbasis_ as a pool that we search
  // to find bases we need:
  for ( GSIZET j=0; j<ppool.size(); j++ ) ppool[j] = gbasis_[j]->getOrder();


  // Set element internal dof from input data:
  GSIZET iwhere ;
  GSIZET nvnodes;   // no. vol nodes
  GSIZET nfnodes;   // no. face nodes
  GSIZET icurr = 0; // current global index
  GSIZET fcurr = 0; // current global face index
  // For each triangle in base mesh owned by this rank...
  for ( GSIZET i=0; i<p.size(1); i++ ) { 
    nvnodes = 1;
    for ( GSIZET j=0; j<GDIM; j++ ) { // set basis from pool
      assert(ppool.contains(p(i,j),iwhere) && "Expansion order not found");
      gb[j] = gbasis_[iwhere];
      nvnodes *= (p(i,j) + 1);
    }
    pelem = new GElem_base(GE_2DEMBEDDED, gb);
    xNodes  = &pelem->xNodes();  // node spatial data
    xiNodes = &pelem->xiNodes(); // node ref interval data
#if 0
    bdy_ind = &pelem->bdy_indices(); // get bdy indices data member
    bdy_typ = &pelem->bdy_types  (); // get bdy types data member
    bdy_ind->clear(); bdy_typ->clear();
#endif

    // Set internal node positions from input data.
    // Note that gxnodes are 'global' and xNodes is
    // element-local:
    for ( GSIZET j=0; j<GDIM; j++ ) {
       gxnodes[j].range(icurr, icurr+nvnodes-1);
      (*xNodes)[j] = gxnodes[j];
    }
    for ( GSIZET j=0; j<GDIM; j++ ) gxnodes[j].range_reset();

    pelem->init(*xNodes);
    gelems_.push_back(pelem);

    assert(nvnodes == gelems_[i]->nnodes() && "Incompatible node count");
    nfnodes = gelems_[i]->nfnodes();
    pelem->igbeg() = icurr;      // beginning global index
    pelem->igend() = icurr+nvnodes-1; // end global index
    pelem->ifbeg() = fcurr;
    pelem->ifend() = fcurr+nfnodes-1; // end global face index
    icurr += nvnodes;
    fcurr += nfnodes;
  } // end of triangle base mesh loop

} // end of method do_elems2d (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : do_elems3d (2)
// DESC   : Build 3d elemental grid. It's assumed that init3d has been called
//          prior to entry, and that the base icos grid has beed set. 
// ARGS   : p      : matrix of size the number of elements X GDIM containing 
//                   the poly expansion order in each direction
//          gxnodes: vector of GDIM vectors containing Cartesian coords of elements
// RETURNS: none.
//**********************************************************************************
void GGridIcos::do_elems3d(GTMatrix<GINT> &p,
                           GTVector<GTVector<GFTYPE>> &gxnodes)
{
  GString                      serr = "GridIcos::do_elems3d (2): ";
  GElem_base                  *pelem;
  GTVector<GTVector<GFTYPE>>  *xNodes;
  GTVector<GTVector<GFTYPE>*> *xiNodes;
  GTVector<GNBasis<GCTYPE,GFTYPE>*>
                               gb(GDIM);
  GTVector<GINT>               ppool(gbasis_.size());

  // Now, treat the gbasis_ as a pool that we search
  // to find bases we need:
  for ( GSIZET j=0; j<ppool.size(); j++ ) ppool[j] = gbasis_[j]->getOrder();


  // Set element internal dof from input data:
  GSIZET iwhere ;
  GSIZET nvnodes;   // no. vol nodes
  GSIZET nfnodes;   // no. face nodes
  GSIZET icurr = 0; // current global index
  GSIZET fcurr = 0; // current global face index
  // For each triangle in base mesh owned by this rank...
  for ( GSIZET i=0; i<p.size(1); i++ ) { 
    nvnodes = 1;
    for ( GSIZET j=0; j<GDIM; j++ ) { // set basis from pool
      assert(ppool.contains(p(i,j),iwhere) && "Expansion order not found");
      gb[j] = gbasis_[iwhere];
      nvnodes *= (p(i,j) + 1);
    }
    pelem = new GElem_base(GE_DEFORMED, gb);
    xNodes  = &pelem->xNodes();  // node spatial data
    xiNodes = &pelem->xiNodes(); // node ref interval data
#if 0
    bdy_ind = &pelem->bdy_indices(); // get bdy indices data member
    bdy_typ = &pelem->bdy_types  (); // get bdy types data member
    bdy_ind->clear(); bdy_typ->clear();
#endif

    // Set internal node positions from input data.
    // Note that gxnodes are 'global' and xNodes is
    // element-local:
    for ( GSIZET j=0; j<GDIM; j++ ) {
       gxnodes[j].range(icurr, icurr+nvnodes-1);
      (*xNodes)[j] = gxnodes[j];
    }
    for ( GSIZET j=0; j<GDIM; j++ ) gxnodes[j].range_reset();

    pelem->init(*xNodes);
    gelems_.push_back(pelem);

    assert(nvnodes == gelems_[i]->nnodes() && "Incompatible node count");
    nfnodes = gelems_[i]->nfnodes();
    pelem->igbeg() = icurr;      // beginning global index
    pelem->igend() = icurr+nvnodes-1; // end global index
    pelem->ifbeg() = fcurr;
    pelem->ifend() = fcurr+nfnodes-1; // end global face index
    icurr += nvnodes;
    fcurr += nfnodes;
  } // end of triangle base mesh loop

} // end of method do_elems3d (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : print
// DESC   : Print final (global triangular) base mesh. It is often the case 
//          that we will print the mesh to a file for visualization, so this is a 
//          utility that allows us to do this easily. A stream operator is 
//          still provided to print in a completely formatted way.
// ARGS   : filename: filename to print to
//          icoord  : GCOORDSYST: GICOS_CART or GICOS_LATLONG, for cartesian or
//                    lat-long where appropriate
// RETURNS: none.
//**********************************************************************************
void GGridIcos::print(const GString &filename, GCOORDSYST icoord)
{
  GString serr = "GridIcos::print: ";
  std::ofstream ios;

  GTPoint<GFTYPE> pt;
  GFTYPE          r, xlat, xlong;

  ios.open(filename);
  if ( icoord == GICOS_LATLONG) { // print in lat-long
    for ( GSIZET i=0; i<tmesh_.size(); i++ ) { // for each triangle
      for ( GSIZET j=0; j<3; j++ ) { // for each vertex of triangle
        pt = *tmesh_[i].v[j];
        r = sqrt(pt.x1*pt.x1 + pt.x2*pt.x2 + pt.x3*pt.x3);
        xlat  = asin(pt.x3/r);
        xlong = atan2(pt.x2,pt.x1);
        xlong = xlong < 0.0 ? 2*PI+xlong : xlong;
#if defined(_G_IS2D)
        ios << xlat << " " <<  xlong << std::endl;
#elif defined(_G_IS3D)
        ios << r << " " << xlat << " " <<  xlong << std::endl;
#endif
      }
    }
  }
  else if ( icoord == GICOS_CART ) { // print in Cartesian
    for ( GSIZET i=0; i<tmesh_.size(); i++ ) { // for each triangle
      for ( GSIZET j=0; j<3; j++ ) { // for each vertex of triangle
        pt = *tmesh_[i].v[j];
        ios << pt.x1 << " " << pt.x2 << " " << pt.x3 << std::endl ;
      }
    }
  }

  ios.close();

} // end of method print


//**********************************************************************************
//**********************************************************************************
// METHOD : project2sphere (1)
// DESC   : Project Cartesian coords to sphere, specified by rad argument, and
//          express in Cartesian coords.  Necessary for 2d grids.
// ARGS   : tmesh: Vector of triangles (vertices), modified to contain 
//                 projected coordinates
// RETURNS: none.
//**********************************************************************************
void GGridIcos::project2sphere(GTVector<GTriangle<GFTYPE>> &tmesh, GFTYPE rad)
{
  GString serr = "GridIcos::project2sphere (1): ";

  GFTYPE r, xlat, xlong;
  GTPoint<GFTYPE> v;

  for ( GSIZET i=0; i<tmesh_.size(); i++ ) { // loop over all triangles in tmesh
    for ( GSIZET j=0; j<3; j++ ) { // guaranteed to have 3 points each
      v = *tmesh[i].v[j];
      r = v.norm();
      xlat  = asin(v.x3/r);
      xlong = atan2(v.x2,v.x1);
      xlong = xlong < 0.0 ? 2*PI+xlong : xlong;
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
void GGridIcos::project2sphere(GTVector<GTPoint<GFTYPE>> &plist, GFTYPE rad)
{
  GString serr = "GridIcos::project2sphere (2): ";

  GFTYPE r, xlat, xlong;
  GTPoint<GFTYPE> v;

  for ( GSIZET i=0; i<plist.size(); i++ ) { // loop over all points
    v = plist[i];
    r = v.norm();
    xlat  = asin(v.x3/r);
    xlong = atan2(v.x2,v.x1);
    xlong = xlong < 0.0 ? 2*PI+xlong : xlong;
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
void GGridIcos::project2sphere(GTVector<GTVector<GFTYPE>> &plist, GFTYPE rad)
{
  GString serr = "GridIcos::project2sphere (3): ";

  GFTYPE r, xlat, xlong, x, y, z;

  for ( GSIZET i=0; i<plist[0].size(); i++ ) { // loop over all points
    x = plist[0][i]; y = plist[1][i]; z = plist[2][i];
    r = sqrt(x*x + y*y + z*z);
    xlat  = asin(z/r);
    xlong = atan2(y,x);
    xlong = xlong < 0.0 ? 2*PI+xlong : xlong;
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
void GGridIcos::spherical2xyz(GTVector<GTPoint<GFTYPE>*> &plist)
{
  GString serr = "GridIcos::spherical2xyz(1): ";

  GFTYPE r, xlat, xlong;

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
void GGridIcos::spherical2xyz(GTVector<GTPoint<GFTYPE>> &plist)
{
  GString serr = "GridIcos::spherical2xyz(2): ";

  GFTYPE r, xlat, xlong;

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
void GGridIcos::spherical2xyz(GTVector<GTVector<GFTYPE>> &plist)
{
  GString serr = "GridIcos::spherical2xyz(3): ";
  assert(plist.size() >= 3 && "Invalid dimensionality on input array");

  GFTYPE r, xlat, xlong;

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
void GGridIcos::xyz2spherical(GTVector<GTPoint<GFTYPE>*> &plist)
{
  GString serr = "GridIcos::xyz2spherica(1): ";

  GFTYPE r, x, y, z;

  for ( GSIZET i=0; i<plist.size(); i++ ) { // loop over all points
   x = plist[i]->x1; y = plist[i]->x2; z = plist[i]->x3;
   r = sqrt(x*x + y*y + z*z);
   plist[i]->x1 = r;
   plist[i]->x2 = asin(z/r);
   plist[i]->x3 = atan2(y,x);
   plist[i]->x3 = plist[i]->x3 < 0.0 ? 2*PI+plist[i]->x3 : plist[i]->x3;
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
void GGridIcos::xyz2spherical(GTVector<GTVector<GFTYPE>> &plist)
{
  GString serr = "GridIcos::xyz2spherica(2): ";

  GFTYPE r, x, y, z;

  for ( GSIZET i=0; i<plist.size(); i++ ) { // loop over all points
   x = plist[i][0]; y = plist[i][1]; z = plist[i][2];
   r = sqrt(x*x + y*y + z*z);
   plist[i][0] = r;
   plist[i][1] = asin(z/r);
   plist[i][2] = atan2(y,x);
   plist[i][2] = plist[i][2] < 0.0 ? 2*PI+plist[i][2] : plist[i][2];
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
void GGridIcos::xyz2spherical(GTVector<GTPoint<GFTYPE>> &plist)
{
  GString serr = "GridIcos::xyz2spherica(3): ";

  GFTYPE r, x, y, z;

  for ( GSIZET i=0; i<plist.size(); i++ ) { // loop over all points
   x = plist[i].x1; y = plist[i].x2; z = plist[i].x3;
   r = sqrt(x*x + y*y + z*z);
   plist[i].x1 = r;
   plist[i].x2 = asin(z/r);
   plist[i].x3 = atan2(y,x);
   plist[i].x3 = plist[i].x3 < 0.0 ? 2*PI+plist[i].x3 : plist[i].x3;
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
void GGridIcos::cart2gnomonic(GTVector<GTPoint<GFTYPE>> &clist, GFTYPE rad, GFTYPE xlatc, GFTYPE xlongc, GTVector<GTPoint<GFTYPE>> &glist)
{
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

  GFTYPE xlatp, xlongp;
  GFTYPE den, r, xlat, xlong;
  for ( GSIZET i=0; i<clist.size(); i++ ) { // loop over all points
    r      = clist[i].norm();
    xlat   = asin(clist[i].x3/r);
    xlong  = atan2(clist[i].x2,clist[i].x1);
    xlong  = xlong < 0.0 ? 2*PI+xlong : xlong;
#if 0
    den    = sin(xlatc)*sin(xlat) + cos(xlatc)*cos(xlat)*cos(xlong-xlongc);  
    xlatp  = asin( cos(xlatc)*sin(xlat) - sin(xlatc)*cos(xlat)*cos(xlong-xlongc) );
    xlongp = atan2( cos(xlat)*sin(xlong-xlongc),den );

    glist[i].x1 = rad*tan(xlongp);
    glist[i].x2 = rad*tan(xlatp)*sec(xlongp);
#else
    den    = sin(xlatc)*sin(xlat) + cos(xlatc)*cos(xlat)*cos(xlong-xlongc); 
    den    = fabs(den) < std::numeric_limits<GFTYPE>::epsilon() ? 0.0 : 1.0/den;
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
void GGridIcos::gnomonic2cart(GTVector<GTVector<GFTYPE>> &glist, GFTYPE rad, GFTYPE xlatc, GFTYPE xlongc, GTVector<GTVector<GFTYPE>> &clist)
{
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

  GFTYPE beta, rho, x, xlat, xlong, y;
  GFTYPE eps = std::numeric_limits<GFTYPE>::epsilon();
  for ( GSIZET i=0; i<glist[0].size(); i++ ) { // loop over all points
    x      = glist[0][i];
    y      = glist[1][i];


    if ( fabs(x) < eps && fabs(y) < eps ) {
      xlat = xlatc;
      xlong = xlongc;
    } else {
      rho    = sqrt(x*x + y*y);
      beta   = atan(rho);
      xlat   = asin( cos(beta)*sin(xlatc) + (y*sin(beta)*cos(xlatc))/rho );
      xlong  = xlongc + atan2( x*sin(beta), 
                               rho*cos(xlatc)*cos(beta)-y*sin(xlatc)*sin(beta) );
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
void GGridIcos::gnomonic2cart(GTVector<GTPoint<GFTYPE>> &glist, GFTYPE rad, GFTYPE xlatc, GFTYPE xlongc, GTVector<GTPoint<GFTYPE>> &clist)
{
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

  GFTYPE beta, rho, x, xlat, xlong, y;
  GFTYPE eps = std::numeric_limits<GFTYPE>::epsilon();
  for ( GSIZET i=0; i<clist.size(); i++ ) { // loop over all points
    x      = glist[i][0];
    y      = glist[i][1];

    if ( fabs(x) < eps && fabs(y) < eps ) {
      xlat = xlatc;
      xlong = xlongc;
    } else {
      rho    = sqrt(x*x + y*y);
      beta   = atan(rho);
      xlat   = asin( cos(beta)*sin(xlatc) + (y*sin(beta)*cos(xlatc))/rho );
      xlong  = xlongc + atan2( x*sin(beta), 
                               rho*cos(xlatc)*cos(beta)-y*sin(xlatc)*sin(beta) );
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
// DESC   : Reorder specified vertices to be consistent with
//          shape functions
// ARGS   : uverts : list of unordered vertices
//          overts : array of ordered vertices, returned
// RETURNS: none.
//**********************************************************************************
void GGridIcos::reorderverts2d(GTVector<GTPoint<GFTYPE>> &uverts, GTVector<GSIZET> &isort, 
                               GTVector<GTPoint<GFTYPE>> &overts)
{
  GString serr = "GridIcos::reorderverts2d: ";

  assert(uverts.size() == 4 && "Incorrect number of vertices");


  GTVector<GFTYPE> x(4);
  GTVector<GSIZET> Ixy(4);


  isort.resize(4);
  for ( GSIZET i=0; i<uverts.size(); i++ ) { 
    x[i] = uverts[i].x1;
  }

  x.sortincreasing(Ixy);

  // Do 'left' -hand vertices:
  if ( uverts[Ixy[0]].x2 < uverts[Ixy[1]].x2 ) {
    isort[0] = Ixy[0];
    isort[3] = Ixy[1];
  } else {
    isort[0] = Ixy[1];
    isort[3] = Ixy[0];
  }

  // Do 'right' -hand vertices:
  if ( uverts[Ixy[2]].x2 < uverts[Ixy[3]].x2 ) {
    isort[1] = Ixy[2];
    isort[2] = Ixy[3];
  } else {
    isort[1] = Ixy[3];
    isort[2] = Ixy[2];
  }

  for ( GSIZET j=0; j<4; j++ ) overts[j] = uverts[isort[j]];
  
} // end of method reorderverts2d


//**********************************************************************************
//**********************************************************************************
// METHOD : order_latlong2d
// DESC   : Order 2d vertices on exit s.t. they roughly define a 'box' 
//          in spherical coords. Used only for quad elements.
// ARGS   : verts : Array of vertices, re-ordered on exit
// RETURNS: none.
//**********************************************************************************
void GGridIcos::order_latlong2d(GTVector<GFPoint> &verts)
{

  assert(verts.size() == 4 && "4 vertices must be provided");

  GString           serr = "GGridIcos::order_latlong2d: ";
  GTVector<GSIZET>  isortlon(4);
  GTVector<GFTYPE>  lon(4), lat(4);
  GTVector<GFPoint> cverts(4);       // copy of input verts
  GTVector<GFPoint> sverts(4);       // sverts in sph coords

  cverts = verts;
  sverts = verts;
  xyz2spherical(sverts); // convert verts to latlon

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
void GGridIcos::order_triangles(GTVector<GTriangle<GFTYPE>> &tmesh)
{

  GString           serr = "GGridIcos::order_triangles: ";
  GTVector<GSIZET>  isortlon(3);
  GTVector<GFTYPE>  lon(3), lat(3);
  GTVector<GFPoint> cverts(3);       // copy of input verts
  GTVector<GFPoint> sverts(3);       // sverts in sph coords

  iup_.resize(tmesh.size());
  iup_ = 0;
  for ( GSIZET i=0; i<tmesh.size(); i++ ) {
    for ( GSIZET j=0; j<3; j++ ) cverts[j] = *tmesh[i].v[j];
    sverts = cverts;
    xyz2spherical(sverts); // convert verts to latlon
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


//**********************************************************************************
//**********************************************************************************
// METHOD : config_bdy
// DESC   : Configure 3d spherical boundaries from ptree
// ARGS   : 
//          ptree : main prop tree 
//          igbdy : For each natural/canonical global boundary face,
//                  gives vector of global bdy ids
//          igbdyt: bdy type ids for each index in igbdy
// RETURNS: none.
//**********************************************************************************
void GGridIcos::config_bdy(const PropertyTree &ptree, 
                           GTVector<GTVector<GSIZET>> &igbdy, 
                           GTVector<GTVector<GBdyType>> &igbdyt)
{
  // Cycle over all geometric boundaries, and configure:

  GBOOL              bret, buniform=FALSE;
  GSIZET             iwhere;
  GTVector<GBOOL>    uniform(2);
  GTVector<GBdyType> bdytype(2);
  GTVector<GBdyType> btmp;
  GTVector<GSIZET>   itmp;
  GTVector<GFTYPE>   rbdy(2);
  GTVector<GString>  bdynames(2);
  GTVector<GString>  confmthd (2);
  GTVector<GString>  bdyupdate(2);
  GString            gname, sbdy, bdyclass;
  PropertyTree       bdytree, gridptree, spectree;

  // Clear input arrays:
  igbdy .clear();
  igbdyt.clear();

  if ( ndim_ == 2 ) return; // no boundaries to configure
 
  bdynames[0] = "bdy_inner";
  bdynames[1] = "bdy_outer";

  gname     = ptree.getValue<GString>("grid_type");
  gridptree = ptree.getPropertyTree(gname);


  bdyupdate = gridptree.getValue<GString>("update_method","");
  rbdy[0] = radiusi_;
  rbdy[1] = radiuso_;

  // Get properties from the main prop tree. 
  // Note: bdys are configured by way of geometry's
  //       natural decomposition: here, by inner and
  //       outer spherical surfaces. But the bdy indices 
  //       and types returned on exist contain info for all bdys:
  for ( auto j=0; j<2; j++ ) { // cycle over 2 spherical surfaces
    sbdy         = gridptree.getValue<GString>(bdynames[j]);
    bdytree      = gridptree.getPropertyTree(sbdy);
    bdyclass     = bdytree.getValue<GString>("bdy_class", "uniform");
    bdytype  [j] = geoflow::str2bdytype(bdytree.getValue<GString>("base_type", "GBDY_NONE"));
    assert(bdytype  [j] == GBDY_PERIODIC && "Invalid boundary condition");
    uniform  [j] = bdyclass == "uniform" ? TRUE : FALSE;
    confmthd [j] = bdytree.getValue<GString>("bdy_config_method","");
    buniform     = buniform || uniform[j];
  }

  // Handle non-uniform (user-configured) bdy types first;
  // Note: If "uniform" not specified for a boundary, then
  //       user MUST supply a method to configure it.
  //       Also, each natural face may be configured independently,
  //       but the bdy indices & corresp. types are concatenated into 
  //       single arrays:
  for ( auto j=0; j<2; j++ ) { 
    // First, find global bdy indices:
    if ( uniform[j] ) continue;
    find_bdy_ind3d(rbdy[j], itmp); 
    spectree  = ptree.getPropertyTree(confmthd[j]);
    bret = GSpecBdyFactory::dospec(spectree, *this, j, itmp, btmp); // get user-defined bdy spec
    assert(bret && "Boundary specification failed");
    igbdy [j].resize(itmp.size()); igbdy [j] = itmp;
    igbdyt[j].resize(itmp.size()); igbdyt[j] = btmp;

  }
  
  // Fill in uniform bdy types:
  for ( auto j=0; j<2; j++ ) { 
    if ( !uniform[j] ) continue;
    // First, find global bdy indices:
    find_bdy_ind3d(rbdy[j], itmp); 
    // Set type for each bdy index:
    for ( auto i=0; i<itmp.size(); i++ ) {
      btmp[i] = bdytype[j]; 
    }
    igbdy [j].resize(itmp.size()); igbdy [j] = itmp;
    igbdyt[j].resize(itmp.size()); igbdyt[j] = btmp;
  }


} // end of method config_bdy


//**********************************************************************************
//**********************************************************************************
// METHOD : find_bdy_ind3d
// DESC   : Find global bdy indices (indices into xNodes_ arrays) that
//          corresp to specified radius
// ARGS   : radius   : radius
//          ibdy     : array of indices into xNodes that comprise this boundary
// RETURNS: none.
//**********************************************************************************
void GGridIcos::find_bdy_ind3d(GFTYPE radius, GTVector<GSIZET> &ibdy)
{

  GFTYPE          eps, r;
  GTPoint<GFTYPE> pt(ndim_);

  ibdy.clear();
  eps = 100*std::numeric_limits<GFTYPE>::epsilon();

  for ( GSIZET i=0; i<xNodes_[0].size(); i++ ) { // face 0
      r = sqrt(pow(xNodes_[0][i],2)+pow(xNodes_[1][i],2)+pow(xNodes_[2][i],2));
      if ( FUZZYEQ(r, radius, eps) ) {
        ibdy.push_back(i);
      }
  }

} // end, method find_bdy_ind3d


//**********************************************************************************
//**********************************************************************************
// METHOD : do_face_normals
// DESC   : Compute normals to each element face
// ARGS   : none 
// RETURNS: none
//**********************************************************************************
void GGridIcos::do_face_normals()
{

  #if defined(_G_IS2D)
    do_face_normals2d();
  #elif defined(_G_IS3D)
    do_face_normals3d();
  #else
    #error Invalid problem dimensionality
  #endif

} // end, method do_face_normals


//**********************************************************************************
//**********************************************************************************
// METHOD : do_face_normals2d
// DESC   : Compute normals to each element face in 2d
// ARGS   : none 
// RETURNS: none
//**********************************************************************************
void GGridIcos::do_face_normals2d()
{

  // Cycle through local elem face indices to set
  // normals. Taken in order, these should correspond


} // end, method do_bdy_normals2d


//**********************************************************************************
//**********************************************************************************
// METHOD : do_face_normals3d
// DESC   : Compute normals to each element face in 3d
// ARGS   : none 
// RETURNS: none
//**********************************************************************************
void GGridIcos::do_face_normals3d()
{


} // end, method do_bdy_normals3d



//**********************************************************************************
//**********************************************************************************
// METHOD : do_bdy_normals
// DESC   : Compute normals to each domain bdy
// ARGS   : none 
// RETURNS: none
//**********************************************************************************
void GGridIcos::do_bdy_normals()
{

  #if defined(_G_IS2D)
    return;
  #elif defined(_G_IS3D)
    do_bdy_normals3d();
  #else
    #error Invalid problem dimensionality
  #endif

} // end, method do_bdy_normals


//**********************************************************************************
//**********************************************************************************
// METHOD : do_bdy_normals3d
// DESC   : Compute normals to each domain bdy in 3d
// ARGS   : none 
// RETURNS: none
//**********************************************************************************
void GGridIcos::do_bdy_normals3d()
{

} // end, method do_bdy_normals3d
