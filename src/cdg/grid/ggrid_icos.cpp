//==================================================================================
// Module       : ggrid_icos
// Date         : 8/31/18 (DLR)
// Description  : Object defining a (global) icosahedral grid, that
//                uses gnomonic projections to locate element vertices.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <omp.h>
#include "gelem_base.hpp"
#include "ggrid_icos.hpp"



//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Instantiate for 2d grid
// ARGS   : traits:
//            ilevel  : GINT refinement level. Level 0 is just the icos generator 
//            radiusi : GFTYPE radius
//            radiuso : GFTYPE outer radius (not used here)
//            bdyTypes: GTVector<GFTYPE>(2) (not used here)
//          b     : vector of basis pointers, of size at least ndim=2
//          nprocs : no. MPI tasks
// RETURNS: none
//**********************************************************************************
GGridIcos::GGridIcos(GGridIcos::Traits &traits, GTVector<GNBasis<GCTYPE,GFTYPE>*> &b, GINT nprocs) :
ilevel_      (traits.ilevel),
ndim_                    (2),
radiusi_    (traits.radiusi),
radiuso_                 (0), // don't need this one
nprocs_             (nprocs),
gdd_               (NULLPTR),
lshapefcn_         (NULLPTR),
bdycallback_       (NULLPTR)
{

  gbasis_.resize(b.size());
  gbasis_ = b;
  lshapefcn_ = new GShapeFcn_linear();

  init2d();
} // end of constructor method (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (2)
// DESC   : Instantiate for 3d grid
// ARGS   : traits:
//            radiusi : GFTYPE inner radius
//            radiuso : GFTYPE outer radius
//            bdyTypes: GTVector<GFTYPE>(2)
//          ne     : no. elements in each coord direction (r,theta,phi)=(r,lat,long)
//          b      : vector of basis pointers, of size at least ndim=2. 
//          nprocs : no. MPI tasks
// RETURNS: none
//**********************************************************************************
GGridIcos::GGridIcos(GGridIcos::Traits &traits, GTVector<GINT> &ne, GTVector<GNBasis<GCTYPE,GFTYPE>*> &b, GINT nprocs) :
ilevel_           (0),
ndim_             (3),
radiusi_    (traits.radiusi),
radiuso_    (traits.radiuso),
nprocs_      (nprocs),
gdd_        (NULLPTR),
lshapefcn_  (NULLPTR),
bdycallback_(NULLPTR)
{
  assert(b.size() == 3 && "Insufficient number of bases");
  gbasis_.resize(b.size());
  gbasis_ = b;
  lshapefcn_ = new GShapeFcn_linear();
  ne_     = ne;

  global_bdy_types_ = traits.bdyTypes;
  assert(!global_bdy_types_.contains(GBDY_PERIODIC) && "Invalid boundary type");

  init3d();
} // end of constructor method (2)


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
// METHOD : set_basis 
// DESC   : Set basis object for coord. direction, idir in (1, 2, 3)
// ARGS   :
// RETURNS: none
//**********************************************************************************
void GGridIcos::set_basis(GTVector<GNBasis<GCTYPE,GFTYPE>*> &b)
{
  GString serr = "GridIcos::set_basis: ";

  gbasis_.resize(b.size());
  gbasis_ = b;
} // end of method set_basis


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

  
  // Normalize vertices to be at radiusi_:
  GFTYPE fact = radiusi_/sqrt(phi*phi + 1.0);
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

  GINT nelems = ne_[0]*ne_[1]*ne_[2];
  hmesh_.resize(nelems);

  GFTYPE dr    = (radiuso_ - radiusi_) / static_cast<GFTYPE>(ne_[0]);  // delta r
  GFTYPE dlat  = (PI/2.0) / static_cast<GFTYPE>(ne_[1]); // delta lat 
  GFTYPE dlong = (2.0*PI) / static_cast<GFTYPE>(ne_[2]); // delta  long

  // Compute vertices of all hexes ('cubes'), starting 
  // at (r, lat, long) = (radiusi_, 0, 0):
  GINT   n=0;
  GTPoint<GFTYPE> v0(3);     // starting sph. point of elem
  GTPoint<GFTYPE> dv(3);     // delta-point
  for ( GSIZET k=0; k<ne_[2]; k++ ) { 
    for ( GSIZET j=0; j<ne_[1]; j++ ) { 
      for ( GSIZET i=0; i<ne_[0]; i++ ) { 
        v0.x1 = i*dr + radiusi_; v0.x2 = j*dlat - 0.5*PI; v0.x3 = k*dlong - 0.5*PI; 
         
        dv.x1 = 0.0; dv.x2 = dlat; dv.x3 = 0.0  ;  hmesh_[n].v1 = v0 + dv;
        dv.x1 = 0.0; dv.x2 = dlat; dv.x3 = dlong;  hmesh_[n].v2 = v0 + dv;
        dv.x1 = 0.0; dv.x2 = 0.0 ; dv.x3 = dlong;  hmesh_[n].v3 = v0 + dv;
        dv.x1 = dr ; dv.x2 = 0.0 ; dv.x3 = 0.0  ;  hmesh_[n].v4 = v0 + dv;
        dv.x1 = dr ; dv.x2 = dlat; dv.x3 = 0.0  ;  hmesh_[n].v5 = v0 + dv;
        dv.x1 = dr ; dv.x2 = dlat; dv.x3 = dlong;  hmesh_[n].v6 = v0 + dv;
        dv.x1 = dr ; dv.x2 = 0.0 ; dv.x3 = dlong;  hmesh_[n].v7 = v0 + dv;

        spherical2xyz(hmesh_[n].v);  
        n++;
      }
    }
  }
  
  // Compute centroids of all hexes ('cubes'):
  GTPoint<GFTYPE> a(3);

  ftcentroids_.clear();
  for ( GSIZET j=0; j<hmesh_.size(); j++ ) { // for each cube
    for ( GSIZET k=0; k<hmesh_[j].v.size(); k++ ) a += *hmesh_[j].v[k];
    a *= (1.0/hmesh_[j].v.size());
    ftcentroids_.push_back(a);
  }

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

  
  // Project all vertices to sphere:
  project2sphere(tmesh_,radiusi_);

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
// METHOD : do_grid
// DESC   : Public entry point for grid computation
// ARGS   : grid: GGrid object
//          rank: MPI rank or partition id
// RETURNS: none.
//**********************************************************************************
void GGridIcos::do_grid(GGrid &grid, GINT irank)
{
  if ( ndim_ == 2 ) do_grid2d(grid, irank);
  if ( ndim_ == 3 ) do_grid3d(grid, irank);
  grid.do_typing();

  // Inititialized global grid quantities:
  grid.init();

} // end, method do_grid


//**********************************************************************************
//**********************************************************************************
// METHOD : do_grid2d
// DESC   : Build 2d elemental grid on base mesh
// ARGS   : grid: GGrid object
//          rank: MPI rank or partition id
// RETURNS: none.
//**********************************************************************************
void GGridIcos::do_grid2d(GGrid &grid, GINT irank)
{
  GString           serr = "GridIcos::do_grid2d: ";
  GFTYPE            fact, xlatc, xlongc;
  GTVector<GFPoint> cverts(4), gverts(4), tverts(4);
  GTPoint<GFTYPE>   c(3), ct(3), v1(3), v2(3), v3(3); // 3d points
  GElem_base        *pelem;
  GElemList         *gelems = &grid.elems();
  GTVector<GINT>    iind;
  GTVector<GINT>    I(1);
  GTVector<GFTYPE>  Ni;
  GTVector<GTVector<GFTYPE>>  *xNodes;
  GTVector<GTVector<GFTYPE>*> *xiNodes;
  GTVector<GTVector<GFTYPE>>   xgtmp(3);


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
      switch (j) {
      case 0: 
      cverts[0] = v1; cverts[1] = (v1+v2)*0.5; cverts[2] = ct; cverts[3] = (v1+v3)*0.5; break;
      case 1: 
      cverts[0] = (v1+v2)*0.5; cverts[1] = v2; cverts[2] = (v2+v3)*0.5; cverts[3] = ct; break;
      case 2: 
      cverts[0] = ct; cverts[1] = (v2+v3)*0.5; cverts[2] = v3; cverts[3] = (v1+v3)*0.5; break;
      }
      
      xNodes  = &pelem->xNodes();  // node spatial data
      xiNodes = &pelem->xiNodes(); // node ref interval data
      Ni.resize(pelem->nnodes());

      project2sphere(cverts, radiusi_); // project verts to sphere     
      c = (cverts[0] + cverts[1] + cverts[2] + cverts[3])*0.25; // elem centroid
      xlatc  = asin(c.x3/radiusi_); // reference lat/long
      xlongc = atan2(c.x2,c.x1);

      cart2gnomonic(cverts, radiusi_, xlatc, xlongc, tverts); // gnomonic vertices of quads
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
      gnomonic2cart(xgtmp, radiusi_, xlatc, xlongc, *xNodes); //
      project2sphere(*xNodes, radiusi_);
      pelem->init(*xNodes);
      gelems->push_back(pelem);

      nfnodes = 0;
      for ( GSIZET j=0; j<(*gelems)[n]->nfaces(); j++ )  // get # face nodes
        nfnodes += (*gelems)[n]->face_indices(j).size();
      pelem->igbeg() = icurr;      // beginning global index
      pelem->igend() = icurr+pelem->nnodes()-1; // end global index
      pelem->ifbeg() = fcurr;
      pelem->ifend() = fcurr+nfnodes-1; // end global face index
      icurr += pelem->nnodes();
      fcurr += nfnodes;

    } // end of element loop for this triangle
  } // end of triangle base mesh loop

  // Can set individual nodes and internal bdy conditions
  // with callback here:
  if ( bdycallback_ != NULLPTR ) {
    (*bdycallback_)(*gelems);
  }


} // end of method do_grid2d


//**********************************************************************************
//**********************************************************************************
// METHOD : do_grid3d
// DESC   : Build 3d elemental grid. It's assumed that init3d has been called
//          prior to entry, and that the hmesh_ has beed set. This arrays
//          of hexagonal 3d elements provides for each hex its vertices in Cartesian
//          coordinates. Also, the centroids of these hexes (in Cartesian coords)
//          should also have been computed. The Cartesian hex vertices will be
//          converted into sherical coordinates, and the GShapeFcn_linear
//          will be used in spherical coords to compute the interior nodes. These
//          will then be converted to Cartesian coords.
// ARGS   : grid: GGrid object
//          rank: MPI rank or partition id
// RETURNS: none.
//**********************************************************************************
void GGridIcos::do_grid3d(GGrid &grid, GINT irank)
{
  GString                      serr = "GridIcos::do_grid3d: ";
  GTVector<GINT>               iind;
  GTVector<GINT>               I(3);
  GTVector<GINT>              *bdy_ind;
  GTVector<GINT>              *face_ind;
  GTVector<GFTYPE>             Ni;
  GElemList                   *gelems = &grid.elems();
  GElem_base                  *pelem;
  GTVector<GTVector<GFTYPE>>  *xNodes;
  GTVector<GTVector<GFTYPE>*> *xiNodes;
  GTVector<GTVector<GFTYPE>>   xgtmp(3);


  assert(gbasis_.size()>0 && "Must set basis first");

  if ( gdd_       == NULLPTR ) gdd_ = new GDD_base(nprocs_);

  // Get indirection indices to create elements
  // only for this task:
  gdd_->doDD(ftcentroids_, irank, iind);

  GSIZET i, n;
  GSIZET nfnodes;   // no. face nodes
  GSIZET icurr = 0; // current global index
  GSIZET fcurr = 0; // current global face index
  for ( GSIZET n=0; n<iind.size(); n++ ) { // for each hex in irank's mesh
    i = iind[n];

    // Transform hex vertices to (r,lat,long):
    xyz2spherical(hmesh_[i].v);

    pelem = new GElem_base(GE_DEFORMED, gbasis_);
    xNodes  = &pelem->xNodes();  // node spatial data
    xiNodes = &pelem->xiNodes(); // node ref interval data
    Ni.resize(pelem->nnodes()); // tensor product shape function
    bdy_ind = &pelem->bdy_indices(); // get bdy indices data member
    bdy_ind->clear();
    for ( GSIZET l=0; l<3; l++ ) { // loop over sph coords 
      (*xNodes)[l] = 0.0;
      for ( GSIZET m=0; m<8; m++ ) { // loop over verts given in sph coords
        I[0] = m;
        lshapefcn_->Ni(I, *xiNodes, Ni);
        (*xNodes)[l] += Ni * (*(hmesh_[i].v[m]))[l];
      }
    }

    // Set element's bdy indices:
    // ... need to check only vertex 0 for this:
    if ( hmesh_[i].v[0]->x1 == radiusi_ || hmesh_[i].v[0]->x1 == radiuso_ ) { 
      face_ind = &pelem->face_indices(0);
      for ( GSIZET k=0; k<face_ind->size(); k++ ) bdy_ind->push_back((*face_ind)[k]); 
    }

    spherical2xyz(*xNodes); // convert nodal coords to Cartesian coords
    pelem->init(*xNodes);
    gelems->push_back(pelem);

    // Set global bdy types at each bdy_ind (this is a coarse
    // application; finer control may be exercised in callback):
    set_global_bdy_3d(*pelem);

    nfnodes = 0;
    for ( GSIZET j=0; j<(*gelems)[n]->nfaces(); j++ )  // get # face nodes
      nfnodes += (*gelems)[n]->face_indices(j).size();
    pelem->igbeg() = icurr;      // beginning global index
    pelem->igend() = icurr + pelem->nnodes()-1; // end global index
    pelem->ifbeg() = fcurr;
    pelem->ifend() = fcurr+nfnodes-1; // end global face index
    icurr += pelem->nnodes();
    fcurr += nfnodes;
  } // end of hex mesh loop

  // If we have a callback function, set the boundary conditions here:
  if ( bdycallback_ != NULLPTR ) {
    (*bdycallback_)(*gelems);
  }

} // end of method do_grid3d


//**********************************************************************************
//**********************************************************************************
// METHOD : set_bdy_callback
// DESC   : Set the callback object and method pointers for setting bdy conditions
// ARGS   : ptr2obj : pointer to object that hosts method callback
//          callback: method name
// RETURNS: none.
//**********************************************************************************
void GGridIcos::set_bdy_callback(std::function<void(GElemList &)> *callback)
{
  bdycallback_  = callback;
} // end of method set_bdy_callback



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
void GGridIcos::xyz2spherical(GTVector<GTPoint<GFTYPE>> &plist)
{
  GString serr = "GridIcos::xyz2spherica(2): ";

  GFTYPE r, x, y, z;

  for ( GSIZET i=0; i<plist.size(); i++ ) { // loop over all points
   x = plist[i].x1; y = plist[i].x2; z = plist[i].x3;
   r = sqrt(x*x + y*y + z*z);
   plist[i].x1 = r;
   plist[i].x2 = asin(z/r);
   plist[i].x3 = atan2(y,x);
  }

} // end of method xyz2spherical (2)


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
// METHOD : set_global_bdy_3d
// DESC   : Set global boundary conditions in 3d
//            NOTE: element node coordinate may change in this call!
// ARGS   : pelem : Element under consideration
// RETURNS: none.
//**********************************************************************************
void GGridIcos::set_global_bdy_3d(GElem_base &pelem)
{
  GTVector<GINT>              *bdy_ind=&pelem.bdy_indices();
  GTVector<GBdyType>          *bdy_typ=&pelem.bdy_types  ();
  GTVector<GTVector<GFTYPE>>  *xNodes =&pelem.xNodes();

  GSIZET ib;
  GFTYPE x, y, z, r;
  for ( GSIZET m=0; m<2; m++ ) { // for each radiusi face
    for ( GSIZET k=0; k<face_ind->size(); k++ ) { // for each bdy index
      ib = (*bdy_ind)[k];
      x = (*xNodes)[0][ib]; y = (*xNodes)[1][ib]; z = (*xNodes)[2][ib];
      r = sqrt( x*x + y*y + z*z );
      if ( r == radiusi_  )
        bdy_typ->push_back( global_bdy_types_[0] );
    }
  }

  for ( GSIZET m=0; m<2; m++ ) { // for each radiuso face
    for ( GSIZET k=0; k<bdy_ind->size(); k++ ) { // for each bdy index
      ib = (*bdy_ind)[k];
      x = (*xNodes)[0][ib]; y = (*xNodes)[1][ib]; z = (*xNodes)[2][ib];
      r = sqrt( x*x + y*y + z*z );
      if ( r == radius0_  )
        bdy_typ->push_back( global_bdy_types_[1] );
    }
  }

} // end, method set_global_bdy_3d
