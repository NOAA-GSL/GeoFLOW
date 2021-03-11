//==================================================================================
// Module       : ggrid_icos.cpp
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


#include <cstdlib>
#include <memory>
#include <cmath>
//#include "omp.h"  // Are we calling API functions ?
#include "gspecbdy_factory.hpp"
#include "ginitstate_factory.hpp"
#include "gupdatebdy_factory.hpp"
#include "ggrid_icos.hpp"
#include "gtpoint.hpp"
#include "gutils.hpp"
#include "tbox/tracer.hpp"

using namespace geoflow::tbox;
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
:                 GGrid(ptree, b, comm),
ilevel_                             (0),
nrows_                              (0),
ndim_                            (GDIM),
radiusi_                          (0.0),
radiuso_                          (0.0), // don't need this one in 2d
gdd_                          (NULLPTR),
lshapefcn_                    (NULLPTR)
{
  GEOFLOW_TRACE();
  assert(b.size() == GDIM && "Basis has incorrect dimensionality");
  
  GString snorm;
  GString gname   = ptree.getValue<GString>("grid_type");
  assert(gname == "grid_icos" || gname == "grid_sphere");
  geoflow::tbox::PropertyTree gridptree = ptree.getPropertyTree(gname);

  gbasis_.resize(b.size());
  gbasis_ = b;
  lshapefcn_ = new GShapeFcn_linear<GTICOS>(2);
  ilevel_  = gridptree.getValue<GINT>("ilevel");
  sreftype_= gridptree.getValue<GString>("refine_type","GICOS_LAGRANGIAN");
  this->cgtraits_.maxit = gridptree.getValue<GDOUBLE>("maxit");
  this->cgtraits_.tol   = gridptree.getValue<GDOUBLE>("tol");
  snorm                 = gridptree.getValue<GString>("norm_type");
  this->cgtraits_.normtype = LinSolverBase<CGTypePack>::str2normtype(snorm);

  
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
	GEOFLOW_TRACE();
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
  GEOFLOW_TRACE();
  str << " radiusi: " << e.radiusi_;
  str << " radiusi: " << e.radiuso_;
  str << " level  : " << e.ilevel_;
  str << " nrows  : " << e.nrows_;
  str << std::endl << " Centroids: " ;
  for ( auto i=0; i<e.ftcentroids_.size(); i++ )
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
void GGridIcos::set_partitioner(GDD_base<GTICOS> *gdd)
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
	GEOFLOW_TRACE();
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
  GTPoint<GTICOS> pt;
  tbase_.resize(ifv0_.size(1));
  for ( auto j=0; j<tbase_.size(); j++ ) tbase_[j].resize(3);
  for ( auto i=0; i<ifv0_.size(1); i++ ) { // for all triangles:
    for ( auto j=0; j<ifv0_.size(2); j++ ) { // for each vertex:
      for ( auto k=0; k<3; k++ ) pt[k] = fv0_(ifv0_(i,j),k);
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
	GEOFLOW_TRACE();
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
	GEOFLOW_TRACE();
  GString serr = "GridIcos::lagrefine: ";
   
  GLLONG ibeg, j, l, m, n, t;

  if      ( "GICOS_BISECTION" == sreftype_ ) {
    // interpret nrows_ as bisection count:
    nrows_ = pow(2,ilevel_)-1; 
  }
  else if ( "GICOS_LAGRANGIAN" == sreftype_ ) {
    // interpret nrows_ as # 'Lagrangian' subdivisions:
    nrows_ = ilevel_;
  }
  else {
    assert(FALSE && "Invalid subdivision type (GICOS_LAGRANGIAN or GICOS_BISECTION");
  }

  // Re-dimension mesh points to be 3d: Expect
  // 20 * (nrows+1)^2 triangles, and 20 * (nrows+1)^2 * 3 quads
  // Total number of points (in horizontal) is
  // 6(
  tmesh_.resize(20*(nrows_*(nrows_+2)+1)); // refined triangular mesh
  for ( j=0; j<tmesh_.size(); j++ ) tmesh_[j].resize(3);

  GTPoint<GTICOS> a(3), b(3), c(3); // base vertices

  n = 0;
  GTVector<GTPoint<GTICOS>> R0(2*(nrows_+1)-1), R1(2*(nrows_+1));
  GTVector<GTPoint<GTICOS>> Rz(2*(nrows_+1)+1); // interleave R0, R1


#pragma omp parallel private (ibeg,l,m,t,a,b,c,R0,R1,Rz) reduction(+: n)
  R0.resize(2*(nrows_+1)-1);
  R1.resize(2*(nrows_+1));
  Rz.resize(2*(nrows_+1)+1);

  // Do refinement of base mesh triangles:
#pragma omp for
  for ( t=0; t<tbase_.size(); t++ ) { // for each base triangle 
    a = tbase_[t].v1; b = tbase_[t].v2; c = tbase_[t].v3;
    for ( l=0; l<nrows_+1; l++ ) { // for each triangle 'row'
      lagvert<GTICOS>(a,b,c,l   ,R0); // get previous row of points
      lagvert<GTICOS>(a,b,c,l+1 ,R1); // get current row of points
      interleave<GTICOS>(R0, R1, l, Rz); // zig-zag from R1 to R0, to R1, etc
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
  project2sphere<GTICOS>(tmesh_,1.0);

  // Order triangles (set iup_ flags):
  order_triangles<GTICOS>(tmesh_);

  // Compute centroids of all triangles:
  ftcentroids_.clear();
  ftcentroids_.resize(tmesh_.size());
  GTICOS fact = 1.0/3.0;
  for ( j=0; j<tmesh_.size(); j++ ) { // for each triangle
    a =  tmesh_[j].v1 + tmesh_[j].v2;
    a += tmesh_[j].v3;
    a *= fact;
    ftcentroids_[j] = a;
  }

} // end of method lagrefine


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
  GEOFLOW_TRACE();
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
  GEOFLOW_TRACE();
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
	GEOFLOW_TRACE();
  GString           serr = "GridIcos::do_elems2d (1): ";
  GFTYPE            fact;
  GTVector<GTPoint<GTICOS>> cverts(4), gverts(4), tverts(4);
  GTPoint<GTICOS>    ct(3), v1(3), v2(3), v3(3); // 3d points
  GElem_base        *pelem;
  GTVector<GINT>    iind;
  GTVector<GINT>    I(1);
  GTVector<GTICOS>  Ni;
  GTVector<GTVector<GFTYPE>>   *xNodes;
  GTVector<GTVector<GFTYPE>*>  *xiNodes;
  GTVector<GTVector<GTICOS>>    xid;
  GTVector<GTVector<GTICOS>*>   pxid;
  GTVector<GTVector<GTICOS>>    xgtmp(3);

  // Do eveything on unit sphere, then project to radiusi_
  // as a final step.

  assert(gbasis_.size()>0 && "Must set basis first");

  if ( gdd_ == NULLPTR ) gdd_ = new GDD_base<GTICOS>(nprocs_);

  // Resize points to appropriate size:
  for ( auto j=0; j<tmesh_.size(); j++ ) tmesh_[j].resize(3);
  for ( auto j=0; j<4; j++ ) {
    cverts[j].resize(3); // is a 3d point
    gverts[j].resize(3); // is only a 2d point
    tverts[j].resize(2); // is only a 2d point
  }

  // Get indirection indices to create elements
  // only for this task:
  gdd_->doDD(ftcentroids_, irank, iind);


  GTVector<GSIZET> isort;


  // When setting elements, must first construct Cartesian
  // coordinates at interior node points. This is done in
  // following steps:
  //   (0) find element vertices from triangular mesh, each set forms a plane
  //   (1) order vertices s.t. they're consistent with reference intervals 
  //   (2) construct interior node points planar vertices and basis/shape fcns
  //   (3) project the node coords to sphere in Cart. coords 
  fact = 1.0/3.0;
  GSIZET i;
  GSIZET nfnodes;   // no. face nodes
  GSIZET icurr = 0; // current global index
  GSIZET fcurr = 0; // current global face index
  // For each triangle in base mesh owned by this rank...
  for ( auto n=0; n<iind.size(); n++ ) { 
    i = iind[n];
    v1 = *tmesh_[i].v[0];
    v2 = *tmesh_[i].v[1];
    v3 = *tmesh_[i].v[2];
    ct = (v1 + (v2 + v3)) * fact;  // triangle centroid; don't overwrite later on
    // Compute element vertices:
    // NOTE: Is this ordering compatible with shape functions 
    // (placement of +/-1 with physical point)?
    for ( auto j=0; j<3; j++ ) { // 3 new elements for each triangle
      pelem = new GElem_base(2, GE_2DEMBEDDED, gbasis_);
      cverts[0] = v1; cverts[1] = (v1+v2)*0.5; cverts[2] = ct; cverts[3] = (v1+v3)*0.5;
#if 0
      order_latlong2d<GTICOS>(cverts);
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
      xid.resize(xiNodes->size());
      pxid.resize(xiNodes->size());
      for ( auto l=0; l<xid.size(); l++ ) pxid[l] = &xid[l];
      copycast<GFTYPE,GTICOS>(*xiNodes, pxid);
      
      Ni.resize(pelem->nnodes());

      project2sphere<GTICOS>(cverts, 1.0); // project verts to unit sphere     
      reorderverts2d<GTICOS>(cverts, tverts, isort, gverts); // reorder vertices consistenet with shape fcn
      for ( auto l=0; l<gverts[0].dim(); l++ ) { // loop over available coords
        xgtmp[l].resizem(pelem->nnodes());
        xgtmp[l] = 0.0;
        for ( auto m=0; m<4; m++ ) { // loop over vertices
          I[0] = m;
          lshapefcn_->Ni(I, pxid, Ni);
          xgtmp[l] += (Ni * (gverts[m][l]*0.25)); // node coordinate
        }
      }
      project2sphere<GTICOS>(xgtmp, radiusi_);
      for ( auto l=0; l<xgtmp.size(); l++ ) {
        assert( xgtmp[l].isfinite() );
        for ( auto k=0; k<xgtmp[l].size(); k++ ) (*xNodes)[l][k] = xgtmp[l][k];
      }
      pelem->init(*xNodes);

      nfnodes = 0;
      for ( auto j=0; j<pelem->nfaces(); j++ )  // get # face nodes
        nfnodes += pelem->face_indices(j).size();
      pelem->igbeg() = icurr;      // beginning global index
      pelem->igend() = icurr+pelem->nnodes()-1; // end global index
      pelem->ifbeg() = fcurr;
      pelem->ifend() = fcurr+nfnodes-1; // end global face index
      icurr += pelem->nnodes();
      fcurr += nfnodes;

      gelems_.push_back(pelem);

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
	GEOFLOW_TRACE();
  GString           serr = "GridIcos::do_elems3d (1): ";
  GSIZET            nxy;
  GTICOS            fact, r0, rdelta, xlatc, xlongc;
  GTVector<GTPoint<GTICOS>>
                    cverts(4), gverts(4), tverts(4);
  GTPoint<GTICOS>   ct(3), v1(3), v2(3), v3(3); // 3d points
  GElem_base        *pelem;
  GElem_base        *pelem2d;
  GTVector<GINT>    iind;
  GTVector<GINT>    I(1);
  GTVector<GTICOS>  Ni;
  GTVector<GTVector<GFTYPE>>  *xNodes;
  GTVector<GTVector<GFTYPE>>   xNodes2d(2);
  GTVector<GTVector<GFTYPE>*>  xiNodes2d(2);
  GTVector<GFTYPE>            *xiNodesr;
  GTVector<GTVector<GTICOS>>   xd, xd2d;
  GTVector<GTVector<GTICOS>>   xid, xid2d;
  GTVector<GTVector<GTICOS>*>  pxid;
  GTVector<GTVector<GTICOS>>   xgtmp(3);

  // Do eveything on unit sphere, then project to radiusi_
  // as a final step.

  assert(gbasis_.size()>0 && "Must set basis first");

  if ( gdd_ == NULLPTR ) gdd_ = new GDD_base<GTICOS>(nprocs_);

  // Resize points to appropriate size:
  for ( auto j=0; j<tmesh_.size(); j++ ) tmesh_[j].resize(3);
  for ( auto j=0; j<4; j++ ) {
    cverts[j].resize(3); // is a 3d point
    gverts[j].resize(3); // is only a 2d point
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
  for ( auto n=0; n<iind.size(); n++ ) { 
    i = iind[n];
    copycast<GTICOS,GTICOS>(*tmesh_[i].v[0], v1);
    copycast<GTICOS,GTICOS>(*tmesh_[i].v[1], v2);
    copycast<GTICOS,GTICOS>(*tmesh_[i].v[2], v3);
    ct = (v1 + (v2 + v3)) * fact;  // triangle centroid; don't overwrite later on
    // Compute element vertices:
    // NOTE: Is this ordering compatible with shape functions 
    // (placement of +/-1 with physical point)?
    for ( auto j=0; j<3; j++ ) { // 3 new elements for each triangle
      pelem2d = new GElem_base(2, GE_2DEMBEDDED, gbasis_);
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
      for ( auto m=0; m<2; m++ ) {
        xiNodes2d[m] = &pelem2d->xiNodes(m);
        xNodes2d [m].resizem(nxy);
      }
      xd2d .resize(xNodes2d.size());
      xid2d.resize(xiNodes2d.size());
      pxid .resize(xiNodes2d.size());
      for ( auto l=0; l<xid2d.size(); l++ ) pxid[l] = &xid2d[l];
      copycast<GFTYPE,GTICOS>(xiNodes2d, pxid);

      project2sphere<GTICOS>(cverts, 1.0); // project verts to unit sphere     
      reorderverts2d<GTICOS>(cverts, tverts, isort, gverts); // reorder vertices consistenet with shape fcn
      for ( auto l=0; l<gverts[0].dim(); l++ ) { // loop over available coords
        xgtmp[l].resizem(nxy);
        xgtmp[l] = 0.0;
        for ( auto m=0; m<4; m++ ) { // loop over vertices
          I[0] = m;
          lshapefcn_->Ni(I, pxid, Ni);
          xgtmp[l] += (Ni * (gverts[m][l]*0.25)); // node coordinate
        }
      }

      // Project to surface of inner sphere:
      project2sphere<GTICOS>(xgtmp, radiusi_);
      xyz2spherical<GTICOS>(xgtmp);

      // Loop over radial elements and build all elements 
      // based on this patch (we're still in sph coords here):
      for ( auto e=0; e<nradelem_; e++ ) {
        pelem = new GElem_base(GDIM, GE_DEFORMED, gbasis_);
        xiNodesr = &pelem->xiNodes(2); // get radial reference nodes
        xNodes   = &pelem->xNodes();  // node spatial data
        r0       = radiusi_ + e*rdelta;
        for ( auto m=0; m<pelem->size(2); m++ ) { // for each radial node
          for ( auto n=0; n<nxy; n++ ) { // find sph coords for each horizontal node
            (*xNodes)[0][n+m*nxy] =  0.5*rdelta*((*xiNodesr)[m] + 1.0);
            (*xNodes)[1][n+m*nxy] =  xgtmp[1][n];
            (*xNodes)[2][n+m*nxy] =  xgtmp[2][n];
          }
        }
        spherical2xyz<GFTYPE>(*xNodes);

        pelem->init(*xNodes);

        nfnodes = 0;
        for ( auto j=0; j<pelem->nfaces(); j++ )  // get # face nodes
          nfnodes += pelem->face_indices(j).size();
        pelem->igbeg() = icurr;      // beginning global index
        pelem->igend() = icurr+pelem->nnodes()-1; // end global index
        pelem->ifbeg() = fcurr;
        pelem->ifend() = fcurr+nfnodes-1; // end global face index
        icurr += pelem->nnodes();
        fcurr += nfnodes;

        gelems_.push_back(pelem);
      }
      delete pelem2d;

    } // end of element loop for this triangle
  } // end of triangle base mesh loop

} // end of method do_elems3d (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : do_elems2d (2)
// DESC   : Build 2d elemental grid on base mesh. Meant to be used for restart,
//          where nodal grid data is already known.
// ARGS   : p      : matrix of size the number of elements X GDIM containing 
//                   the poly expansion order in each direction
//          gxnodes: vector of GDIM vectors containing Cartesian coords of elements
// RETURNS: none.
//**********************************************************************************
void GGridIcos::do_elems2d(GTMatrix<GINT> &p,
                           GTVector<GTVector<GFTYPE>> &gxnodes)
{
	GEOFLOW_TRACE();
  GString                     serr = "GridIcos::do_elems2d (2): ";
  GElem_base                  *pelem;
  GTVector<GTVector<GFTYPE>>  *xNodes;
  GTVector<GNBasis<GCTYPE,GFTYPE>*>
                               gb(GDIM);
  GTVector<GINT>               ppool(gbasis_.size());

  // Now, treat the gbasis_ as a pool that we search
  // to find bases we need:
  for ( auto j=0; j<ppool.size(); j++ ) ppool[j] = gbasis_[j]->getOrder();


  // Set element internal dof from input data:
  GSIZET iwhere ;
  GSIZET nvnodes;   // no. vol nodes
  GSIZET nfnodes;   // no. face nodes
  GSIZET icurr = 0; // current global index
  GSIZET fcurr = 0; // current global face index
  // For each triangle in base mesh owned by this rank...
  for ( auto i=0; i<p.size(1); i++ ) { 
    nvnodes = 1;
    for ( auto j=0; j<GDIM; j++ ) { // set basis from pool
      assert(ppool.contains(p(i,j),iwhere) && "Expansion order not found");
      gb[j] = gbasis_[iwhere];
      nvnodes *= (p(i,j) + 1);
    }
    pelem = new GElem_base(GDIM, GE_2DEMBEDDED, gb);
    xNodes  = &pelem->xNodes();  // node spatial data

    // Set internal node positions from input data.
    // Note that gxnodes are 'global' and xNodes is
    // element-local:
    for ( auto j=0; j<GDIM; j++ ) {
       gxnodes[j].range(icurr, icurr+nvnodes-1);
      (*xNodes)[j] = gxnodes[j];
    }
    for ( auto j=0; j<GDIM; j++ ) gxnodes[j].range_reset();

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
//          prior to entry, and that the base icos grid has beed set. Meant
//          to be used on restart, where nodal grid data is already known.
// ARGS   : p      : matrix of size the number of elements X GDIM containing 
//                   the poly expansion order in each direction
//          gxnodes: vector of GDIM vectors containing Cartesian coords of elements
// RETURNS: none.
//**********************************************************************************
void GGridIcos::do_elems3d(GTMatrix<GINT> &p,
                           GTVector<GTVector<GFTYPE>> &gxnodes)
{
	GEOFLOW_TRACE();
  GString                      serr = "GridIcos::do_elems3d (2): ";
  GElem_base                  *pelem;
  GTVector<GTVector<GFTYPE>>  *xNodes;
  GTVector<GTVector<GFTYPE>*> *xiNodes;
  GTVector<GNBasis<GCTYPE,GFTYPE>*>
                               gb(GDIM);
  GTVector<GINT>               ppool(gbasis_.size());

  // Now, treat the gbasis_ as a pool that we search
  // to find bases we need:
  for ( auto j=0; j<ppool.size(); j++ ) ppool[j] = gbasis_[j]->getOrder();


  // Set element internal dof from input data:
  GSIZET iwhere ;
  GSIZET nvnodes;   // no. vol nodes
  GSIZET nfnodes;   // no. face nodes
  GSIZET icurr = 0; // current global index
  GSIZET fcurr = 0; // current global face index
  // For each triangle in base mesh owned by this rank...
  for ( auto i=0; i<p.size(1); i++ ) { 
    nvnodes = 1;
    for ( auto j=0; j<GDIM; j++ ) { // set basis from pool
      assert(ppool.contains(p(i,j),iwhere) && "Expansion order not found");
      gb[j] = gbasis_[iwhere];
      nvnodes *= (p(i,j) + 1);
    }
    pelem = new GElem_base(GDIM, GE_DEFORMED, gb);
    xNodes  = &pelem->xNodes();  // node spatial data

    // Set internal node positions from input data.
    // Note that gxnodes are 'global' and xNodes is
    // element-local:
    for ( auto j=0; j<GDIM; j++ ) {
       gxnodes[j].range(icurr, icurr+nvnodes-1);
      (*xNodes)[j] = gxnodes[j];
    }
    for ( auto j=0; j<GDIM; j++ ) gxnodes[j].range_reset();

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
  GEOFLOW_TRACE();
  GString serr = "GridIcos::print: ";
  std::ofstream ios;

  GTPoint<GTICOS> pt;
  GTICOS          r, xlat, xlong;

  ios.open(filename);
  if ( icoord == GICOS_LATLONG) { // print in lat-long
    for ( auto i=0; i<tmesh_.size(); i++ ) { // for each triangle
      for ( auto j=0; j<3; j++ ) { // for each vertex of triangle
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
    for ( auto i=0; i<tmesh_.size(); i++ ) { // for each triangle
      for ( auto j=0; j<3; j++ ) { // for each vertex of triangle
        pt = *tmesh_[i].v[j];
        ios << pt.x1 << " " << pt.x2 << " " << pt.x3 << std::endl ;
      }
    }
  }

  ios.close();

} // end of method print


//**********************************************************************************
//**********************************************************************************
// METHOD : config_gbdy
// DESC   : Configure 3d spherical boundaries from ptree
// ARGS   : 
//          ptree : main prop tree 
//          igbdyf : For each natural/canonical global boundary face,
//                  gives vector of global bdy ids. Allocated here.
//          igbdyft: bdy type ids for each index in igbdyf. Allocated here.
//          igbdy  : 'flat' version of igbdyf
// RETURNS: none.
//**********************************************************************************
void GGridIcos::config_gbdy(const geoflow::tbox::PropertyTree &ptree, 
                           GTVector<GTVector<GSIZET>>   &igbdyf, 
                           GTVector<GTVector<GBdyType>> &igbdyft,
                           GTVector<GSIZET>             &igbdy)
{
	GEOFLOW_TRACE();
  // Cycle over all geometric boundaries, and configure:

  GBOOL              bret;
  GSIZET             nind;
  GTVector<GSIZET>   itmp;
  GTVector<GFTYPE>   rbdy(2);
  GTVector<GString>  bdynames(2);
  std::vector<GString>  svec;
  GString            gname, sbdy, bdyclass;
  PropertyTree       bdytree, gridptree;
  UpdateBasePtr      base_ptr;
  stBdyBlock         stblock;

  // Clear input arrays:
  igbdyf .clear();
  igbdyft.clear();

  if ( ndim_ == 2 ) return; // no boundaries to configure
 
  bdynames[0] = "bdy_inner";
  bdynames[1] = "bdy_outer";

  if ( !ptree.isValue<GString>("grid_type") ) {
    cout << "GGridIcos::config_gbdy: grid_type not set" << endl;
    assert(FALSE);
  }
  gname     = ptree.getValue<GString>("grid_type");

  if ( !ptree.isPropertyTree(gname) ) {
    cout << "GGridIcos::config_gbdy: grid_type block " << gname << " not found" << endl;
    assert(FALSE);
  }
  gridptree = ptree.getPropertyTree(gname);

  rbdy[0] = radiusi_;
  rbdy[1] = radiuso_;

  igbdyf.resize(2); // 2 canonical bdys
  igbdyft.resize(2); // 2 canonical bdys

  bdy_update_list_.resize(2);
 

  // Get properties from the main prop tree. 
  // Note: bdys are configured by way of geometry's
  //       natural decomposition: here, by inner and
  //       outer spherical surfaces. But the bdy indices 
  //       and types returned on exit contain info for all bdys:
  for ( auto j=0; j<2; j++ ) { // cycle over 2 spherical surfaces
    sbdy         = gridptree.getValue<GString>(bdynames[j]);
    bdytree      = ptree.getPropertyTree(sbdy);
    bdyclass     = bdytree.getValue<GString>("bdy_class", "uniform");
    find_gbdy_ind3d(rbdy[j], itmp); 
    do_gbdy_normals(j, FALSE, dXdXi_, itmp, bdyNormals_, idepClomp_); // all bdy nodes 

    igbdyf [j].resize(itmp.size()); igbdyf [j] = itmp;
    igbdyft[j].resize(itmp.size()); igbdyft[j] = GBDY_NONE;
    nind = 0;
    for ( auto k=0; k<igbdyf.size(); k++ ) nind += igbdyf[k].size();
    igbdy  .resize(nind); // vol indices of bdy nodes in base; bdy update needs this
    nind = 0;
    for ( auto i=0; i<igbdyf[j].size(); i++ ) igbdy[nind++] = igbdyf[j][i];

    geoflow::get_bdy_block(bdytree, stblock);
    if ( "uniform" == bdyclass ) { // uniform bdy conditions
      assert(!stblock.tbdy.contains(GBDY_PERIODIC) && "Invalid boundary condition");
      // May have different uniform bdys for different state comps:
      for ( auto k=0; k<stblock.tbdy.size(); k++ ) { 
        base_ptr = GUpdateBdyFactory<BdyTypePack>::build(ptree, sbdy, *this,  j, 
                                            stblock.tbdy[k], stblock.istate[k], stblock.value[k], itmp);
        if ( stblock.tbdy[k] != GBDY_NONE ) igbdyft[j] = stblock.tbdy[k];
        bdy_update_list_[j].push_back(base_ptr);
      }
    }
    else if ( "mixed" == bdyclass ) { // mixed bdy conditions
      assert( bdytree.isArray<GString>("bdy_blocks") && "No bdy_blocks specified"); 
      svec = bdytree.getArray<GString>("bdy_blocks");
      for ( auto i=0; i<svec.size(); i++ ) {
        bdytree = ptree.getPropertyTree(svec[i]);
        geoflow::get_bdy_block(bdytree, stblock);
        assert(!stblock.tbdy.contains(GBDY_PERIODIC) && "Invalid boundary condition");
        GSpecBdyFactory::dospec(bdytree, *this, j, itmp);
        for ( auto k=0; k<svec.size(); k++ ) { // for each sub-block
          base_ptr = GUpdateBdyFactory<BdyTypePack>::build(ptree, svec[k], *this, j, 
                                              stblock.tbdy[k], stblock.istate[k], stblock.value[k], itmp);
          
          for ( auto m=0; m<itmp.size(); m++ ) {
            if ( igbdyf[j].contains(itmp[m]) ) igbdyft[j][m] = stblock.tbdy[k];
          }
          if ( stblock.tbdy[k] != GBDY_NONE ) igbdyft[j] = stblock.tbdy[k];
          bdy_update_list_[j].push_back(base_ptr);
        }
      } 
    }
    else {
      assert(FALSE && "Invalid bdy_class");
    }
  } // end, canonical bdy loop

} // end of method config_gbdy


//**********************************************************************************
//**********************************************************************************
// METHOD : find_gbdy_ind3d
// DESC   : Find global bdy indices (indices into xNodes_ arrays) that
//          corresp to specified radius
// ARGS   : radius   : radius
//          ibdy     : array of indices into xNodes that comprise this boundary
// RETURNS: none.
//**********************************************************************************
void GGridIcos::find_gbdy_ind3d(GFTYPE radius, GTVector<GSIZET> &ibdy)
{
	GEOFLOW_TRACE();
  GFTYPE          eps, r;
  GTPoint<GFTYPE> pt(ndim_);

  ibdy.clear();
  eps = 1.0e4*std::numeric_limits<GFTYPE>::epsilon();

  for ( auto i=0; i<xNodes_[0].size(); i++ ) { 
      r = sqrt(pow(xNodes_[0][i],2)+pow(xNodes_[1][i],2)+pow(xNodes_[2][i],2));
      if ( FUZZYEQ(r, radius, eps) ) {
        ibdy.push_back(i);
      }
  }

} // end, method find_gbdy_ind3d


//**********************************************************************************
//**********************************************************************************
// METHOD : elem_bdy_data
// DESC   : Compute normals to each element face
// ARGS   : 
//          dXdXi     : matrix of dX_i/dXi_j matrix elements, s.t.
//                      dXdX_i(i,j) = dx^j/dxi^i
//          gieface   : vector of face indices into global volume fields 
//                      for all facase
//          gdeface   : description for each face node
//          face_mass : mass on eleme faces
//          normals   : vector of normal components
//          depComp   : vector index dependent on the other indices (first 
//                      component index whose normal component is nonzero)
// RETURNS: none
//**********************************************************************************
void GGridIcos::elem_bdy_data(const GTMatrix<GTVector<GFTYPE>> &dXdXi,
                              GTVector<GSIEZT>                 &igeface,
                              GTVector<GFTYPE>                 &face_mass,
                              GTVector<GTVector<GFTYPE>>       &normals)
{

	GEOFLOW_TRACE();
  #if defined(_G_IS2D)
    elem_bdy_data2d(dXdXi, gieface, gdeface, face_mass, normals, depComp);
  #elif defined(_G_IS3D)
    elem_bdy_data3d(dXdXi, gieface, gdeface, face_mass, normals, depComp);
  #else
    #error Invalid problem dimensionality
  #endif

} // end, method elem_bdy_data


//**********************************************************************************
//**********************************************************************************
// METHOD : elem_bdy_data2d
// DESC   : Compute normals to each element face in 2d
// ARGS   : 
//          dXdXi     : matrix of dX_i/dXi_j matrix elements, s.t.
//                      dXdX_i(i,j) = dx^j/dxi^i
//          gieface   : vector of face indices into global volume fields 
//                      for all facase
//          gdeface   : description for each face node
//          face_mass : mass on eleme faces
//          normals   : vector of normal components
//          depComp   : vector index dependent on the other indices (first 
//                      component index whose normal component is nonzero)
// RETURNS: none
//**********************************************************************************
void GGridIcos::elem_bdy_data2d(GTMatrix<GTVector<GFTYPE>> &dXdXi,
                                  GTVector<GSIZET>           &gieface,
                                  GTVector<GUINT>            &gdeface,
                                  GTVector<GFTYPE>           &face_mass,
                                  GTVector<GTVector<GFTYPE>> &normals,
                                  GTVector<GINT>             &depComp)
{
	GEOFLOW_TRACE();
   GINT              ib, ic, id,  ip;
   GFTYPE            tiny;
   GFTYPE            xm;
   GTPoint<GFTYPE>   kp(3), xp(3), p1(3), p2(3);

   tiny  = 100.0*std::numeric_limits<GFTYPE>::epsilon();
   kp    = 1.0/sqrt(3.0); // unit vector in radial direction


   if ( this->gtype_ == GE_DEFORMED
    ||  this->gtype_ == GE_2DEMBEDDED ) {
     // Bdy normal is hat{k} X dvec{X} / dxi_iedge,
     // for face:
     for ( auto j=0; j<gieface.size(); j++ ) { // all points on iedge
       ib = gieface  [j];
       id = NDHOST(gdeface[j]);
       ip = id%2;
       xm = id == 2 || id == 3 ? 1.0 : -1.0;
       for ( auto i=0; i<dXdXi.size(2); i++ ) { // over _X_
         p1[i] = dXdXi(ip,i)[ib];
       }
       kp.cross(p1, xp);   // xp =  p1 X k
       face_mass  [j] *= xp.mag(); // |k X d_X_/dxi_ip| is face Jac
       xp *= xm; xp.unit();
       for ( ic=0; ic<GDIM && fabs(xp[ic]) < tiny; ic++ );
       assert(ic < GDIM); // no normal components > 0
       for ( auto i=0; i<normals.size(); i++ ) normals[i][j] = xp[i];
     }
   }
   else {
     assert(FALSE && "Invalid grid type");
   }

} // end, method elem_bdy_data2d


//**********************************************************************************
//**********************************************************************************
// METHOD : elem_bdy_data3d
// DESC   : Compute normals to each element face in 3d
// ARGS   : 
//          dXdXi     : matrix of dX_i/dXi_j matrix elements, s.t.
//                      dXdX_i(i,j) = dx^j/dxi^i
//          gieface   : vector of face indices into global volume fields 
//                      for all facase
//          face_mass : mass on eleme faces
//          gdeface   : description for each face node
//          normals   : vector of normal components
//          depComp   : vector index dependent on the other indices (first 
//                      component index whose normal component is nonzero)
// RETURNS: none
//**********************************************************************************
void GGridIcos::elem_bdy_data3d(GINT                             ibdy,
                                  GBOOL                            bunique,
                                  const GTMatrix<GTVector<GFTYPE>> &dXdXi,
                                  const GTVector<GSIZET>           &igbdy,
                                  GTVector<GTVector<GFTYPE>>       &normals,
                                  GTVector<GINT>                   &idepComp)
{
   GEOFLOW_TRACE();

   GINT            ib, ic, id;
   GINT            ixi[6][2] = { {0,2}, {1,2}, {2,0},
                                 {2,1}, {1,0}, {0,1} };
   GFTYPE          tiny;
   GTPoint<GFTYPE> xp(3), p1(3), p2(3);

   tiny  = 100.0*std::numeric_limits<GFTYPE>::epsilon();

   if ( this->gtype_ == GE_DEFORMED
    ||  this->gtype_ == GE_2DEMBEDDED ) {
     // Bdy normal is dvec{X} / dxi_xi X dvec{X} / dxi_eta
     for ( auto j=0; j<gieface.size(); j++ ) { // all points on iedge
       ib = gieface[j];
       id = gdeface[j];
       xm = id == 2 || id == 3 ? 1.0 : -1.0;
       // Find derivs of _X_ wrt face's reference coords;
       // the cross prod of these vectors is the normal:
       for ( auto i=0; i<dXdXi.size(2); i++ ) { // d_X_/dXi
         p1[i] = dXdXi(ixi[j][0],i)[ib]; // d_X_/dxi
         p2[i] = dXdXi(ixi[j][1],i)[ib]; // d_X_/deta
       }
       p1.cross(p2, xp);   // xp = p1 X p2
       face_mass  [j] *= xp.mag(); // d_X_/dxi X d_X_/deta| is face Jac
       xp.unit();
       for ( ic=0; ic<GDIM && fabs(xp[ic]) < tiny; ic++ );
       assert(ic < GDIM); // no normal components > 0
       for ( auto i=0; i<normals.size(); i++ ) normals[i][j] = xp[i];
     }
   }
   else {
     assert(FALSE && "Invalid grid type");
   }

} // end, method elem_bdy_data3d


//**********************************************************************************
//**********************************************************************************
// METHOD : do_gbdy_normals
// DESC   : Compute normals to each domain bdy 
// ARGS   : 
//          igbdy    : vector of bdy indices into global volume fields 
//          debdy    : array of node 'descriptions', with dimension of igbdy
//          dXdXi    : matrix of dX_i/dXi_j matrix elements, s.t.
//                     dXdX_i(i,j) = dx^j/dxi^i
//          normals  : vector of normal components, each of dim of igbdy
//          depComp  : vector index dependent on the other indices (first 
//                     component index whose normal component is nonzero)

// RETURNS: none
//**********************************************************************************
void GGridIcos::do_gbdy_normals(const GTVector<GSIZET>           &igbdy,
                                const GTVector<GUINT>            &debdy,
                                const GTMatrix<GTVector<GFTYPE>> &dXdXi,
                                GTVector<GTVector<GFTYPE>>       &normals,
                                GTVector<GINT>                   &idepComp)
{
  GEOFLOW_TRACE();
  GSIZET nbdy;

  #if defined(_G_IS2D)
    return;
  #if defined(_G_IS3D)
    nbdy = igbdy.size();
    idepComp.resize(nbdy);
    for ( auto j=0; j<normals.size(); j++ ) normals[j].resize(nbdy);

    do_gbdy_normals3d(igbdy, debdy, dXdXi, normals, idepComp);

  #else
    #error Invalid problem dimensionality
  #endif



} // end, method do_gbdy_normals


//**********************************************************************************
//**********************************************************************************
// METHOD : do_gbdy_normals3d
// DESC   : Compute normals to each domain bdy in 2d. Also
//          provide the vector component index for the dependent
//          component, so that we can easily enforce, say, the constraint
//                hat(n) \cdot \vec{u} = 0
//          for the correct \vec{u} component. We must
//          keep in mind that there may be terrain on the boundary,
//          so we cannot rely on alignment of bdy surface with
//          coordinate directions. This method should be called 
//          after terrain is added.
// ARGS   : 
//          igbdy    : vector of bdy indices into global volume fields 
//          debdy    : array of node 'descriptions', with dimension of igbdy
//          dXdXi    : matrix of dX_i/dXi_j matrix elements, s.t.
//                     dXdX_i(i,j) = dx^j/dxi^i
//          normals  : vector of normal components, each of dim of igbdy
//          depComp  : vector index dependent on the other indices (first 
//                     component index whose normal component is nonzero)

// RETURNS: none
//**********************************************************************************
void GGridIcos::do_gbdy_normals(const GTVector<GSIZET>           &igbdy,
                                const GTVector<GUINT>            &debdy,
                                const GTMatrix<GTVector<GFTYPE>> &dXdXi,
                                GTVector<GTVector<GFTYPE>>       &normals,
                                GTVector<GINT>                   &depComp)
{
  GEOFLOW_TRACE();
  GSIZET          ib, ic, ip;
  GUINT           id;
  GFTYPE          tiny;
  GFTYPE          xm;
  GTPoint<GFTYPE> xp(3), p1(3), p2(3);

  tiny  = 100.0*std::numeric_limits<GFTYPE>::epsilon();

  if ( this->gtype_ == GE_DEFORMED ) {
     // Bdy normal is dvec{X} / dxi_xi X dvec{X} / dxi_eta
     for ( auto j=0; j<igbdy.size(); j++ ) { // all points on face
       ib = igbdy[j];
       id = GET_NDHOST(debdy[j]); // host face id
       xm = id == 1 || id == 2 || id == 5 ? 1.0 : -1.0;
       for ( auto i=0; i<dXdXi.size(2); i++ ) { // over _X_
         p1[i] = dXdXi(0,i)[ib]; // d_X_/dxi
         p2[i] = dXdXi(1,i)[ib]; // d_X_/deta
       }
       p1.cross(p2, xp);   // xp = p1 X p2
       xp *= xm;
       xp.unit();
       for ( ic=0; ic<xp.dim(); ic++ ) if ( fabs(xp[ic]) > tiny ) break;
       assert(ic >= GDIM); // no normal components > 0
       for ( auto i=0; i<normals.size(); i++ ) normals[i][j] = xp[i];
       depComp[j] = ic;  // dependent component
     }
   }
   else {
     assert(FALSE && "Invalid grid type");
   }

} // end, method do_gbdy_normals3d

