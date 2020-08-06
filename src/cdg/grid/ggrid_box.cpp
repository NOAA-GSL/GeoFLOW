//==================================================================================
// Module       : ggrid_box
// Date         : 11/11/18 (DLR)
// Description  : Object defining a (global) 2d or 3d box grid. Builds
//                elements that base class then uses to build global
//                computational data structures.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#include <cstdlib>
#include <memory>
#include <cmath>
#include "gcomm.hpp"
#include "ggrid_box.hpp"
#include "gspecbdy_factory.hpp"
#include "gupdatebdy_factory.hpp"
#include "ginitstate_factory.hpp"
#include "gmtk.hpp"
#include "gutils.hpp"
#include "tbox/mpixx.hpp"
#include "tbox/global_manager.hpp"


using namespace geoflow;
using namespace geoflow::tbox;
using namespace std;

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Instantiate with box lengths (x, y, z), directional basis, and
//          domain decompoisition object. This generates a 2d or 3d grid
//          depending on whether b.size = 2 or 3..
// ARGS   : ptree  : main proptery tree
//          b      : vector of basis pointers. Number of elements in b is
//                   what _must_ be in L, and ne.
//          comm   : communicator
// RETURNS: none
//**********************************************************************************
GGridBox::GGridBox(const geoflow::tbox::PropertyTree &ptree, GTVector<GNBasis<GCTYPE,GFTYPE>*> &b, GC_COMM &comm)
:          GGrid(ptree, b, comm),
ndim_                     (GDIM),
gdd_                   (NULLPTR),
lshapefcn_             (NULLPTR)
{
  assert((b.size() == GDIM ) 
        && "Basis has incorrect dimensionalilty");

  irank_  = GComm::WorldRank(comm_);
  nprocs_ = GComm::WorldSize(comm_);

  GINT    inorm;
  GString gname   = ptree.getValue<GString>("grid_type");
  GString tname   = ptree.getValue<GString>("terrain_type");
  GString snorm;
  assert(gname == "grid_box");
  geoflow::tbox::PropertyTree gridptree = ptree.getPropertyTree(gname);

  // If terrain is being used, elements may not be
  // GE_REGULAR:
  this->gtype_ = GE_REGULAR;
  if ( "none" != tname 
   &&  ""     != tname ) this->gtype_ = GE_DEFORMED;

  gbasis_.resize(GDIM);
  gbasis_ = b;
  Lbox_.resize(GDIM);
  ne_.resize(GDIM);

  GTPoint<GFTYPE> gp(ndim_); 
  std::vector<GFTYPE> spt(3); // tmp array
  std::vector  <int> sne   ; // tmp array
  spt = gridptree.getArray<GFTYPE>("xyz0");
  P0_.resize(GDIM);
  P1_.resize(GDIM);
  dP_.resize(GDIM);
  for ( auto j=0; j<GDIM; j++ ) P0_[j] = spt[j];
  spt = gridptree.getArray<GFTYPE>("delxyz");
  sne = gridptree.getArray<int>("num_elems");
  this->cgtraits_.maxit = gridptree.getValue<GDOUBLE>("maxit");
  this->cgtraits_.tol   = gridptree.getValue<GDOUBLE>("tol");
  snorm                 = gridptree.getValue<GString>("norm_type");
  this->cgtraits_.normtype = LinSolverBase<CGTypePack>::str2normtype(snorm);

  eps_ = 100*std::numeric_limits<GFTYPE>::epsilon();
  // compute global bdy range, and global vertices:
  for ( auto j=0; j<GDIM; j++ ) dP_[j] = spt[j];
  P1_ = P0_ + dP_;
  gverts_.resize(pow(2,ndim_));
  for ( auto j=0; j<gverts_.size(); j++ ) gverts_[j].resize(GDIM);
  if ( ndim_ == 2 ) {
    gverts_[0] = P0_; 
    gp = P0_; gp.x1 += dP_.x1; gverts_[1] = gp; 
    gverts_[2] = P1_; 
    gp = P0_; gp.x2 += dP_.x2; gverts_[3] = gp;
  }
  else if ( ndim_ == 3 ) {
    gverts_[0] = P0_; 
    gp = P0_; gp.x1 += dP_.x1; gverts_[1] = gp; 
    gverts_[2] = P1_; 
    gp = P0_; gp.x2 += dP_.x2; gverts_[3] = gp;
    gp = P0_; gp.x3 += dP_.x3; gverts_[4] = gp; 
    gp.x1 += dP_.x1; gverts_[5] = gp; 
    gverts_[6] = P1_; 
    gp = P0_; gp.x3 += dP_.x3; gp.x2 += dP_.x2; gverts_[7] = gp;
  }

  ne_.resize(b.size());
  for( auto j=0; j<b.size(); j++ ) {
    Lbox_[j] = fabs(dP_[j]);
    ne_  [j] = sne[j];
  }

  lshapefcn_ = new GShapeFcn_linear<GFTYPE>();
  if ( GDIM == 2 ) {
    init2d();
  }
  else if ( GDIM == 3 ) {
    init3d();
  }

} // end of constructor method (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
GGridBox::~GGridBox()
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
std::ostream &operator<<(std::ostream &str, GGridBox &e)
{
  
  str << "    Lbox: " << e.Lbox_;
  str << std::endl << " Centroids: " ;
  for( auto i=0; i<e.ftcentroids_.size(); i++ )
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
void GGridBox::set_partitioner(GDD_base<GFTYPE> *gdd)
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
void GGridBox::init2d()
{

  find_subdomain();

} // end of method init2d


//**********************************************************************************
//**********************************************************************************
// METHOD : init3d
// DESC   : Initialize for 3d elements
// ARGS   : none.
// RETURNS: none.
//**********************************************************************************
void GGridBox::init3d()
{

  find_subdomain();

} // end, method init3d



//**********************************************************************************
//**********************************************************************************
// METHOD : do_elems (1)
// DESC   : Public entry point for grid element computation
// ARGS   : none.
// RETURNS: none.
//**********************************************************************************
void GGridBox::do_elems()
{
  if ( ndim_ == 2 ) do_elems2d();
  if ( ndim_ == 3 ) do_elems3d();

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
void GGridBox::do_elems(GTMatrix<GINT> &p,
                        GTVector<GTVector<GFTYPE>> &xnodes)
{
  if ( ndim_ == 2 ) do_elems2d(p, xnodes);
  if ( ndim_ == 3 ) do_elems3d(p, xnodes);

} // end, method do_elems (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : do_elems2d (1)
// DESC   : Build 2d element list. It's assumed that init2d has been called
//          prior to entry, and that the qmesh_ has beed set. This arrays
//          of hexagonal 3d elements provides for each hex its vertices in Cartesian
//          coordinates. Also, the centroids of these hexes (in Cartesian coords)
//          should also have been computed. The Cartesian hex vertices will be
//          converted into sherical coordinates, and the GShapeFcn_linear
//          will be used in to compute the interior nodes. 
// ARGS   : 
//          rank: MPI rank or partition id
// RETURNS: none.
//**********************************************************************************
void GGridBox::do_elems2d()
{
  assert(gbasis_.size()>0 && "Must set basis first");
  assert(ndim_ == 2 && "Dimension must be 2");

  GString                      serr = "GGridBox::do_elems2d (1): ";
  GTPoint<GFTYPE>              cent;
  GTVector<GINT>               iind;
  GTVector<GINT>               I(3);
  GTVector<GINT>              *bdy_ind;
  GTVector<GBdyType>          *bdy_typ;
  GTVector<GINT>              *face_ind;
  GTVector<GFTYPE>             Ni;
  GElem_base                  *pelem;
  GTVector<GTVector<GFTYPE>>  *xNodes;
  GTVector<GTVector<GFTYPE>*> *xiNodes;
  GTVector<GTVector<GFTYPE>>   xgtmp(3);


#if 0
  if ( gdd_       == NULLPTR ) gdd_ = new GDD_base(nprocs_);

  // Get indirection indices to create elements
  // only for this task:
  gdd_->doDD(ftcentroids_, irank_, iind);
#endif

  GSIZET nvnodes;   // no. vol nodes
  GSIZET nfnodes;   // no. face nodes
  GSIZET nbnodes;   // no. bdy nodes
  GLONG  icurr = 0; // current global index
  GLONG  fcurr = 0; // current global face index
  GLONG  bcurr = 0; // current global bdy index

  for ( auto i=0; i<qmesh_.size(); i++ ) { // for each quad in irank's mesh
    pelem = new GElem_base(this->gtype_, gbasis_);
    xNodes  = &pelem->xNodes();  // node spatial data
    xiNodes = &pelem->xiNodes(); // node ref interval data
    Ni.resize(pelem->nnodes()); // tensor product shape function
    bdy_ind = &pelem->bdy_indices(); // get bdy indices data member
    bdy_typ = &pelem->bdy_types  (); // get bdy types data member
    bdy_ind->clear(); bdy_typ->clear();
    for ( auto l=0; l<ndim_; l++ ) { // loop over element Cart coords
      (*xNodes)[l] = 0.0;
      for ( auto m=0; m<pow(2,ndim_); m++ ) { // loop over verts given in Cart coords
        I[0] = m;
        lshapefcn_->Ni(I, *xiNodes, Ni);
        (*xNodes)[l] += Ni * ( (*(qmesh_[i].v[m]))[l] * 0.25 );
      }
    }

    pelem->init(*xNodes);

    // With face/edge centroids computed, compute 
    // global boundary nodes:
    for ( auto j=0; j<2*ndim_; j++ ) { // cycle over all edges
      cent = pelem->edgeCentroid(j);
      face_ind = NULLPTR;
      if ( FUZZYEQ(P0_.x2,cent.x2,eps_) ) face_ind = &pelem->edge_indices(0);
      if ( FUZZYEQ(P1_.x1,cent.x1,eps_) ) face_ind = &pelem->edge_indices(1);
      if ( FUZZYEQ(P1_.x2,cent.x2,eps_) ) face_ind = &pelem->edge_indices(2);
      if ( FUZZYEQ(P0_.x1,cent.x1,eps_) ) face_ind = &pelem->edge_indices(3);
      // For now, we don't allow corner nodes to be repeated, 
      // and we'll have to choose the best way to define the 
      // normal vectors at the these nodes:
      for ( auto k=0; face_ind != NULLPTR && k<face_ind->size(); k++ ) {
        if ( !bdy_ind->contains((*face_ind)[k]) ) {
          bdy_ind->push_back((*face_ind)[k]); 
          bdy_typ->push_back(GBDY_NONE); 
        }
      }
    }


    // Find global global interior and bdy start & stop indices represented 
    // locally within element:
    nvnodes = pelem->nnodes();
    nfnodes = pelem->nfnodes();
    nbnodes = pelem->bdy_indices().size();
    pelem->igbeg() = icurr;           // beg global vol index
    pelem->igend() = icurr+nvnodes-1; // end global vol index
    pelem->ifbeg() = fcurr;           // beg global face index
    pelem->ifend() = fcurr+nfnodes-1; // end global face index
    pelem->ibbeg() = bcurr;           // beg global bdy index
    pelem->ibend() = bcurr+nbnodes-1; // end global bdy index
    icurr += nvnodes;
    fcurr += nfnodes;
    bcurr += nbnodes;

    gelems_.push_back(pelem);
  } // end of quad mesh loop

} // end of method do_elems2d (1)



//**********************************************************************************
//**********************************************************************************
// METHOD : do_elems3d (1)
// DESC   : Build 3d element list. It's assumed that init3d has been called
//          prior to entry, and that the qmesh_ has beed set. This arrays
//          of hexagonal 3d elements provides for each hex its vertices in Cartesian
//          coordinates. Also, the centroids of these hexes (in Cartesian coords)
//          should also have been computed. The Cartesian hex vertices will be
//          converted into sherical coordinates, and the GShapeFcn_linear
//          will be used to compute the interior nodes. 
// ARGS   : 
// RETURNS: none.
//**********************************************************************************
void GGridBox::do_elems3d()
{

  assert(gbasis_.size()>0 && "Must set basis first");
  assert(ndim_ == 3 && "Dimension must be 3");

  GString                      serr = "GGridBox::do_elems3d (1): ";
  GTPoint<GFTYPE>              cent;
  GTVector<GINT>               iind;
  GTVector<GINT>               I(3);
  GTVector<GINT>              *bdy_ind;
  GTVector<GBdyType>          *bdy_typ;
  GTVector<GINT>              *face_ind;
  GTVector<GFTYPE>             Ni;
  GElem_base                  *pelem;
  GTVector<GTVector<GFTYPE>>  *xNodes;
  GTVector<GTVector<GFTYPE>*> *xiNodes;
  GTVector<GTVector<GFTYPE>>   xgtmp(3);


#if 0
  if ( gdd_       == NULLPTR ) gdd_ = new GDD_base(nprocs_);

  // Get indirection indices to create elements
  // only for this task:
  gdd_->doDD(ftcentroids_, irank_, iind);
#endif

  GSIZET nvnodes;   // no. vol indices
  GSIZET nfnodes;   // no. face indices
  GSIZET nbnodes;   // no. bdy indices
  GLONG  icurr = 0; // current global index
  GLONG  fcurr = 0; // current global face index
  GLONG  bcurr = 0; // current global bdy index
  for ( auto i=0; i<hmesh_.size(); i++ ) { // for each hex in irank's mesh

    pelem = new GElem_base(this->gtype_, gbasis_);
    xNodes  = &pelem->xNodes();  // node spatial data
    xiNodes = &pelem->xiNodes(); // node ref interval data
    Ni.resize(pelem->nnodes()); // tensor product shape function
    bdy_ind = &pelem->bdy_indices(); // get bdy indices data member
    bdy_typ = &pelem->bdy_types  (); // get bdy types data member
    bdy_ind->clear(); bdy_typ->clear();
    pelem->igbeg() = icurr;      // beginning global index
    pelem->igend() = icurr + pelem->nnodes()-1; // end global index
    for ( auto l=0; l<ndim_; l++ ) { // loop over element Cart coords
      (*xNodes)[l] = 0.0;
      for ( auto m=0; m<pow(2,ndim_); m++ ) { // loop over verts given in Cart coords
        I[0] = m;
        lshapefcn_->Ni(I, *xiNodes, Ni);
        (*xNodes)[l] += Ni * ( (*(hmesh_[i].v[m]))[l] * 0.125 );
      }
    }

    pelem->init(*xNodes);

#if 0
    // With edge/face centroids set, compute global bdy_nodes:
    for ( auto j=0; j<2*ndim_; j++ ) { // cycle over faces
      cent = pelem->faceCentroid(j);
      face_ind = NULLPTR;
      if ( FUZZYEQ(P0_.x2,cent.x2,eps_) ) face_ind = &pelem->edge_indices(0);
      if ( FUZZYEQ(P1_.x1,cent.x1,eps_) ) face_ind = &pelem->edge_indices(1);
      if ( FUZZYEQ(P1_.x2,cent.x2,eps_) ) face_ind = &pelem->edge_indices(2);
      if ( FUZZYEQ(P0_.x1,cent.x1,eps_) ) face_ind = &pelem->edge_indices(3);
      if ( FUZZYEQ(P0_.x3,cent.x3,eps_) ) face_ind = &pelem->edge_indices(4);
      if ( FUZZYEQ(P1_.x3,cent.x3,eps_) ) face_ind = &pelem->edge_indices(5);
//    face_ind = &pelem->face_indices(0);
      // For now, we don't allow corner nodes to be repeated, 
      // and we'll have to choose the best way to define the 
      // normal vectors at the 'corner' nodes:
      for ( auto k=0; face_ind != NULLPTR && k<face_ind->size(); k++ ) {
        if ( !bdy_ind->contains((*face_ind)[k]) ) {
          bdy_ind->push_back((*face_ind)[k]); 
          bdy_typ->push_back(GBDY_NONE); // default always to GBDY_NONE 
        }
      }
    }
#endif

    gelems_.push_back(pelem);

    nvnodes = pelem->nnodes();
    nfnodes = pelem->nfnodes();
    nbnodes = pelem->bdy_indices().size();
    pelem->igbeg() = icurr;           // beg global vol index
    pelem->igend() = icurr+nvnodes-1; // end global vol index
    pelem->ifbeg() = fcurr;           // beg global face index
    pelem->ifend() = fcurr+nfnodes-1; // end global face index
    pelem->ibbeg() = bcurr;           // beg global bdy index
    pelem->ibend() = bcurr+nbnodes-1; // end global bdy index
    icurr += nvnodes;
    fcurr += nfnodes;
    bcurr += nbnodes;

  } // end of hex mesh loop

} // end of method do_elems3d (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : do_elems2d (2)
// DESC   : Build 2d element list from input data.
// ARGS   : p      : matrix of size the number of elements X GDIM containing 
//                   the poly expansion order in each direction
//          gxnodes: vector of GDIM vectors containing Cartesian coords of elements
//                   for each node point
// RETURNS: none.
//**********************************************************************************
void GGridBox::do_elems2d(GTMatrix<GINT> &p,
                          GTVector<GTVector<GFTYPE>> &gxnodes)
{
  assert(gbasis_.size()>0 && "Must set basis pool first");
  assert(ndim_ == 2 && "Dimension must be 2");

  GString                      serr = "GGridBox::do_elems2d (2): ";
  GTPoint<GFTYPE>              cent;
  GTVector<GINT>               I(3);
  GTVector<GINT>              *bdy_ind;
  GTVector<GBdyType>          *bdy_typ;
  GTVector<GINT>              *face_ind;
  GElem_base                  *pelem;
  GTVector<GTVector<GFTYPE>>  *xNodes;
  GTVector<GTVector<GFTYPE>*> *xiNodes;
  GTVector<GTVector<GFTYPE>>   xgtmp(3);
  GTVector<GNBasis<GCTYPE,GFTYPE>*>
                               gb(GDIM);
  GTVector<GINT>               ppool(gbasis_.size());

  assert(gbasis_.size()>0 && "Must set basis first");

  // Now, treat the gbasis_ as a pool that we search
  // to find bases we need:
  for ( auto j=0; j<ppool.size(); j++ ) ppool[j] = gbasis_[j]->getOrder();

  GSIZET iwhere, n;
  GSIZET nvnodes;   // no. vol nodes
  GSIZET nfnodes;   // no. face nodes
  GSIZET nbnodes;   // no. bdy nodes
  GSIZET icurr = 0; // current global index
  GSIZET fcurr = 0; // current global face index
  GSIZET bcurr = 0; // current global bdy index
  for ( auto i=0; i<p.size(1); i++ ) { // for each element
    nvnodes = 1;
    for ( auto j=0; j<GDIM; j++ ) { // set basis from pool
      assert(ppool.contains(p(i,j),iwhere) && "Expansion order not found");
      gb[j] = gbasis_[iwhere];
      nvnodes *= (p(i,j) + 1);
    }
    pelem = new GElem_base(this->gtype_, gb);
    xNodes  = &pelem->xNodes();  // node spatial data
    xiNodes = &pelem->xiNodes(); // node ref interval data
    bdy_ind = &pelem->bdy_indices(); // get bdy indices data member
    bdy_typ = &pelem->bdy_types  (); // get bdy types data member
    bdy_ind->clear(); bdy_typ->clear();

    // Set internal node positions from input data.
    // Note that gxnodes are 'global' and xNodes is
    // element-local:
    for ( auto j=0; j<GDIM; j++ ) {
       gxnodes[j].range(icurr, icurr+nvnodes-1);
      (*xNodes)[j] = gxnodes[j];
    }
    for ( auto j=0; j<GDIM; j++ ) gxnodes[j].range_reset();

    pelem->init(*xNodes);


    // With face/edge centroids computed, compute 
    // global boundary nodes:
    for ( auto j=0; j<2*ndim_; j++ ) { // cycle over all edges
      cent = pelem->edgeCentroid(j);
      face_ind = NULLPTR;
      if ( FUZZYEQ(P0_.x2,cent.x2,eps_) ) face_ind = &pelem->edge_indices(0);
      if ( FUZZYEQ(P1_.x1,cent.x1,eps_) ) face_ind = &pelem->edge_indices(1);
      if ( FUZZYEQ(P1_.x2,cent.x2,eps_) ) face_ind = &pelem->edge_indices(2);
      if ( FUZZYEQ(P0_.x1,cent.x1,eps_) ) face_ind = &pelem->edge_indices(3);
      // For now, we don't allow corner nodes to be repeated, 
      // and we'll have to choose the best way to define the 
      // normal vectors at the 'corner' nodes:
      for ( auto k=0; face_ind != NULLPTR && k<face_ind->size(); k++ ) {
        if ( !bdy_ind->contains((*face_ind)[k]) ) {
          bdy_ind->push_back((*face_ind)[k]); 
          bdy_typ->push_back(GBDY_NONE); 
        }
      }
    }

    gelems_.push_back(pelem);

    // Find global global interior and bdy start & stop indices represented 
    // locally within element:
    assert(nvnodes == gelems_[n]->nnodes() && "Incompatible node count");
    nfnodes = pelem->nfnodes();
    nbnodes = pelem->bdy_indices().size();
    pelem->igbeg() = icurr;           // beg global vol index
    pelem->igend() = icurr+nvnodes-1; // end global vol index
    pelem->ifbeg() = fcurr;           // beg global face index
    pelem->ifend() = fcurr+nfnodes-1; // end global face index
    pelem->ibbeg() = bcurr;           // beg global bdy index
    pelem->ibend() = bcurr+nbnodes-1; // end global bdy index
    icurr += nvnodes;
    fcurr += nfnodes;
    bcurr += nbnodes;
  } // end of quad mesh loop

} // end of method do_elems2d (2)



//**********************************************************************************
//**********************************************************************************
// METHOD : do_elems3d (2)
// DESC   : Build 3d element list from input data. 
// ARGS   : p      : matrix of size the number of elements X GDIM containing 
//                   the poly expansion order in each direction
//          gxnodes: vector of GDIM vectors containing Cartesian coords of elements
//                   for each node point
// RETURNS: none.
//**********************************************************************************
void GGridBox::do_elems3d(GTMatrix<GINT> &p,
                          GTVector<GTVector<GFTYPE>> &gxnodes)
{
  assert(gbasis_.size()>0 && "Must set basis pool first");
  assert(ndim_ == 3 && "Dimension must be 3");

  GString                      serr = "GGridBox::do_elems3d (2): ";
  GTPoint<GFTYPE>              cent;
  GTVector<GINT>               iind;
  GTVector<GINT>               I(3);
  GTVector<GINT>              *bdy_ind;
  GTVector<GBdyType>          *bdy_typ;
  GTVector<GINT>              *face_ind;
  GTVector<GFTYPE>             Ni;
  GElem_base                  *pelem;
  GTVector<GTVector<GFTYPE>>  *xNodes;
  GTVector<GTVector<GFTYPE>*> *xiNodes;
  GTVector<GTVector<GFTYPE>>   xgtmp(3);
  GTVector<GNBasis<GCTYPE,GFTYPE>*>
                               gb(GDIM);
  GTVector<GINT>               ppool(gbasis_.size());

  assert(gbasis_.size()>0 && "Must set basis first");

  // Now, treat the gbasis_ as a pool that we search
  // to find bases we need:
  for ( auto j=0; j<ppool.size(); j++ ) ppool[j] = gbasis_[j]->getOrder();

  GSIZET i, iwhere, n;
  GSIZET nvnodes;   // no. vol indices
  GSIZET nfnodes;   // no. face indices
  GSIZET nbnodes;   // no. bdy indices
  GLONG  icurr = 0; // current global index
  GLONG  fcurr = 0; // current global face index
  GLONG  bcurr = 0; // current global bdy index
  for ( auto i=0; i<p.size(1); i++ ) { // for each hex in irank's mesh
    nvnodes = 1; 
    for ( auto j=0; j<GDIM; j++ ) { // set basis from pool
      assert(ppool.contains(p(i,j),iwhere) && "Expansion order not found");
      gb[j] = gbasis_[iwhere];
      nvnodes *= (p(i,j) + 1);
    }    

    pelem = new GElem_base(this->gtype_, gbasis_);
    xNodes  = &pelem->xNodes();  // node spatial data
    xiNodes = &pelem->xiNodes(); // node ref interval data
    bdy_ind = &pelem->bdy_indices(); // get bdy indices data member
    bdy_typ = &pelem->bdy_types  (); // get bdy types data member
    bdy_ind->clear(); bdy_typ->clear();
    pelem->igbeg() = icurr;      // beginning global index
    pelem->igend() = icurr + pelem->nnodes()-1; // end global index

    // Set internal node positions from input data.
    // Note that gxnodes are 'global' and xNodes is
    // element-local:
    for ( auto j=0; j<GDIM; j++ ) {
       gxnodes[j].range(icurr, icurr+nvnodes-1);
      (*xNodes)[j] = gxnodes[j];
    }
    for ( auto j=0; j<GDIM; j++ ) gxnodes[j].range_reset();

    pelem->init(*xNodes);

#if 0
    // With edge/face centroids set, compute bdy_nodes:
    for ( auto j=0; j<2*ndim_; j++ ) { // cycle over faces
      cent = pelem->faceCentroid(j);
      face_ind = NULLPTR;
      if ( FUZZYEQ(P0_.x2,cent.x2,eps_) ) face_ind = &pelem->edge_indices(0);
      if ( FUZZYEQ(P1_.x1,cent.x1,eps_) ) face_ind = &pelem->edge_indices(1);
      if ( FUZZYEQ(P1_.x2,cent.x2,eps_) ) face_ind = &pelem->edge_indices(2);
      if ( FUZZYEQ(P0_.x1,cent.x1,eps_) ) face_ind = &pelem->edge_indices(3);
      if ( FUZZYEQ(P0_.x3,cent.x3,eps_) ) face_ind = &pelem->edge_indices(4);
      if ( FUZZYEQ(P1_.x3,cent.x3,eps_) ) face_ind = &pelem->edge_indices(5);
      face_ind = &pelem->face_indices(0);

      // For now, we don't allow corner nodes to be repeated, 
      // and we'll have to choose the best way to define the 
      // normal vectors at the 'corner' nodes:
      for ( auto k=0; face_ind != NULLPTR && k<face_ind->size(); k++ ) {
        if ( !bdy_ind->contains((*face_ind)[k]) ) {
          bdy_ind->push_back((*face_ind)[k]); 
          bdy_typ->push_back(GBDY_NONE); // default always to GBDY_NONE 
        }
      }
    }
#endif

    gelems_.push_back(pelem);

    assert(nvnodes == gelems_[i]->nnodes() && "Incompatible node count");
    nfnodes = pelem->nfnodes();
    nbnodes = pelem->bdy_indices().size();
    pelem->igbeg() = icurr;           // beg global vol index
    pelem->igend() = icurr+nvnodes-1; // end global vol index
    pelem->ifbeg() = fcurr;           // beg global face index
    pelem->ifend() = fcurr+nfnodes-1; // end global face index
    pelem->ibbeg() = bcurr;           // beg global bdy index
    pelem->ibend() = bcurr+nbnodes-1; // end global bdy index
    icurr += nvnodes;
    fcurr += nfnodes;
    bcurr += nbnodes;
  } // end of hex mesh loop

} // end of method do_elems3d (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : print
// DESC   : Print final mesh points. It is often the case 
//          that we will print the mesh to a file for visualization, so this is a 
//          utility that allows us to do this easily. A stream operator is 
//          still provided to print in a completely formatted way.
// ARGS   : filename: filename to print to
// RETURNS: none.
//**********************************************************************************
void GGridBox::print(const GString &filename)
{
  GString serr = "GGridBox::print: ";
  std::ofstream ios;

  GFTYPE x, y, z;
  GTPoint<GFTYPE> pt;

  ios.open(filename);
  if ( ndim_ == 2 ) {
    for( auto i=0; i<qmesh_.size(); i++ ) { // for each quad/hex
      for( auto j=0; j<pow(2,ndim_); j++ ) { // for each vertex of reg polygon
          pt = *qmesh_[i].v[j];
          ios << pt.x1 << " " <<  pt.x2 << std::endl;
          ios << pt.x1 << " " << pt.x2 << " " <<  pt.x3 << std::endl;
      }
    }
  } 

  if ( ndim_ == 3 ) {
    for( auto i=0; i<hmesh_.size(); i++ ) { // for each quad/hex
      for( auto j=0; j<pow(2,ndim_); j++ ) { // for each vertex of reg polygon
          pt = *hmesh_[i].v[j];
          ios << pt.x1 << " " <<  pt.x2 << std::endl;
          ios << pt.x1 << " " << pt.x2 << " " <<  pt.x3 << std::endl;
      }
    }
  }
  ios.close();

} // end of method print


//**********************************************************************************
//**********************************************************************************
// METHOD : periodize
// DESC   : Periodize grid coordinates, by making x_2 = x_1 for all
//          periodic coordinate directions, x.
// ARGS   : none
// RETURNS: none.
//**********************************************************************************
void GGridBox::periodize()
{
  assert(bInitialized_ && "Object not initialized");

  
  GTVector<GFTYPE>  x(xNodes_.size()); // coord values to set to


  periodicids_.clear();
  periodicdirs_.clear();

  // Coords set to correspond to bottom-most domain point:
  for( auto k=0; k<x.size(); k++ ) x[k] = P0_[k];

  GUINT  bit;
  GSIZET id, n, num=0;
  for ( auto k=0; k<igbdy_binned_.size(); k++ ) {
    num += igbdy_binned_[k][GBDY_PERIODIC].size();
  }
  periodicids_ .resize(num);
  periodicdirs_.resize(num);

  n = 0;
  for( auto k=0; k<igbdy_binned_.size(); k++ ) { // for each global face
    for( auto j=0; j<igbdy_binned_[k][GBDY_PERIODIC].size(); j++, n++ ) { // for each global bdy node
      id = igbdy_binned_[k][GBDY_PERIODIC][j];
      periodicids_ [n] = id;       
      periodicdirs_[n] = 0;
      for( auto i=0; i<xNodes_.size(); i++ ) { // for x, y, z dirs
        if ( FUZZYEQ(P1_[i],xNodes_[i][id],eps_) ) { // right/top-mosstblock.tbdy[k];i coord will change
          periodicdirs_[n] |= 1U << i;  // position right-most direction bit  
        }
      }
    }
  }

  // Now, cycle through periodic nodes and periodize coordinates:
  for( auto k=0; k<periodicids_.size(); k++ ) { // for each periodic node
    id = periodicids_[k];
    for( auto i= 0; i<xNodes_.size(); i++ ) { // coord direction
      // Set coord in this direction if corresp bit is set:
      bit = (periodicdirs_[k] >> i) & 1; 
      if ( bit ) xNodes_[i][id] = x[i];
    }
  }

} // end of method periodize


//**********************************************************************************
//**********************************************************************************
// METHOD : unperiodize
// DESC   : Un-periodize grid, if of appropriate type
// ARGS   : none
// RETURNS: none.
//**********************************************************************************
void GGridBox::unperiodize()
{
  // Use data from 'periodize' method to unset change in
  // coordinates:

  GTVector<GFTYPE>  x(xNodes_.size()); // coord values to set to

  // Coords to set to correspond to top-most domain point:
  for( auto k=0; k<x.size(); k++ ) x[k] = P1_[k];

  // Cycle through periodic nodes and un-periodize coordinates:
  GUINT  bit;
  GSIZET id;
  for( auto k=0; k<periodicids_.size(); k++ ) { // for each periodic node
    id = periodicids_[k];
    for( auto i= 0; i<xNodes_.size(); i++ ) { // coord direction
      // Set coord in this direction if corresp bit is set:
      bit = (periodicdirs_[k] >> i) & 1; 
      if ( bit ) xNodes_[i][id] = x[i];
    }
  }

  periodicids_.clear();
  periodicdirs_.clear();

} // end of method unperiodize


//**********************************************************************************
//**********************************************************************************
// METHOD : find_subdomain
// DESC   : Find this rank's default portion of global domain, and
//          store in qmesh member data.
// ARGS   : none.
// RETURNS: none.
//**********************************************************************************
void GGridBox::find_subdomain()
{

  GSIZET          n=0, nglobal, nperrank, nthisrank, ntot, nxy;
  GLONG           beg_lin, end_lin;
  GLONG           ib, ie, jb, je, kb, ke;
  GTPoint<GFTYPE> v0(ndim_);     // starting sph. point of elem
  GTPoint<GFTYPE> dv(ndim_);     // delta-point
  GTPoint<GFTYPE> dx(ndim_);

  nglobal = 1; // total num elements in global grid
  for( auto k=0; k<ne_.size(); k++ ) nglobal *= ne_[k];
 
  nperrank = nglobal / nprocs_; // #elems per rank
  nthisrank = irank_ != nprocs_-1 ? nperrank : nglobal - (nprocs_-1)*nperrank;


 // Get uniform element sizes:
  for( auto k=0; k<ndim_; k++ ) {
    dx[k] = Lbox_[k] / static_cast<GFTYPE>(ne_[k]);
  }

  // Compute vertices of all hexes ('cubes')
  // for this task:
  if ( ndim_ == 2 ) {

    qmesh_.resize(nthisrank);
    beg_lin = nperrank*irank_; end_lin = beg_lin + nthisrank - 1;
    jb = beg_lin/ne_[0]; je = end_lin/ne_[0];
    for ( GLONG j=jb; j<=je; j++ ) {
      ib = MAX(beg_lin-static_cast<GLONG>(j*ne_[0]),0); 
      ie = MIN(end_lin-static_cast<GLONG>(j*ne_[0]),ne_[0]-1); 
      for ( GLONG i=ib; i<=ie; i++ ) {
        for ( auto l=0; l<4; l++ ) qmesh_[n][l].resize(ndim_);
        v0.x1 = P0_.x1+i*dx.x1; v0.x2 = P0_.x2+j*dx.x2;
                                         qmesh_[n].v1 = v0;
        dv.x1 = dx.x1 ; dv.x2 = 0.0  ;   qmesh_[n].v2 = v0 + dv;
        dv.x1 = dx.x1 ; dv.x2 = dx.x2;   qmesh_[n].v3 = v0 + dv;
        dv.x1 = 0.0   ; dv.x2 = dx.x2;   qmesh_[n].v4 = v0 + dv;
        n++;
      }
    }

  } // end, ndim==2 test
  else if ( ndim_ == 3 ) {
    hmesh_.resize(nthisrank);
    nxy = ne_[0] * ne_[1];
    beg_lin = nperrank*irank_; end_lin = beg_lin + nthisrank - 1;
    kb = beg_lin/nxy; ke = end_lin/nxy;
    for ( GLONG k=kb; k<=ke; k++ ) { 
      jb = MAX((beg_lin-static_cast<GLONG>(k*nxy))/ne_[0],0); 
      je = MIN((end_lin-static_cast<GLONG>(k*nxy))/ne_[0],ne_[1]-1);
      for ( GLONG j=jb; j<=je; j++ ) { 
        ib = MAX(beg_lin-static_cast<GLONG>(k*nxy+j*ne_[0]),0); 
        ie = MIN(end_lin-static_cast<GLONG>(k*nxy+j*ne_[0]),ne_[0]-1); 
        for ( GLONG i=ib; i<=ie; i++ ) { 
          v0.x1 = P0_.x1+i*dx.x1; v0.x2 = P0_.x2+j*dx.x2; v0.x3 = P0_.x3+k*dx.x3; 
                                                          hmesh_[n].v1 = v0;
          dv.x1 = dx.x1 ; dv.x2 = 0.0  ; dv.x3 = 0.0   ;  hmesh_[n].v2 = v0 + dv;
          dv.x1 = dx.x1 ; dv.x2 = dx.x2; dv.x3 = 0.0   ;  hmesh_[n].v3 = v0 + dv;
          dv.x1 = 0.0   ; dv.x2 = dx.x2; dv.x3 = 0.0   ;  hmesh_[n].v4 = v0 + dv;
          dv.x1 = 0.0   ; dv.x2 = 0.0  ; dv.x3 = dx.x3 ;  hmesh_[n].v5 = v0 + dv;
          dv.x1 = dx.x1 ; dv.x2 = 0.0  ; dv.x3 = dx.x3 ;  hmesh_[n].v6 = v0 + dv;
          dv.x1 = dx.x1 ; dv.x2 = dx.x2; dv.x3 = dx.x3 ;  hmesh_[n].v7 = v0 + dv;
          dv.x1 = 0.0   ; dv.x2 = dx.x2; dv.x3 = dx.x3 ;  hmesh_[n].v8 = v0 + dv;
          n++;
        }
      }
    }

  } // end, ndim==3 test


} // end of method find_subdomain


//**********************************************************************************
//**********************************************************************************
// METHOD : config_bdy
// DESC   : Configure 2d & 3d box boundary from ptree
// ARGS   : 
//          ptree : main prop tree 
//          igbdy : For each natural/canonical global boundary face,
//                  gives vector of global bdy ids. Allocated here.
//          igbdyt: bdy type ids for each index in igbdy. Allocated here.
// RETURNS: none.
//**********************************************************************************
void GGridBox::config_bdy(const PropertyTree &ptree, 
                          GTVector<GTVector<GSIZET>> &igbdy, 
                          GTVector<GTVector<GBdyType>> &igbdyt)
{
  // Cycle over all geometric boundaries, and configure:

  GBOOL              bret, bperiodic=FALSE;
  GSIZET             iwhere;
  GTVector<GBOOL>    buniform(2*GDIM);
  GTVector<GBdyType> bdytype(2*GDIM);
  GTVector<GSIZET>   itmp;
  GTVector<GString>  bdynames (2*GDIM);
  std::vector<GString>
                     svec;
  GString            gname, sbdy, bdyclass;
  PropertyTree       bdytree, gridptree;
  stBdyBlock         stblock;
  UpdateBasePtr      base_ptr;

  bdynames[0] = "bdy_x_0";
  bdynames[1] = "bdy_x_1";
  bdynames[2] = "bdy_y_0";
  bdynames[3] = "bdy_y_1";
  if ( GDIM == 3 ) {
    bdynames[4] = "bdy_z_0";
    bdynames[5] = "bdy_z_1";
  }

  gname     = ptree.getValue<GString>("grid_type");
  assert(gname == "grid_box");
  gridptree = ptree.getPropertyTree(gname);


  // Clear input arrays:
  igbdy .resize(2*GDIM);
  igbdyt.resize(2*GDIM);

  bdy_update_list_.resize(2*GDIM);

  // Handle uniform, nonuniform bdy conditions:
  // Note: If "uniform" not specified for a boundary, then
  //       user MUST supply a method to configure it.
  for ( auto j=0; j<2*GDIM; j++ ) { 
    sbdy         = gridptree.getValue<GString>(bdynames[j]);
cout << "config_bdy: getting bdy tree: " << sbdy << endl;
    bdytree      = ptree.getPropertyTree(sbdy);
    bdyclass     = bdytree.getValue<GString>("bdy_class", "uniform");
    if ( ndim_ == 2 ) 
      find_bdy_ind2d(j, TRUE, itmp);
    if ( ndim_ == 3 ) 
      find_bdy_ind3d(j, TRUE, itmp);
    igbdy [j].resize(itmp.size()); igbdy [j] = itmp;
    igbdyt[j].resize(itmp.size()); igbdyt[j] = GBDY_NONE;
    if ( "uniform" == bdyclass ) { // uniform bdy conditions
cout << "config_bdy: extracting data from bdy tree: " << sbdy << endl;
      geoflow::get_bdy_block(bdytree, stblock);
      if ( stblock.tbdy.contains(GBDY_PERIODIC) ) {
        assert(stblock.tbdy.onlycontains(GBDY_PERIODIC) && "All variables must be GBDY_PERIODIC");
        bdytype  [j] = GBDY_PERIODIC;
        igbdyt   [j] = GBDY_PERIODIC;
        bperiodic    = bperiodic || bdytype[j] == GBDY_PERIODIC;
      
      }
      // May have different uniform bdys for different state comps:
      for ( auto k=0; k<stblock.tbdy.size() && !bperiodic; k++ ) {
cout << "config_bdy: building bc for bdy cond " << k << endl;
        base_ptr = GUpdateBdyFactory<BdyTypePack>::build(ptree, sbdy, *this,  j,
                                            stblock.tbdy[k], stblock.istate[k], itmp);
        igbdyt[j] = stblock.tbdy[k];
        bdy_update_list_[j].push_back(base_ptr);
      }
    }
    else if ( "mixed" == bdyclass ) { // mixed bdy conditions
      assert( bdytree.isArray<GString>("bdy_blocks") && "no bdy_blocks specified");
      svec = bdytree.getArray<GString>("bdy_blocks");
      for ( auto i=0; i<svec.size(); i++ ) { // loop over bdy blocks
        assert( ptree.isPropertyTree(svec[i]) && "no component bdy_blocks specified");
        bdytree = ptree.getPropertyTree(svec[i]);
        geoflow::get_bdy_block(bdytree, stblock);
        assert(!stblock.tbdy.contains(GBDY_PERIODIC) && "GBDY_PERIODIC bdys must be uniform");
        GSpecBdyFactory::dospec(bdytree, *this, j, itmp);
        for ( auto k=0; k<svec.size(); k++ ) { // for each sub-block
          base_ptr = GUpdateBdyFactory<BdyTypePack>::build(ptree, svec[k], *this,  j,
                                              stblock.tbdy[k], stblock.istate[k], itmp);

          for ( auto m=0; m<itmp.size(); m++ ) {
            if ( igbdy[j].contains(itmp[m]) ) igbdyt[j][m] = stblock.tbdy[k];
          }
          if ( stblock.tbdy[k] != GBDY_NONE ) igbdyt[j] = stblock.tbdy[k];
          bdy_update_list_[j].push_back(base_ptr);
        }
      }
    }
    else {
      assert(FALSE && "Invalid bdy_class");
    }


  } // end, global bdy face loop

  if ( bperiodic ) {
    if ( ndim_ == 2 ) {
      assert( (  (bdytype[0] == GBDY_PERIODIC && bdytype[2] == GBDY_PERIODIC)
             ||  (bdytype[3] == GBDY_PERIODIC && bdytype[1] == GBDY_PERIODIC) )
             &&  "Incompatible GBDY_PERIODIC boundary specification");
    }
    else if ( ndim_ == 3 ) {
      assert( (  (bdytype[0] == GBDY_PERIODIC && bdytype[2] == GBDY_PERIODIC)
             ||  (bdytype[3] == GBDY_PERIODIC && bdytype[1] == GBDY_PERIODIC)  
             ||  (bdytype[4] == GBDY_PERIODIC && bdytype[5] == GBDY_PERIODIC) )
             &&  "Incompatible GBDY_PERIODIC boundary specification");
    }
  }


} // end of method config_bdy


//**********************************************************************************
//**********************************************************************************
// METHOD : find_bdy_ind2d
// DESC   : Find global bdy indices (indices into xNodes_ arrays) that
//          corresp to specified bdy in 2d
// ARGS   : bdyid    : box-boundary id
//          incl_vert: include vertex points on bdy that are shared 
//                     with another bdy
//          ibdy     : array of indices into xNodes that comprise this boundary
// RETURNS: none.
//**********************************************************************************
void GGridBox::find_bdy_ind2d(GINT bdyid, GBOOL incl_vert, GTVector<GSIZET> &ibdy)
{
  assert(bdyid >=0 && bdyid < 2*GDIM);

  GTPoint<GFTYPE> pt(ndim_);
  
  ibdy.clear();
  pt.setBracket(eps_);

  switch ( bdyid ) {
    case 0: // lower horiz bdy:
      for( auto i=0; i<xNodes_[0].size(); i++ ) { 
        if ( FUZZYEQ(P0_.x2,xNodes_[1][i],eps_) ) {
          pt.assign(xNodes_, i);
          if ( incl_vert  || !is_global_vertex(pt) ) ibdy.push_back(i);
        }
      }
      break;

    case 1: // right vert bdy:
      for( auto i=0; i<xNodes_[0].size(); i++ ) { // bdy 1
        if ( FUZZYEQ(P1_.x1,xNodes_[0][i],eps_) ) {
          pt.assign(xNodes_, i);
          if ( incl_vert  || !is_global_vertex(pt) ) ibdy.push_back(i);
        }
      }
      break;

    case 2: // top horiz bdy:
      for( auto i=0; i<xNodes_[0].size(); i++ ) { // bdy 2
        if ( FUZZYEQ(P1_.x2,xNodes_[1][i],eps_) ) {
          pt.assign(xNodes_, i);
          if ( incl_vert  || !!is_global_vertex(pt) ) ibdy.push_back(i);
        }
      }
      break;

    case 3: // left vert bdy:
      for( auto i=0; i<xNodes_[0].size(); i++ ) { // bdy 3
        if ( FUZZYEQ(P0_.x1,xNodes_[0][i],eps_) ) {
          pt.assign(xNodes_, i);
          if ( incl_vert  || !is_global_vertex(pt) ) ibdy.push_back(i);
        }
      }
      break;

    default : // error:
      assert(FALSE && "Invalid 2D global bdy id");
    
  } // end, switch 

} // end, method find_bdy_ind2d


//**********************************************************************************
//**********************************************************************************
// METHOD : find_bdy_ind3d
// DESC   : Find global bdy indices (indices into xNodes_ arrays) that
//          corresp to specified bdy in 3d
// ARGS   : bdyid    : box-boundary id
//          incl_edge: include edge points on bdy that are shared 
//                     with another bdy
//          ibdy     : array of indices into xNodes that comprise this boundary
// RETURNS: none.
//**********************************************************************************
void GGridBox::find_bdy_ind3d(GINT bdyid, GBOOL incl_edge, GTVector<GSIZET> &ibdy)
{

  assert(bdyid >=0 && bdyid < 2*GDIM);

  GSIZET nbdy, n;
  GTPoint<GFTYPE> pt(ndim_);

  ibdy.clear();
  nbdy = n = 0;
  switch ( bdyid ) {
    case 0: // southern vert face:
      for( auto i=0; i<xNodes_[0].size(); i++ ) { // face 0
        if ( FUZZYEQ(P0_.x2,xNodes_[1][i],eps_) ) {
          pt.assign(xNodes_, i);
          if ( incl_edge || !on_global_edge(0,pt) ) nbdy++;
        }
      }
      ibdy.resize(nbdy);
      for( auto i=0; i<xNodes_[0].size(); i++ ) { // face 0
        if ( FUZZYEQ(P0_.x2,xNodes_[1][i],eps_) ) {
          pt.assign(xNodes_, i);
          if ( incl_edge || !on_global_edge(0,pt) ) ibdy[n++] = i;
        }
      }
      break;

    case 1: // eastern vert. face:
      for( auto i=0; i<xNodes_[0].size(); i++ ) { // face 1
        if ( FUZZYEQ(P1_.x1,xNodes_[0][i],eps_) ) {
          pt.assign(xNodes_, i);
          if ( incl_edge || !on_global_edge(1,pt) ) nbdy++;
        }
      }
      ibdy.resize(nbdy);
      for( auto i=0; i<xNodes_[0].size(); i++ ) { // face 1
        if ( FUZZYEQ(P1_.x1,xNodes_[0][i],eps_) ) {
          pt.assign(xNodes_, i);
          if ( incl_edge || !on_global_edge(1,pt) ) ibdy[n++] = i;
        }
      }
      break;

    case 2: // northern vert. face:
      for( auto i=0; i<xNodes_[0].size(); i++ ) { // face 2
        if ( FUZZYEQ(P1_.x2,xNodes_[1][i],eps_) ) {
          pt.assign(xNodes_, i);
          if ( incl_edge || !on_global_edge(2,pt) ) nbdy++;
        }
      }
      ibdy.resize(nbdy);
      for( auto i=0; i<xNodes_[0].size(); i++ ) { // face 2
        if ( FUZZYEQ(P1_.x2,xNodes_[1][i],eps_) ) {
          pt.assign(xNodes_, i);
          if ( incl_edge || !on_global_edge(2,pt) ) ibdy[n++] = i;
        }
      }
      break;

    case 3: // western vertical face:
      for( auto i=0; i<xNodes_[0].size(); i++ ) { // face 3
        if ( FUZZYEQ(P0_.x1,xNodes_[0][i],eps_) ) {
          pt.assign(xNodes_, i);
          if ( incl_edge || !on_global_edge(3,pt) ) nbdy++;
        }
      }
      ibdy.resize(nbdy);
      for( auto i=0; i<xNodes_[0].size(); i++ ) { // face 3
        if ( FUZZYEQ(P0_.x1,xNodes_[0][i],eps_) ) {
          pt.assign(xNodes_, i);
          if ( incl_edge || !on_global_edge(3,pt) ) ibdy[n++] = i;
        }
      }
      break;

    case 4: // bottom horiz face:
      for( auto i=0; i<xNodes_[0].size(); i++ ) { // face 4
        if ( FUZZYEQ(P0_.x3,xNodes_[2][i],eps_) ) {
          pt.assign(xNodes_, i);
          if ( incl_edge || !on_global_edge(4,pt) ) nbdy++;
        }
      }
      ibdy.resize(nbdy);
      for( auto i=0; i<xNodes_[0].size(); i++ ) { // face 4
        if ( FUZZYEQ(P0_.x3,xNodes_[2][i],eps_) ) {
          pt.assign(xNodes_, i);
          if ( incl_edge || !on_global_edge(4,pt) ) ibdy[n++] = i;
        }
      }
      break;

    case 5: // top horiz face:
      for( auto i=0; i<xNodes_[0].size(); i++ ) { // face 5
        if ( FUZZYEQ(P1_.x3,xNodes_[2][i],eps_) ) {
          pt.assign(xNodes_, i);
          if ( incl_edge || !on_global_edge(5,pt) ) nbdy++;
        }
      }
      ibdy.resize(nbdy);
      for( auto i=0; i<xNodes_[0].size(); i++ ) { // face 5
        if ( FUZZYEQ(P1_.x3,xNodes_[2][i],eps_) ) {
          pt.assign(xNodes_, i);
          if ( incl_edge || !on_global_edge(5,pt) ) ibdy[n++] = i;
        }
      }
      break;

    default : // error:
      assert(FALSE && "Invalid 3D globa bdy id");

  } // end, switch

} // end, method find_bdy_ind3d


//**********************************************************************************
//**********************************************************************************
// METHOD : is_global_vertex
// DESC   : Utilitiy method to determine of specified point is one of the 
//          global 2d boundary vertices
// ARGS   : pt : point to check
// RETURNS: TRUE on yes, else FALSE
//**********************************************************************************
GBOOL GGridBox::is_global_vertex(GTPoint<GFTYPE> &pt)
{
  GBOOL           bret = FALSE;

  for( auto j=0; j<pow(2,ndim_) && !bret; j++ ) {
    bret = bret || ( pt == gverts_[j] ); // There is fuzziness in ==
  }

  return bret;

} // end, method is_global_vertex


//**********************************************************************************
//**********************************************************************************
// METHOD : on_global_edge
// DESC   : Utilitiy method to determine of specified point is on 
//          edge of global 3d boundary
// ARGS   : iface: face index to check
//          pt   : point to check
// RETURNS: TRUE on yes, else FALSE
//**********************************************************************************
GBOOL GGridBox::on_global_edge(GINT iface, GTPoint<GFTYPE> &pt)
{

  assert( iface >=0 && iface <=5 && "Invalid face ID specification");

  GBOOL           bret = FALSE;
  GINT            nface=0;
  GTVector<GINT>  face(3); // at most 3 faces that point _can_ belong to
//GTPoint<GFTYPE> pt(ndim_);

  // Find faces point belongs to:
  for ( GINT j=0; j<ndim_; j++ ) {
    if     ( FUZZYEQ(pt.x1,P0_.x1,eps_) && !face.containsn(0,nface) ) 
      { face[nface] = 0; nface++; }
    else if( FUZZYEQ(pt.x2,P1_.x2,eps_) && !face.containsn(1,nface) ) 
      { face[nface] = 1; nface++; } 
    else if( FUZZYEQ(pt.x1,P1_.x1,eps_) && !face.containsn(2,nface) ) 
      { face[nface] = 2; nface++; }
    else if( FUZZYEQ(pt.x2,P0_.x2,eps_) && !face.containsn(3,nface) ) 
      { face[nface] = 3; nface++; }
    else if( FUZZYEQ(pt.x3,P0_.x3,eps_) && !face.containsn(4,nface) ) 
      { face[nface] = 4; nface++; }
    else if( FUZZYEQ(pt.x3,P1_.x3,eps_) && !face.containsn(5,nface) ) 
      { face[nface] = 5; nface++; }
  }

  if ( nface == 0 ) return FALSE; // in volume somewhere

  if ( nface == 1 ) return FALSE; // on some face, not on an edge or vertex

  GINT iedges[][4][2] = { // for each face, faces comprising each edge 
                         { {0,4},{0,1},{0,5},{2,3} },
                         { {1,4},{1,2},{1,5},{0,1} },
                         { {2,4},{2,3},{2,5},{1,2} },
                         { {3,4},{3,0},{3,5},{2,3} },
                         { {0,4},{1,4},{2,4},{3,4} },
                         { {0,5},{1,5},{2,5},{3,5} },
                       };
   
  // Find which edge, if any, edge-point sits in:
  for ( GINT j=0; j<4 && !bret; j++ ) {
    bret = bret || ( face.contains(iedges[iface][j][0]) && 
                     face.contains(iedges[iface][j][1]) );
  }
  
  return bret;
} // end, method on_global_edge


//**********************************************************************************
//**********************************************************************************
// METHOD : do_face_normals
// DESC   : Compute normals to each element face
// ARGS   : none 
// RETURNS: none
//**********************************************************************************
void GGridBox::do_face_normals()
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
void GGridBox::do_face_normals2d()
{

#if 0
  // Cycle through local elem face indices to set
  // normals. Taken in order, these should correspond
   GSIZET m=0, nn=0;
   GSIZET ibeg, iend;   // beg, end indices for global arrays
   GSIZET ibbeg, ibend; // beg, end indices for global arrays for bdy quantities
   GSIZET ifbeg, ifend; // beg, end indices for global arrays for face quantities
   GTVector<GTVector<GINT>>   *ieface ; // domain face indices
   GTVector<GSIZET>            gieface; // global element face indices
   GTVector<GINT>             *iverts ; // elem vertex indices
   gieface.resize(gieface_.size());
   for( auto e=0; e<gelems_.size(); e++ ) {
     ibeg   = gelems_[e]->igbeg(); iend  = gelems_[e]->igend();
     ifbeg  = gelems_[e]->ifbeg(); ifend = gelems_[e]->ifend();
     ibbeg  = gelems_[e]->ibbeg(); ibend = gelems_[e]->ibend();
     ieface = &gelems_[e]->face_indices();
  

     for( auto j=0; j<ieface->size(); j++ ) { // cycle over all elem faces
       iverts = &gelems_[e]->vert_indices(j);
       for( auto k=0; k<(*ieface)[j].size(); k++ ) {
         ig = nn + (*ieface)[j][k];
         if ( !gieface.containsn(ig, m) ) { // don't include repeated face ind
           gieface[m] = ig;
           m++;
           if      ( j == 0 ) {
             if ( iverts->contains((*ieface)[j][k]:w

           }
           else if ( j == 1 ) {
           }
           else if ( j == 1 ) {
           }
           else if ( j == 1 ) {
           }
         }
       }
     }
     nn += gelems_[e]->nnodes();


   } // end, element loop
#endif


} // end, method do_face_normals2d


//**********************************************************************************
//**********************************************************************************
// METHOD : do_face_normals3d
// DESC   : Compute normals to each element face in 3d
// ARGS   : none 
// RETURNS: none
//**********************************************************************************
void GGridBox::do_face_normals3d()
{

  assert(FALSE);

} // end, method do_face_normals3d


//**********************************************************************************
//**********************************************************************************
// METHOD : do_bdy_normals
// DESC   : Compute normals to each domain bdy 
// ARGS   : 
//          dXdXi     : matrix of dX_i/dXi_j matrix elements, s.t.
//                      dXdX_i(i,j) = dx^j/dxi^i
//          igbdy_face: vector of bdy indices into global volume fields 
//                      for each face
//          normals   : vector of normal components
//          idepComp  : vector index dependent on the other indices (first 
//                      component index whose normal component is nonzero)
// RETURNS: none
//**********************************************************************************
void GGridBox::do_bdy_normals(GTMatrix<GTVector<GFTYPE>>    &dXdXi,
                              GTVector<GTVector<GSIZET>>    &igbdy_face,
                              GTVector<GTVector<GFTYPE>>    &normals,
                              GTVector<GINT>                &idepComp)
{

  GSIZET icurr, nbdy, nface;

  nbdy = 0;
  for ( auto j=0; j<igbdy_face.size(); j++ ) {
    nbdy += igbdy_face[j].size();
  }
  idepComp.resize(nbdy);
  for ( auto j=0; j<normals.size(); j++ ) normals[j].resize(nbdy);


  // Compute global boundary normals and associated data:

  #if defined(_G_IS2D)

  icurr = 0;
  for ( auto j=0; j<2*GDIM; j++ ) { // for each global bdy face 
    nface = igbdy_face[j].size();   // # bdy nodes on this face
    idepComp.range(icurr,icurr+nface-1);
    for ( auto i=0; i<normals.size(); i++ ) normals[i].range(icurr,icurr+nface-1);
    do_bdy_normals2d(dXdXi, igbdy_face[j], j, normals, idepComp);
    icurr += nface;
  }

  #elif defined(_G_IS3D)

  icurr = 0;
  for ( auto j=0; j<2*GDIM; j++ ) { // for each global bdy face 
    nface = igbdy_face[j].size();   // # bdy nodes on this face
    idepComp.range(icurr,icurr+nface-1);
    for ( auto i=0; i<normals.size(); i++ ) normals[i].range(icurr,icurr+nface-1);
    do_bdy_normals3d(dXdXi, igbdy_face[j], j, normals, idepComp);
    icurr += nface;
  }

  #else
    #error Invalid problem dimensionality
  #endif

  // Reset vector ranges:
  idepComp.range_reset();
  for ( auto j=0; j<normals.size(); j++ ) normals[j].range_reset();


} // end, method do_bdy_normals


//**********************************************************************************
//**********************************************************************************
// METHOD : do_bdy_normals2d
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
//          dXdXi   : matrix of dX_i/dXi_j matrix elements, s.t.
//                    dXdX_i(i,j) = dx^j/dxi^i
//          igbdy   : vector of bdy indices into global volume fields
//          iface   : which global face igbdy list represents
//          normals : vector of normal components
//          idepComp: vector index dependent on the other indices (first 
//                    component index whose normal component is nonzero)
//          Note:
//          dXdXi_  : matrix of dX_i/dXi_j matrix elements, s.t.
//                    dXdX_i(i,j) = dx^j/dxi^i,
//                    must be computed prior to entry.
// RETURNS: none
//**********************************************************************************
void GGridBox::do_bdy_normals2d(GTMatrix<GTVector<GFTYPE>>    &dXdXi,
                                GTVector<GSIZET>              &igbdy,
                                GINT                           iedge,
                                GTVector<GTVector<GFTYPE>>    &normals,
                                GTVector<GINT>               &idepComp)
{
   GSIZET          ib, ic, ip; 
   GFTYPE          tiny;
   GFTYPE          xm;
   GTPoint<GFTYPE> kp(3), xp(3), p1(3), p2(3);

   tiny  = 100.0*std::numeric_limits<GFTYPE>::epsilon(); 
   kp    = 0.0;
   kp[2] = 1.0; // k-vector

   // Normals depend on element type:
   if ( this->gtype_ == GE_REGULAR ) {
     // All normal components are 0, except the one
     // perp. to iedge:
     ip = (iedge+1)%2; // perp component
     for ( auto j=0; j<igbdy.size(); j++ ) { // all points on iedge
       xm = iedge == 1 || iedge == 2 ? -1.0 : 1.0;
       for ( auto i=0; i<normals.size(); i++ ) normals[i][j] = 0.0; 
       normals[ip][j] = xm;
       idepComp[j]    = ip; // dependent component
     }
   }
   else if ( this->gtype_ == GE_DEFORMED ) {
     // Bdy normal is hat{k} X dvec{X} / dxi_iedge,
     // for edge iedge:
     for ( auto j=0; j<igbdy.size(); j++ ) { // all points on iedge
       ib = igbdy[j];
       xm = iedge == 1 || iedge == 2 ? -1.0 : 1.0;
       for ( auto i=0; i<dXdXi.size(2); i++ ) { // over _X_
         p1[i] = dXdXi(iedge%2,i)[ib]; 
       }
       kp.cross(p1, xp);   // xp = k X p1
       xp.unit();
       for ( ic=0; ic<GDIM; ic++ ) if ( fabs(xp[ic]) > tiny ) break;
       assert(ic >= GDIM); // no normal components > 0
       for ( auto i=0; i<normals.size(); i++ ) normals[i][j] = xp[i];
       idepComp[j] = ic;  // dependent component
     }
   }
   else if ( this->gtype_ == GE_2DEMBEDDED ) {
     // Bdy normal is dvec{X} / dxi_xi X dvec{X} / dxi_eta
     for ( auto j=0; j<igbdy.size(); j++ ) { // all points on iedge
       ib = igbdy[j];
       for ( auto i=0; i<dXdXi.size(2); i++ ) { // d_X_/dXi
         p1[i] = dXdXi(0,i)[ib]; // d_X_/dxi
         p2[i] = dXdXi(1,i)[ib]; // d_X_/deta
       }
       p1.cross(p2, xp);   // xp = p1 X p2
       xp.unit();
       for ( ic=0; ic<xp.dim(); ic++ ) if ( fabs(xp[ic]) > tiny ) break;
       assert(ic >= GDIM); // no normal components > 0
       for ( auto i=0; i<normals.size(); i++ ) normals[i][j] = xp[i];
       idepComp[j] = ic;  // dependent component
     }
   }

} // end, method do_bdy_normals2d


//**********************************************************************************
//**********************************************************************************
// METHOD : do_bdy_normals3d
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
//          dXdXi   : matrix of dX_i/dXi_j matrix elements, s.t.
//                    dXdX_i(i,j) = dx^j/dxi^i
//          igbdy   : vector of bdy indices into global volume fields
//          iface   : which global face igbdy list represents
//          normals : vector of normal components
//          idepComp: vector index dependent on the other indices (first 
//                    component index whose normal component is nonzero)
// RETURNS: none
//**********************************************************************************
void GGridBox::do_bdy_normals3d(GTMatrix<GTVector<GFTYPE>>    &dXdXi,
                                GTVector<GSIZET>              &igbdy,
                                GINT                           iface,
                                GTVector<GTVector<GFTYPE>>   &normals,
                                GTVector<GINT>             &idepComp)
{
   GSIZET          ib, ic, ip; 
   GINT            ixi[6][2] = { {0,2}, {1,2}, {0,2}, 
                               {1,2}, {0,1}, {0,1} };
   GFTYPE          xsgn  [] = { 1.0, -1.0, -1.0, 1.0, 1.0, -1.0};
   GFTYPE          tiny;
   GFTYPE          xm;
   GTPoint<GFTYPE> xp(3), p1(3), p2(3);
   tiny  = 100.0*std::numeric_limits<GFTYPE>::epsilon(); 

   // Normals depend on element type:
   if ( this->gtype_ == GE_REGULAR ) {
     // All normal components are 0, except the one
     // perp. to face:
     ip = iface < 4 ? (iface+1)%2 : 2; // perp component
     xm = iface == 1 || iface == 2 || iface == 5 ? -1.0 : 1.0;
     for ( auto j=0; j<igbdy.size(); j++ ) { // all points on face
       for ( auto i=0; i<normals.size(); i++ ) {
         normals[i][j] = 0.0; 
       }
       normals[ip][j] = xm;
       idepComp[j]    = ip; // dependent component
     }
   }
   else if ( this->gtype_ == GE_DEFORMED ) {
     // Bdy normal is dvec{X} / dxi_xi X dvec{X} / dxi_eta
     for ( auto j=0; j<igbdy.size(); j++ ) { // all points on iedge
       ib = igbdy[j];
       for ( auto i=0; i<dXdXi.size(2); i++ ) { // over _X_
         p1[i] = dXdXi(ixi[iface][0],i)[ib]; // d_X_/dxi
         p2[i] = dXdXi(ixi[iface][1],i)[ib]; // d_X_/deta
       }
       p1.cross(p2, xp);   // xp = p1 X p2
       xp.unit(); 
       for ( ic=0; ic<xp.dim(); ic++ ) if ( fabs(xp[ic]) > tiny ) break;
       assert(ic >= GDIM); // no normal components > 0
       for ( auto i=0; i<normals.size(); i++ ) normals[i][j] = xp[i];
       idepComp[j] = ic;  // dependent component
     }
   }

} // end, method do_bdy_normals3d


