//==================================================================================
// Module       : ggrid_box
// Date         : 11/11/18 (DLR)
// Description  : Object defining a (global) 2d or 3d box grid. Builds
//                elements that base class then uses to build global
//                computational data structures.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include "gcomm.hpp"
#include "geoflow.hpp"
#include "ggrid_box.hpp"
#include "tbox/mpixx.hpp"
#include "tbox/global_manager.hpp"


using namespace std;
using namespace geoflow;
using namespace toolbox;

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

  GString gname   = ptree.getValue<GString>("grid_type");
  assert(gname == "grid_box");
  geoflow::tbox::PropertyTree gridptree = ptree.getPropertyTree(gname);


  gbasis_.resize(GDIM);
  gbasis_ = b;
  Lbox_.resize(GDIM);
  ne_.resize(GDIM);

  GTPoint<GFTYPE> gp(ndim_); 
  std::vector<GFTYPE> spt(3); // tmp array
  std::vector  <int> sne   ; // tmp array
  spt = gridptree.getArray<GFTYPE>("xyz0");
  P0_ = spt;
  spt = gridptree.getArray<GFTYPE>("delxyz");
  sne = gridptree.getArray<int>("num_elems");

  eps_ = 100*std::numeric_limits<GFTYPE>::epsilon();
  // compute global bdy range, and global vertices:
  dP_  = spt;
  P1_ = P0_ + dP_;
  gverts_.resize(ndim_);
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

    gp = P0_; gp.x3 += dP_.x3; gverts_[0] = gp; 
    gp.x1 += dP_.x1; gverts_[1] = gp; 
    gverts_[2] = P1_; 
    gp = P0_; gp.x3 += dP.x3; gp.x2 += dP_.x2; gverts_[3] = gp;
  }

  ne_.resize(b.size());
  for ( GSIZET j=0; j<b.size(); j++ ) {
    Lbox_[j] = fabs(dP_[j]);
    ne_  [j] = sne[j];
  }

  lshapefcn_ = new GShapeFcn_linear();
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
void GGridBox::set_partitioner(GDD_base *gdd)
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
  GSIZET icurr = 0; // current global index
  GSIZET fcurr = 0; // current global face index
  GSIZET bcurr = 0; // current global bdy index

  for ( GSIZET i=0; i<qmesh_.size(); i++ ) { // for each quad in irank's mesh
    pelem = new GElem_base(GE_REGULAR, gbasis_);
    xNodes  = &pelem->xNodes();  // node spatial data
    xiNodes = &pelem->xiNodes(); // node ref interval data
    Ni.resize(pelem->nnodes()); // tensor product shape function
    bdy_ind = &pelem->bdy_indices(); // get bdy indices data member
    bdy_typ = &pelem->bdy_types  (); // get bdy types data member
    bdy_ind->clear(); bdy_typ->clear();
    for ( GSIZET l=0; l<ndim_; l++ ) { // loop over element Cart coords
      (*xNodes)[l] = 0.0;
      for ( GSIZET m=0; m<pow(2,ndim_); m++ ) { // loop over verts given in Cart coords
        I[0] = m;
        lshapefcn_->Ni(I, *xiNodes, Ni);
        (*xNodes)[l] += Ni * ( (*(qmesh_[i].v[m]))[l] * 0.25 );
      }
    }

    pelem->init(*xNodes);

    // With face/edge centroids computed, compute 
    // global boundary nodes:
    for ( GSIZET j=0; j<2*ndim_; j++ ) { // cycle over all edges
      cent = pelem->edgeCentroid(j);
      face_ind = NULLPTR;
      if ( FUZZYEQ(P0_.x2,cent.x2,eps_) ) face_ind = &pelem->edge_indices(0);
      if ( FUZZYEQ(P1_.x1,cent.x1,eps_) ) face_ind = &pelem->edge_indices(1);
      if ( FUZZYEQ(P1_.x2,cent.x2,eps_) ) face_ind = &pelem->edge_indices(2);
      if ( FUZZYEQ(P0_.x1,cent.x1,eps_) ) face_ind = &pelem->edge_indices(3);
      // For now, we don't allow corner nodes to be repeated, 
      // and we'll have to choose the best way to define the 
      // normal vectors at the 'corner' nodes:
      for ( GSIZET k=0; face_ind != NULLPTR && k<face_ind->size(); k++ ) {
        if ( !bdy_ind->contains((*face_ind)[k]) ) {
          bdy_ind->push_back((*face_ind)[k]); 
          bdy_typ->push_back(GBDY_NONE); 
        }
      }
    }

    gelems_.push_back(pelem);

    // Set global bdy types at each bdy_ind (this is a coarse 
    // application; finer control may be exercised in callback):
    init_periodic_bdytypes_2d(*pelem);

    // Find global global interior and bdy start & stop indices represented 
    // locally within element:
    nvnodes = gelems_[i]->nnodes();
    nfnodes = gelems_[i]->nfnodes();
    nbnodes = gelems_[i]->bdy_indices().size();
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

#if 0
  // Can set individual nodes and internal bdy conditions
  // with callback here: NOTE: Must re-globalize bdy_indices!!!
  if ( bdycallback_ != NULLPTR ) {
    assert(FALSE && "Re-globalization of bdy data not done");
    (*bdycallback_)(gelems_);
  }
#endif

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

  GSIZET i;
  GSIZET nvnodes;   // no. vol indices
  GSIZET nfnodes;   // no. face indices
  GSIZET nbnodes;   // no. bdy indices
  GSIZET icurr = 0; // current global index
  GSIZET fcurr = 0; // current global face index
  GSIZET bcurr = 0; // current global bdy index
  for ( GSIZET i=0; i<hmesh_.size(); i++ ) { // for each hex in irank's mesh

    pelem = new GElem_base(GE_REGULAR, gbasis_);
    xNodes  = &pelem->xNodes();  // node spatial data
    xiNodes = &pelem->xiNodes(); // node ref interval data
    Ni.resize(pelem->nnodes()); // tensor product shape function
    bdy_ind = &pelem->bdy_indices(); // get bdy indices data member
    bdy_typ = &pelem->bdy_types  (); // get bdy types data member
    bdy_ind->clear(); bdy_typ->clear();
    pelem->igbeg() = icurr;      // beginning global index
    pelem->igend() = icurr + pelem->nnodes()-1; // end global index
    for ( GSIZET l=0; l<ndim_; l++ ) { // loop over element Cart coords
      (*xNodes)[l] = 0.0;
      for ( GSIZET m=0; m<pow(2,ndim_); m++ ) { // loop over verts given in Cart coords
        I[0] = m;
        lshapefcn_->Ni(I, *xiNodes, Ni);
        (*xNodes)[l] += Ni * ( (*(hmesh_[i].v[m]))[l] * 0.125 );
      }
    }

    pelem->init(*xNodes);

    // With edge/face centroids set, compute bdy_nodes:
    for ( GSIZET j=0; j<2*ndim_; j++ ) { // cycle over faces
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
      for ( GSIZET k=0; face_ind != NULLPTR && k<face_ind->size(); k++ ) {
        if ( !bdy_ind->contains((*face_ind)[k]) ) {
          bdy_ind->push_back((*face_ind)[k]); 
          bdy_typ->push_back(GBDY_NONE); // default always to GBDY_NONE 
        }
      }
    }

    gelems_.push_back(pelem);

    // Set global bdy types at each bdy_ind (this is a coarse 
    // application; finer control may be exercised in callback):
    init_periodic_bdytypes_3d(*pelem);

    nvnodes = gelems_[i]->nnodes();
    nfnodes = gelems_[i]->nfnodes();
    nbnodes = gelems_[i]->bdy_indices().size();
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

#if 0
  // Can set individual nodes and internal bdy conditions
  // with callback here:
  if ( bdycallback_ != NULLPTR ) {
    (*bdycallback_)(gelems_);
  }
#endif

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
  for ( GSIZET j=0; j<ppool.size(); j++ ) ppool[j] = gbasis_[j]->getOrder();

  GSIZET iwhere, n;
  GSIZET nvnodes;   // no. vol nodes
  GSIZET nfnodes;   // no. face nodes
  GSIZET nbnodes;   // no. bdy nodes
  GSIZET icurr = 0; // current global index
  GSIZET fcurr = 0; // current global face index
  GSIZET bcurr = 0; // current global bdy index
  for ( GSIZET i=0; i<p.size(1); i++ ) { // for each element
    nvnodes = 1;
    for ( GSIZET j=0; j<GDIM; j++ ) { // set basis from pool
      assert(ppool.contains(p(i,j),iwhere) && "Expansion order not found");
      gb[j] = gbasis_[iwhere];
      nvnodes *= (p(i,j) + 1);
    }
    pelem = new GElem_base(GE_REGULAR, gb);
    xNodes  = &pelem->xNodes();  // node spatial data
    xiNodes = &pelem->xiNodes(); // node ref interval data
    bdy_ind = &pelem->bdy_indices(); // get bdy indices data member
    bdy_typ = &pelem->bdy_types  (); // get bdy types data member
    bdy_ind->clear(); bdy_typ->clear();

    // Set internal node positions from input data.
    // Note that gxnodes are 'global' and xNodes is
    // element-local:
    for ( GSIZET j=0; j<GDIM; j++ ) {
       gxnodes[j].range(icurr, icurr+nvnodes-1);
      (*xNodes)[j] = gxnodes[j];
    }
    for ( GSIZET j=0; j<GDIM; j++ ) gxnodes[j].range_reset();

    pelem->init(*xNodes);


    // With face/edge centroids computed, compute 
    // global boundary nodes:
    for ( GSIZET j=0; j<2*ndim_; j++ ) { // cycle over all edges
      cent = pelem->edgeCentroid(j);
      face_ind = NULLPTR;
      if ( FUZZYEQ(P0_.x2,cent.x2,eps_) ) face_ind = &pelem->edge_indices(0);
      if ( FUZZYEQ(P1_.x1,cent.x1,eps_) ) face_ind = &pelem->edge_indices(1);
      if ( FUZZYEQ(P1_.x2,cent.x2,eps_) ) face_ind = &pelem->edge_indices(2);
      if ( FUZZYEQ(P0_.x1,cent.x1,eps_) ) face_ind = &pelem->edge_indices(3);
      // For now, we don't allow corner nodes to be repeated, 
      // and we'll have to choose the best way to define the 
      // normal vectors at the 'corner' nodes:
      for ( GSIZET k=0; face_ind != NULLPTR && k<face_ind->size(); k++ ) {
        if ( !bdy_ind->contains((*face_ind)[k]) ) {
          bdy_ind->push_back((*face_ind)[k]); 
          bdy_typ->push_back(GBDY_NONE); 
        }
      }
    }


    gelems_.push_back(pelem);

    // Set global bdy types at each bdy_ind (this is a coarse 
    // application; finer control may be exercised in callback):
    init_periodic_bdytypes_2d(*pelem);

    // Find global global interior and bdy start & stop indices represented 
    // locally within element:
    assert(nvnodes == gelems_[n]->nnodes() && "Incompatible node count");
    nfnodes = gelems_[i]->nfnodes();
    nbnodes = gelems_[i]->bdy_indices().size();
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


#if 0
  // Can set individual nodes and internal bdy conditions
  // with callback here: NOTE: Must re-globalize bdy_indices!!!
  if ( bdycallback_ != NULLPTR ) {
    assert(FALSE && "Re-globalization of bdy data not done");
    (*bdycallback_)(gelems_);
  }
#endif

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
  for ( GSIZET j=0; j<ppool.size(); j++ ) ppool[j] = gbasis_[j]->getOrder();

  GSIZET i, iwhere, n;
  GSIZET nvnodes;   // no. vol indices
  GSIZET nfnodes;   // no. face indices
  GSIZET nbnodes;   // no. bdy indices
  GSIZET icurr = 0; // current global index
  GSIZET fcurr = 0; // current global face index
  GSIZET bcurr = 0; // current global bdy index
  for ( GSIZET i=0; i<p.size(1); i++ ) { // for each hex in irank's mesh
    nvnodes = 1; 
    for ( GSIZET j=0; j<GDIM; j++ ) { // set basis from pool
      assert(ppool.contains(p(i,j),iwhere) && "Expansion order not found");
      gb[j] = gbasis_[iwhere];
      nvnodes *= (p(i,j) + 1);
    }    

    pelem = new GElem_base(GE_REGULAR, gbasis_);
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
    for ( GSIZET j=0; j<GDIM; j++ ) {
       gxnodes[j].range(icurr, icurr+nvnodes-1);
      (*xNodes)[j] = gxnodes[j];
    }
    for ( GSIZET j=0; j<GDIM; j++ ) gxnodes[j].range_reset();

    pelem->init(*xNodes);

    // With edge/face centroids set, compute bdy_nodes:
    for ( GSIZET j=0; j<2*ndim_; j++ ) { // cycle over faces
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
      for ( GSIZET k=0; face_ind != NULLPTR && k<face_ind->size(); k++ ) {
        if ( !bdy_ind->contains((*face_ind)[k]) ) {
          bdy_ind->push_back((*face_ind)[k]); 
          bdy_typ->push_back(GBDY_NONE); // default always to GBDY_NONE 
        }
      }
    }

    gelems_.push_back(pelem);

    // Set global bdy types at each bdy_ind (this is a coarse 
    // application; finer control may be exercised in callback):
    init_periodic_bdytypes_3d(*pelem);

    assert(nvnodes == gelems_[i]->nnodes() && "Incompatible node count");
    nfnodes = gelems_[i]->nfnodes();
    nbnodes = gelems_[i]->bdy_indices().size();
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

#if 0
  // Can set individual nodes and internal bdy conditions
  // with callback here:
  if ( bdycallback_ != NULLPTR ) {
    (*bdycallback_)(gelems_);
  }
#endif

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
    for ( GSIZET i=0; i<qmesh_.size(); i++ ) { // for each quad/hex
      for ( GSIZET j=0; j<pow(2,ndim_); j++ ) { // for each vertex of reg polygon
          pt = *qmesh_[i].v[j];
          ios << pt.x1 << " " <<  pt.x2 << std::endl;
          ios << pt.x1 << " " << pt.x2 << " " <<  pt.x3 << std::endl;
      }
    }
  } 

  if ( ndim_ == 3 ) {
    for ( GSIZET i=0; i<hmesh_.size(); i++ ) { // for each quad/hex
      for ( GSIZET j=0; j<pow(2,ndim_); j++ ) { // for each vertex of reg polygon
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
  for ( GSIZET k=0; k<x.size(); k++ ) x[k] = P0_[k];

  GUINT  bit;
  GSIZET id, n=0;
  for ( GSIZET k=0; k<igbdy_[GBDY_PERIODIC].size(); k++ ) { // for each blobal bdy node
    id = igbdy_[GBDY_PERIODIC][k];
    periodicids_.push_back(id);       
    periodicdirs_.push_back(0);
    for ( GSIZET i=0; i<xNodes_.size(); i++ ) { // for x, y, z dirs
      if ( FUZZYEQ(P1_[i],xNodes_[i][id],eps_) ) { // right/top-most coord will change
        periodicdirs_[n] |= 1U << i;  // position right-most direction bit  
      }
    }
    n++;
  }

  // Now, cycle through periodic nodes and periodize coordinates:
  for ( GSIZET k=0; k<periodicids_.size(); k++ ) { // for each periodic node
    id = periodicids_[k];
    for ( GSIZET i= 0; i<xNodes_.size(); i++ ) { // coord direction
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
  for ( GSIZET k=0; k<x.size(); k++ ) x[k] = P1_[k];

  // Cycle through periodic nodes and un-periodize coordinates:
  GUINT  bit;
  GSIZET id;
  for ( GSIZET k=0; k<periodicids_.size(); k++ ) { // for each periodic node
    id = periodicids_[k];
    for ( GSIZET i= 0; i<xNodes_.size(); i++ ) { // coord direction
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
  for ( GSIZET k=0; k<ne_.size(); k++ ) nglobal *= ne_[k];
 
  nperrank = nglobal / nprocs_; // #elems per rank
  nthisrank = irank_ != nprocs_-1 ? nperrank : nglobal - (nprocs_-1)*nperrank;


 // Get uniform element sizes:
  for ( GSIZET k=0; k<ndim_; k++ ) {
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
//          igbdy : which nodes indices represent global boundaries. This contains
//                  indices for all bdys this task owns
//          igbdyt: bdy type ids for each index in igbdy
// RETURNS: none.
//**********************************************************************************
void GGridBox::config_bdy(const PropertyTree &ptree, 
                          GTVector<GSIZET> &igbdy, 
                          GTVector<GSIZET> &igbdyt)
{
  // Cycle over all geometric boundaries, and configure:

  GBOOL              bperiodic=FALSE, buniform=FALSE, buseinit;
  GSIZET             iwhere;
  GTVector<GBOOL>    uniform(2*GDIM);
  GTVector<GBdyType> bdytype(2*GDIM);
  GTVector<GBdyType> btmp;
  GTVector<GSIZET>   itmp;
  GTVector<GString>  bdynames (2*GDIM);
  GTVector<GString>  confmthd (2*GDIM);
  GTVector<GString>  bdyupdate(2*GDIM);
  GString            gname, sbdy, bdyclass, bdyinit;
  PropertyTree       bdytree, gridptree, spectree;

  assert(gname == "grid_box");

  bdynames[0] = "bdy_x_0";
  bdynames[1] = "bdy_x_1";
  bdynames[2] = "bdy_y_0";
  bdynames[3] = "bdy_y_1";
  bdynames[4] = "bdy_z_0";
  bdynames[5] = "bdy_z_1";

  gname     = ptree.getValue<GString>("grid_type");
  gridptree = ptree.getPropertyTree(gname);


  // Clear input arrays:
  igbdy .clear();
  igbdyt.clear();

  bdyupdate = gridtree.getValue<GString>("update_method","");
  bdyinit   = gridtree.getValue<GString>("bdy_init_method","");
  buseinit  = gridtree.getValue<GBOOL>  ("use_state_init_method",FALSE);


  // Get properties from the main prop tree. 
  // Note: bdys are configured by way of geometry's
  //       natural decomposition: here, by face (3d) or
  //       edge (2d). But the bdy indices and types
  //       returned on exist contain info for all bdys:
  for ( auto j=0; j<2*GDIM; j++ ) { // cycle over faces
    sbdy         = gridptree.getValue<GString>(bdynames[j]);
    bdytree      = gridtree.getPropertyTree(sbdy);
    bdyclass     = bdytree.getValue<GString>("bdy_class", "uniform");
    bdytype  [j] = geoflow::str2bdytype(bdytree.getValue<GString>("base_type", GBDY_NONE));
    uniform  [j] = bdyclass == "uniform" ? TRUE : FALSE;
    confmthd  [j] = bdytree.getValue<GString>("bdy_config_method","");
    bperiodic    = bperiodic || bdytype[j] == GBDY_PERIODIC;
    buniform     = buniform || uniform[j];
    assert(bperiodic && !buniform[j] && "GBDY_PERIODIC boundary must have bdy_class = uniform");
  }

  assert( (  (bdytype[0] == GBDY_PERIODIC && bdytype[2] != GBDY_PERIODIC)
         ||  (bdytype[3] == GBDY_PERIODIC && bdytype[1] != GBDY_PERIODIC) )
         &&  "Incompatible GBDY_PERIODIC boundary specification");
       
  // Handle non-uniform (user-configured) bdy types first;
  // Note: If "uniform" not specified for a boundary, then
  //       user MUST supply a method to configure it.
  //       Also, each natural face may be configured independently,
  //       but the bdy indices & corresp. types are concatenated into 
  //       single arrays:
  for ( auto j=0; j<2*GDIM; j++ ) { 
    // First, find global bdy indices:
    if ( buniform[j] ) continue;
    if ( ndim_ == 2 ) {
      find_bdy_indices2d(j, TRUE, itmp); // include vertices
    }
    else {
      find_bdy_indices3d(j, TRUE, itmp); // include edges
    }
    spectree  = ptree->getPropertyTree(confmthd[j]);
    bret = GSpecBFactory::dospec(spectree, *this, j, itmp, btmp); // get user-defined bdy spec
    assert(bret && "Boundary specification failed");
    igbdy .concat(itmp.data(), itmp.size());
    igbdyt.concat(btmp.data(), btmp.size());
  }
  
  // Fill in uniform bdy types:
  for ( auto j=0; j<2*GDIM; j++ ) { 
    if ( !buniform[j] ) continue;
    // First, find global bdy indices:
    if ( bperiodic && bdytype[j] != GBDY_PERIODIC  ) {
      if ( ndim_ == 2 ) {
        find_bdy_indices2d(j, FALSE, itmp); // doesn't include vertices
      }
      else {
        find_bdy_indices3d(j, FALSE, itmp); // doesn't include edges
      }

    }
    else {
      if ( ndim_ == 2 ) {
        find_bdy_indices2d(j, TRUE, itmp); // include vertices
      }
      else {
        find_bdy_indices3d(j, TRUE, itmp); // include edges
      }
    }
    // Set type for each bdy index:
    for ( auto i=0; i<itmp.size(); i++ ) {
      btmp[i] = bdytype[j]; 
    }
    igbdy .concat(itmp.data(), itmp.size());
    igbdyt.concat(btmp.data(), btmp.size());
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
      for ( GSIZET i=0; i<xNodes[0].size(); i++ ) { 
        if ( FUZZYEQ(P0_.x2,xNodes[1][i],eps_) ) {
          pt.assign(xNodes_, i);
          if ( incl_vert  || !is_global_vertex(pt) ) ibdy.push_back(i);
        }
      }
      break;

    case 1: // right vert bdy:
      for ( GSIZET i=0; i<xNodes[0].size(); i++ ) { // bdy 1
        if ( FUZZYEQ(P1_.x1,xNodes[0][i],eps_) ) {
          pt.assign(xNodes_, i);
          if ( incl_vert  || !is_global_vertex(pt) ) ibdy.push_back(i);
        }
      }
      break;

    case 2: // top horiz bdy:
      for ( GSIZET i=0; i<xNodes[0].size(); i++ ) { // bdy 2
        if ( FUZZYEQ(P1_.x2,xNodes[1][i],eps_) ) {
          pt.assign(xNodes_, i);
          if ( incl_vert  || !!is_global_vertex(pt) ) ibdy.push_back(i);
        }
      }
      break;

    case 3: // left vert bdy:
      for ( GSIZET i=0; i<xNodes[0].size(); i++ ) { // bdy 3
        if ( FUZZYEQ(P0_.x1,xNodes[0][i],eps_) ) {
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

  ibdy.clear();

  GTPoint<GFTYPE> pt(ndim_);

  pt.setBracket(eps_);

  switch ( bdyid ) {
    case 0: // southern vert face:
      for ( GSIZET i=0; i<xNodes[0].size(); i++ ) { // face 0
        if ( FUZZYEQ(P0_.x1,xNodes[0][i],eps_) ) {
          pt.assign(xNodes_, i);
          if ( incl_edge || !on_global_edge(0,pt) ) ibdy.push_back(i);
        }
      }
      break;

    case 1: // eastern vert. face:
      for ( GSIZET i=0; i<xNodes[0].size(); i++ ) { // face 1
        if ( FUZZYEQ(P1_.x1,xNodes[1][i],eps_) ) {
          pt.assign(xNodes_, i);
          if ( incl_edge || !on_global_edge(1,pt) ) ibdy.push_back(i);
        }
      }
      break;

    case 2: // northern vert. face:
      for ( GSIZET i=0; i<xNodes[0].size(); i++ ) { // face 2
        if ( FUZZYEQ(P1_.x1,xNodes[0][i],eps_) ) {
          pt.assign(xNodes_, i);
          if ( incl_edge || !on_global_edge(2,pt) ) ibdy.push_back(i);
        }
      }
      break;

    case 3: // western vertical face:
      for ( GSIZET i=0; i<xNodes[0].size(); i++ ) { // face 3
        if ( FUZZYEQ(P0_.x1,xNodes[1][i],eps_) ) {
          pt.assign(xNodes_, i);
          if ( incl_edge || !on_global_edge(3,pt) ) ibdy.push_back(i);
        }
      }
      break;

    case 4: // bottom horiz face:
      for ( GSIZET i=0; i<xNodes[0].size(); i++ ) { // face 4
        if ( FUZZYEQ(P0_.x1,xNodes[2][i],eps_) ) {
          pt.assign(xNodes_, i);
          if ( incl_edge || !on_global_edge(4,pt) ) ibdy.push_back(i);
        }
      }
      break;

    case 5: // top horiz face:
      for ( GSIZET i=0; i<xNodes[0].size(); i++ ) { // face 5
        if ( FUZZYEQ(P1_.x1,xNodes[2][i],eps_) ) {
          pt.assign(xNodes_, i);
          if ( incl_edge || !on_global_edge(5,pt) ) ibdy.push_back(i);
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

  for ( GSIZET j=0; j<pow(2,ndim_) && !bret; j++ ) {
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
  GTPoint<GFTYPE> pt(ndim_);

  // Find faces point belongs to:
  for ( GINT j=0; j<ndim_; j++ ) {
    if     ( FUZZYEQ(P0_.x1,eps_) && !face.containsn(0,nface) ) 
      { face[nface] = 0; nface++; }
    else if( FUZZYEQ(P1_.x2,eps_) && !face.containsn(1,nface) ) 
      { face[nface] = 1; nface++; } 
    else if( FUZZYEQ(P1_.x1,eps_) && !face.containsn(2,nface) ) 
      { face[nface] = 2; nface++; }
    else if( FUZZYEQ(P0_.x2,eps_) && !face.containsn(3,nface) ) 
      { face[nface] = 3; nface++; }
    else if( FUZZYEQ(P0_.x3,eps_) && !face.containsn(4,nface) ) 
      { face[nface] = 4; nface++; }
    else if( FUZZYEQ(P1_.x3,eps_) && !face.containsn(5,nface) ) 
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
