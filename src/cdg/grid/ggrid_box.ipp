//==================================================================================
// Module       : ggrid_box
// Date         : 11/11/18 (DLR)
// Description  : Object defining a (global) 2d or 3d box grid. Builds
//                elements that base class then uses to build global
//                computational data structures.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================


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
template<typename Types>
GGridBox<Types>::GGridBox(const geoflow::tbox::PropertyTree &ptree, GTVector<GNBasis<GCTYPE,Ftype>*> &b, GC_COMM &comm)
:   GGrid<Types>(ptree, b, comm),
ndim_                     (GDIM),
gdd_                   (NULLPTR),
lshapefcn_             (NULLPTR)
{
  GEOFLOW_TRACE();
  assert((b.size() == GDIM ) 
        && "Basis has incorrect dimensionalilty");

  this->irank_  = GComm::WorldRank(this->comm_);
  this->nprocs_ = GComm::WorldSize(this->comm_);

  GString gname   = ptree.getValue<GString>("grid_type");
  GString tname   = ptree.getValue<GString>("terrain_type","");
  assert(gname == "grid_box");
  geoflow::tbox::PropertyTree gridptree = ptree.getPropertyTree(gname);

  // A hack: determine if ptree wants us to compute
  // bdy test data:
  this->do_gbdy_test_ = ptree.getValue<GBOOL>("gbdy_test_data",FALSE);

  // If terrain is being used, elements may not be
  // GE_REGULAR:
  this->gtype_ = GE_REGULAR;
  if ( "none" != tname 
   &&  ""     != tname ) this->gtype_ = GE_DEFORMED;

  gbasis_.resize(GDIM);
  gbasis_ = b;
  Lbox_.resize(GDIM);
  ne_.resize(GDIM);

  GTPoint<Ftype> gp(ndim_); 
  std::vector<Ftype> spt(3); // tmp array
  std::vector  <int> sne   ; // tmp array
  spt = gridptree.getArray<Ftype>("xyz0");
  P0_.resize(GDIM);
  P1_.resize(GDIM);
  dP_.resize(GDIM);
  for ( auto j=0; j<GDIM; j++ ) P0_[j] = spt[j];
  spt = gridptree.getArray<Ftype>("delxyz");
  for ( auto j=0; j<GDIM; j++ ) dP_[j] = spt[j];
  sne = gridptree.getArray<int>("num_elems");

  // Compute global bdy range, and global vertices:
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

  lshapefcn_ = new GShapeFcn_linear<Ftype>(GDIM);
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
template<typename Types>
GGridBox<Types>::~GGridBox()
{
  GEOFLOW_TRACE();
  if ( lshapefcn_ != NULLPTR ) delete lshapefcn_;
  if ( gdd_       != NULLPTR ) delete gdd_;
} // end, destructor


#if 0
//**********************************************************************************
//**********************************************************************************
// METHOD :  << operator method (1)
// DESC   : output stream operator
// ARGS   :
// RETURNS: ostream &
//**********************************************************************************
template<typename Types>
std::ostream &operator<<(std::ostream &str, GGridBox<Types> &e)
{
  GEOFLOW_TRACE();
  
  str << "    Lbox: " << e.Lbox_;
  str << std::endl << " Centroids: " ;
  for( auto i=0; i<e.ftcentroids_.size(); i++ )
     str << (e.ftcentroids_[i]) << " " ;
  str << std::endl;

  return str;
} // end of operator <<
#endif


//**********************************************************************************
//**********************************************************************************
// METHOD : set_partitioner
// DESC   : Set domain decomposition object
// ARGS   : GDD_base pointer
// RETURNS: none
//**********************************************************************************
template<typename Types>
void GGridBox<Types>::set_partitioner(GDD_base<Ftype> *gdd)
{
  GEOFLOW_TRACE();

  gdd_ = gdd;

} // end of method set_partitioner


//**********************************************************************************
//**********************************************************************************
// METHOD : init2d
// DESC   : Initialize base state/base icosahedron
// ARGS   : none.
// RETURNS: none.
//**********************************************************************************
template<typename Types>
void GGridBox<Types>::init2d()
{
  GEOFLOW_TRACE();

  find_rank_subdomain();

} // end of method init2d


//**********************************************************************************
//**********************************************************************************
// METHOD : init3d
// DESC   : Initialize for 3d elements
// ARGS   : none.
// RETURNS: none.
//**********************************************************************************
template<typename Types>
void GGridBox<Types>::init3d()
{
  GEOFLOW_TRACE();

  find_rank_subdomain();

} // end, method init3d



//**********************************************************************************
//**********************************************************************************
// METHOD : do_elems (1)
// DESC   : Public entry point for grid element computation
// ARGS   : none.
// RETURNS: none.
//**********************************************************************************
template<typename Types>
void GGridBox<Types>::do_elems()
{
  GEOFLOW_TRACE();
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
template<typename Types>
void GGridBox<Types>::do_elems(GTMatrix<GINT> &p,
                        GTVector<GTVector<Ftype>> &xnodes)
{
  GEOFLOW_TRACE();
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
template<typename Types>
void GGridBox<Types>::do_elems2d()
{
  GEOFLOW_TRACE();
  assert(gbasis_.size()>0 && "Must set basis first");
  assert(ndim_ == 2 && "Dimension must be 2");

  GString                      serr = "GGridBox<Types>::do_elems2d (1): ";
  GTPoint<Ftype>              cent;
  GTVector<GINT>               iind;
  GTVector<GINT>               I(3);
  GTVector<GINT>              *bdy_ind;
  GTVector<GBdyType>          *bdy_typ;
  GTVector<GINT>              *face_ind;
  GTVector<Ftype>             Ni;
  GElem_base                  *pelem;
  GTVector<GTVector<Ftype>>  *xNodes;
  GTVector<GTVector<Ftype>*> *xiNodes;
  GTVector<GTVector<Ftype>>   xgtmp(3);


#if 0
  if ( gdd_       == NULLPTR ) gdd_ = new GDD_base(this->nprocs_);

  // Get indirection indices to create elements
  // only for this task:
  gdd_->doDD(ftcentroids_, this->irank_, iind);
#endif

  GSIZET nvnodes;   // no. vol nodes
  GSIZET nfnodes;   // no. face nodes
  GSIZET nbnodes;   // no. bdy nodes
  GLONG  icurr = 0; // current global index
  GLONG  fcurr = 0; // current global face index
  GLONG  bcurr = 0; // current global bdy index

  for ( auto i=0; i<qmesh_.size(); i++ ) { // for each quad in irank's mesh
    pelem = new GElem_base(GDIM, this->gtype_, gbasis_);
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

#if 0
    // With face/edge centroids computed, compute 
    // global boundary nodes:
    for ( auto j=0; j<2*ndim_; j++ ) { // cycle over all edges
      cent = pelem->edgeCentroid(j);
      face_ind = NULLPTR;
      if ( FUZZYEQ(P0_.x2,cent.x2,this->eps_) ) face_ind = &pelem->edge_indices(0);
      if ( FUZZYEQ(P1_.x1,cent.x1,this->eps_) ) face_ind = &pelem->edge_indices(1);
      if ( FUZZYEQ(P1_.x2,cent.x2,this->eps_) ) face_ind = &pelem->edge_indices(2);
      if ( FUZZYEQ(P0_.x1,cent.x1,this->eps_) ) face_ind = &pelem->edge_indices(3);
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
#endif


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

    this->gelems_.push_back(pelem);
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
template<typename Types>
void GGridBox<Types>::do_elems3d()
{
  GEOFLOW_TRACE();
  assert(gbasis_.size()>0 && "Must set basis first");
  assert(ndim_ == 3 && "Dimension must be 3");

  GString                      serr = "GGridBox<Types>::do_elems3d (1): ";
  GTPoint<Ftype>              cent;
  GTVector<GINT>               iind;
  GTVector<GINT>               I(3);
  GTVector<GINT>              *bdy_ind;
  GTVector<GBdyType>          *bdy_typ;
  GTVector<GINT>              *face_ind;
  GTVector<Ftype>             Ni;
  GElem_base                  *pelem;
  GTVector<GTVector<Ftype>>  *xNodes;
  GTVector<GTVector<Ftype>*> *xiNodes;
  GTVector<GTVector<Ftype>>   xgtmp(3);


#if 0
  if ( gdd_       == NULLPTR ) gdd_ = new GDD_base(this->nprocs_);

  // Get indirection indices to create elements
  // only for this task:
  gdd_->doDD(ftcentroids_, this->irank_, iind);
#endif

  GSIZET nvnodes;   // no. vol indices
  GSIZET nfnodes;   // no. face indices
  GSIZET nbnodes;   // no. bdy indices
  GLONG  icurr = 0; // current global index
  GLONG  fcurr = 0; // current global face index
  GLONG  bcurr = 0; // current global bdy index
  for ( auto i=0; i<hmesh_.size(); i++ ) { // for each hex in irank's mesh

    pelem = new GElem_base(GDIM, this->gtype_, gbasis_);
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
      if ( FUZZYEQ(P0_.x2,cent.x2,this->eps_) ) face_ind = &pelem->edge_indices(0);
      if ( FUZZYEQ(P1_.x1,cent.x1,this->eps_) ) face_ind = &pelem->edge_indices(1);
      if ( FUZZYEQ(P1_.x2,cent.x2,this->eps_) ) face_ind = &pelem->edge_indices(2);
      if ( FUZZYEQ(P0_.x1,cent.x1,this->eps_) ) face_ind = &pelem->edge_indices(3);
      if ( FUZZYEQ(P0_.x3,cent.x3,this->eps_) ) face_ind = &pelem->edge_indices(4);
      if ( FUZZYEQ(P1_.x3,cent.x3,this->eps_) ) face_ind = &pelem->edge_indices(5);
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

    this->gelems_.push_back(pelem);

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
template<typename Types>
void GGridBox<Types>::do_elems2d(GTMatrix<GINT> &p,
                          GTVector<GTVector<Ftype>> &gxnodes)
{
  GEOFLOW_TRACE();
  assert(gbasis_.size()>0 && "Must set basis pool first");
  assert(ndim_ == 2 && "Dimension must be 2");

  GString                      serr = "GGridBox<Types>::do_elems2d (2): ";
  GTPoint<Ftype>              cent;
  GTVector<GINT>               I(3);
  GTVector<GINT>              *bdy_ind;
  GTVector<GBdyType>          *bdy_typ;
  GTVector<GINT>              *face_ind;
  GElem_base                  *pelem;
  GTVector<GTVector<Ftype>>  *xNodes;
  GTVector<GTVector<Ftype>*> *xiNodes;
  GTVector<GTVector<Ftype>>   xgtmp(3);
  GTVector<GNBasis<GCTYPE,Ftype>*>
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
    pelem = new GElem_base(GDIM, this->gtype_, gb);
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

#if 0
    // With face/edge centroids computed, compute 
    // global boundary nodes:
    for ( auto j=0; j<2*ndim_; j++ ) { // cycle over all edges
      cent = pelem->edgeCentroid(j);
      face_ind = NULLPTR;
      if ( FUZZYEQ(P0_.x2,cent.x2,this->eps_) ) face_ind = &pelem->edge_indices(0);
      if ( FUZZYEQ(P1_.x1,cent.x1,this->eps_) ) face_ind = &pelem->edge_indices(1);
      if ( FUZZYEQ(P1_.x2,cent.x2,this->eps_) ) face_ind = &pelem->edge_indices(2);
      if ( FUZZYEQ(P0_.x1,cent.x1,this->eps_) ) face_ind = &pelem->edge_indices(3);
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
#endif

    this->gelems_.push_back(pelem);

    // Find global global interior and bdy start & stop indices represented 
    // locally within element:
    assert(nvnodes == this->gelems_[n]->nnodes() && "Incompatible node count");
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
template<typename Types>
void GGridBox<Types>::do_elems3d(GTMatrix<GINT> &p,
                          GTVector<GTVector<Ftype>> &gxnodes)
{
  GEOFLOW_TRACE();
  assert(gbasis_.size()>0 && "Must set basis pool first");
  assert(ndim_ == 3 && "Dimension must be 3");

  GString                      serr = "GGridBox<Types>::do_elems3d (2): ";
  GTPoint<Ftype>              cent;
  GTVector<GINT>               iind;
  GTVector<GINT>               I(3);
  GTVector<GINT>              *bdy_ind;
  GTVector<GBdyType>          *bdy_typ;
  GTVector<GINT>              *face_ind;
  GTVector<Ftype>             Ni;
  GElem_base                  *pelem;
  GTVector<GTVector<Ftype>>  *xNodes;
  GTVector<GTVector<Ftype>*> *xiNodes;
  GTVector<GTVector<Ftype>>   xgtmp(3);
  GTVector<GNBasis<GCTYPE,Ftype>*>
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

    pelem = new GElem_base(GDIM, this->gtype_, gbasis_);
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
      if ( FUZZYEQ(P0_.x2,cent.x2,this->eps_) ) face_ind = &pelem->edge_indices(0);
      if ( FUZZYEQ(P1_.x1,cent.x1,this->eps_) ) face_ind = &pelem->edge_indices(1);
      if ( FUZZYEQ(P1_.x2,cent.x2,this->eps_) ) face_ind = &pelem->edge_indices(2);
      if ( FUZZYEQ(P0_.x1,cent.x1,this->eps_) ) face_ind = &pelem->edge_indices(3);
      if ( FUZZYEQ(P0_.x3,cent.x3,this->eps_) ) face_ind = &pelem->edge_indices(4);
      if ( FUZZYEQ(P1_.x3,cent.x3,this->eps_) ) face_ind = &pelem->edge_indices(5);
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

    this->gelems_.push_back(pelem);

    assert(nvnodes == this->gelems_[i]->nnodes() && "Incompatible node count");
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
template<typename Types>
void GGridBox<Types>::print(const GString &filename)
{
  GEOFLOW_TRACE();
  GString serr = "GGridBox<Types>::print: ";
  std::ofstream ios;

  Ftype x, y, z;
  GTPoint<Ftype> pt;

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
template<typename Types>
void GGridBox<Types>::periodize()
{
  GEOFLOW_TRACE();
  assert(this->bInitialized_ && "Object not initialized");

  
  GTVector<Ftype>  x(this->xNodes_.size()); // coord values to set to

  assert(this->minnodedist_ > 0.0);
  this->eps_ = 0.00125*this->minnodedist_;

  periodicids_.clear();
  periodicdirs_.clear();

  // Coords set to correspond to bottom-left-most domain point:
  for( auto k=0; k<x.size(); k++ ) x[k] = P0_[k];

  GUINT  bit;
  GSIZET id, n, num=0;
  for ( auto k=0; k<this->igbdy_binned_.size(); k++ ) {
    num += this->igbdy_binned_[k][GBDY_PERIODIC].size();
  }
  periodicids_ .resize(num);
  periodicdirs_.resize(num);

  n = 0;
  for( auto k=0; k<this->igbdy_binned_.size(); k++ ) { // for each global face
    for( auto j=0; j<this->igbdy_binned_[k][GBDY_PERIODIC].size(); j++, n++ ) { // for each global bdy node
      id = this->igbdy_binned_[k][GBDY_PERIODIC][j];
      periodicids_ [n] = id;       
      periodicdirs_[n] = 0;
      for( auto i=0; i<this->xNodes_.size(); i++ ) { // for x, y, z dirs
        if ( FUZZYEQ(P1_[i],this->xNodes_[i][id],this->eps_) ) { // right/top-most block.tbdy[k];i coord will change
          periodicdirs_[n] |= 1U << i;  // position right-most direction bit  
        }
      }
    }
  }

  // Now, cycle through periodic nodes and periodize coordinates:
  for( auto k=0; k<periodicids_.size(); k++ ) { // for each periodic node
    id = periodicids_[k];
    for( auto i= 0; i<this->xNodes_.size(); i++ ) { // coord direction
      // Set coord in this direction if corresp bit is set:
      bit = (periodicdirs_[k] >> i) & 1; 
      if ( bit ) this->xNodes_[i][id] = x[i];
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
template<typename Types>
void GGridBox<Types>::unperiodize()
{
  GEOFLOW_TRACE();
  // Use data from 'periodize' method to unset change in
  // coordinates:

  GTVector<Ftype>  x(this->xNodes_.size()); // coord values to set to

  // Coords to set to correspond to top-most domain point:
  for( auto k=0; k<x.size(); k++ ) x[k] = P1_[k];

  // Cycle through periodic nodes and un-periodize coordinates:
  GUINT  bit;
  GSIZET id;
  for( auto k=0; k<periodicids_.size(); k++ ) { // for each periodic node
    id = periodicids_[k];
    for( auto i= 0; i<this->xNodes_.size(); i++ ) { // coord direction
      // Set coord in this direction if corresp bit is set:
      bit = (periodicdirs_[k] >> i) & 1; 
      if ( bit ) this->xNodes_[i][id] = x[i];
    }
  }

  periodicids_.clear();
  periodicdirs_.clear();

} // end of method unperiodize

//**********************************************************************************
//**********************************************************************************
// METHOD : max_duplicates
// DESC   : Maximum duplicate points found within the grid
// ARGS   : none.
// RETURNS: Maximum duplicate
//**********************************************************************************
template<typename Types>
std::size_t GGridBox<Types>::max_duplicates() const {
  return std::pow(2,ndim_);
}

//**********************************************************************************
//**********************************************************************************
// METHOD : find_rank_subdomain
// DESC   : Find this rank's default portion of global domain, and
//          store in qmesh member data.
// ARGS   : none.
// RETURNS: none.
//**********************************************************************************
template<typename Types>
void GGridBox<Types>::find_rank_subdomain()
{
  GEOFLOW_TRACE();
  GSIZET          n=0, nglobal, nperrank, nthisrank, ntot, nxy;
  GLONG           beg_lin, end_lin;
  GLONG           ib, ie, jb, je, kb, ke;
  GTPoint<Ftype> v0(ndim_);     // starting sph. point of elem
  GTPoint<Ftype> dv(ndim_);     // delta-point
  GTPoint<Ftype> dx(ndim_);

  nglobal = 1; // total num elements in global grid
  for( auto k=0; k<ne_.size(); k++ ) nglobal *= ne_[k];
 
  nperrank = nglobal / this->nprocs_; // #elems per rank
  nthisrank = this->irank_ != this->nprocs_-1 ? nperrank : nglobal - (this->nprocs_-1)*nperrank;


 // Get uniform element sizes:
  for( auto k=0; k<ndim_; k++ ) {
    dx[k] = Lbox_[k] / static_cast<Ftype>(ne_[k]);
  }

  if ( this->do_gbdy_test_ ) {
    this->testty_.resize(nthisrank); // elem type
    this->testid_.resize(nthisrank);  // id for type
  }

  // Compute vertices of all hexes ('cubes')
  // for this task:
  if ( ndim_ == 2 ) {

    qmesh_.resize(nthisrank);
    beg_lin = nperrank*this->irank_; end_lin = beg_lin + nthisrank - 1;
    jb = beg_lin/ne_[0]; je = end_lin/ne_[0];
    for ( auto j=jb; j<=je; j++ ) {
      ib = MAX(beg_lin-static_cast<GLONG>(j*ne_[0]),0); 
      ie = MIN(end_lin-static_cast<GLONG>(j*ne_[0]),ne_[0]-1); 
      for ( auto i=ib; i<=ie; i++ ) {
        for ( auto l=0; l<4; l++ ) qmesh_[n][l].resize(ndim_);
        v0.x1 = P0_.x1+i*dx.x1; v0.x2 = P0_.x2+j*dx.x2;
                                         qmesh_[n].v1 = v0;
        dv.x1 = dx.x1 ; dv.x2 = 0.0  ;   qmesh_[n].v2 = v0 + dv;
        dv.x1 = dx.x1 ; dv.x2 = dx.x2;   qmesh_[n].v3 = v0 + dv;
        dv.x1 = 0.0   ; dv.x2 = dx.x2;   qmesh_[n].v4 = v0 + dv;

if ( this->do_gbdy_test_ ) {
if      ( i == 0        && j == 0        ) {
  this->testty_[n] = 0; // corner elem
  this->testid_[n] = 0; 
}
else if ( i == ne_[0]-1 && j == 0        ) {
  this->testty_[n] = 0; // corner elem
  this->testid_[n] = 1; 
}
else if ( i == ne_[0]-1 && j == ne_[1]-1 ) {
  this->testty_[n] = 0; // corner elem
  this->testid_[n] = 2; 
}
else if ( i == 0        && j == ne_[1]-1 ) {
  this->testty_[n] = 0; // corner elem
  this->testid_[n] = 3; 
}
else if ( j == 0 ) {
  this->testty_[n] = 1; // non-corner edge
  this->testid_[n] = 0; 
}
else if ( i == ne_[0]-1 ) {
  this->testty_[n] = 1; // non-corner edge elem
  this->testid_[n] = 1; 
}
else if ( j == ne_[1]-1 ) {
  this->testty_[n] = 1; // non-corner edge elem
  this->testid_[n] = 2; 
}
else if ( i == 0 ) {
  this->testty_[n] = 1; // non-corner edge elem
  this->testid_[n] = 3; 
}
else {
  this->testty_[n] = 2; // interior elem
  this->testid_[n] = -1; // interior id
}
} // end, do_gbdy_test_
        n++;
      }
    }
  } // end, ndim==2 test
  else if ( ndim_ == 3 ) {
    hmesh_.resize(nthisrank);
    nxy = ne_[0] * ne_[1];
    beg_lin = nperrank*this->irank_; end_lin = beg_lin + nthisrank - 1;
    kb = beg_lin/nxy; ke = end_lin/nxy;
    for ( auto k=kb; k<=ke; k++ ) { 
      jb = MAX((beg_lin-static_cast<GLONG>(k*nxy))/ne_[0],0); 
      je = MIN((end_lin-static_cast<GLONG>(k*nxy))/ne_[0],ne_[1]-1);
      for ( auto j=jb; j<=je; j++ ) { 
        ib = MAX(beg_lin-static_cast<GLONG>(k*nxy+j*ne_[0]),0); 
        ie = MIN(end_lin-static_cast<GLONG>(k*nxy+j*ne_[0]),ne_[0]-1); 
        for ( auto i=ib; i<=ie; i++ ) { 
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


} // end of method find_rank_subdomain


//**********************************************************************************
//**********************************************************************************
// METHOD : config_gbdy
// DESC   : Configure 2d & 3d box boundary from ptree
// ARGS   : 
//          ptree   : main prop tree 
//          bterrain: is this a restart (TRUE) or not FALSE)?
//          igbdyf  : For each natural/canonical global boundary face,
//                    gives vector of global bdy ids. Allocated here.
//          igbdyft : bdy type ids for each index in igbdyf. Allocated here.
//          igbdy   : 'flat' version of igbdyf
//          degbdy  : 'node descriptor' for each point in igdby
// RETURNS: none.
//**********************************************************************************
template<typename Types>
void GGridBox<Types>::config_gbdy(const PropertyTree           &ptree, 
                           GBOOL                         bterrain,
                           GTVector<GTVector<GSIZET>>   &igbdyf, 
                           GTVector<GTVector<GBdyType>> &igbdyft,
                           GTVector<GSIZET>             &igbdy,
                           GTVector<GUINT>              &degbdy)
                          
{
  GEOFLOW_TRACE();

  // Cycle over all geometric boundaries, and configure:
  GBOOL              bret, bperiodic=FALSE;
  GINT               iret, k;
  GSIZET             igbdy_start, nind;
  GTVector<GBOOL>    buniform(2*GDIM);
  GTVector<GBdyType> bdytype(2*GDIM);
  GTVector<GUINT>    utmp;
  GTVector<GSIZET>   ikeep, itmp;
  GTVector<GString>  bdynames (2*GDIM);
  std::vector<GString>
                     svec;
  GString            gname, sbdy, bdyclass;
  PropertyTree       bdytree, gridptree;
  stBdyBlock         bcblock;
  UpdateBasePtr      base_ptr;

  assert(this->minnodedist_ > 0.0);
  this->eps_ = 0.00125*this->minnodedist_;


  this->bdyNormals_.resize(GDIM);

  // If doing re-doing with terrain, assume that the incomming 
  // bdy spec (igbdy* and degbdy) is done, and just compute
  // the normals:
  // NOTE: Later, we'll have to check that PERIODIC bdy conditions
  //       & terrain are consistent, but for now, we just assume
  //       user has specified these properly
  if  ( bterrain ) {
    for ( auto j=0; j<2*GDIM; j++ ) { // over each canonical bdy
      do_gbdy_normals(this->dXdXi_, igbdy,  degbdy, this->bdyNormals_, this->bdyTangents_); // all bdy nodes 
    }
    return;
  }



  bdynames[0] = "bdy_y_0";
  bdynames[1] = "bdy_x_1";
  bdynames[2] = "bdy_y_1";
  bdynames[3] = "bdy_x_0";
  if ( GDIM == 3 ) {
    bdynames[4] = "bdy_z_0";
    bdynames[5] = "bdy_z_1";
  }

  gname     = ptree.getValue<GString>("grid_type");
  assert(gname == "grid_box");
  gridptree = ptree.getPropertyTree(gname);


  // Clear input arrays:
  igbdyf .resize(2*GDIM);
  igbdyft.resize(2*GDIM);

  this->bdy_update_list_.resize(2*GDIM);

  // Handle uniform, nonuniform bdy conditions:
  // Note: If "uniform" not specified for a boundary, then
  //       user MUST supply a method to configure it.
  igbdy_start = 0;
  for ( auto j=0; j<2*GDIM; j++ ) { // over each canonical bdy
    sbdy         = gridptree.getValue<GString>(bdynames[j]);
    if ( "none" == sbdy || "" == sbdy ) {
      continue;
    }
    bdytree      = ptree.getPropertyTree(sbdy);
    bdyclass     = bdytree.getValue<GString>("bdy_class", "uniform");
    find_gbdy_ind  (j, FALSE, ikeep, itmp, utmp); // bdy node ids only
    igbdyf [j].resize(itmp.size()); igbdyf [j] = itmp;
    igbdyft[j].resize(itmp.size()); igbdyft[j] = GBDY_NONE;
    nind = igbdy.size(); // current size
    igbdy.reserve(nind+itmp.size());
    degbdy.reserve(nind+itmp.size());
//cout << "config_gbdy: itmp[" << j << "]=" << itmp << endl;
    for ( auto k=0; k<itmp.size(); k++ ) {
      igbdy.push_back(itmp[k]);
      degbdy.push_back(utmp[k]);
    }
    bperiodic = FALSE;
    if ( "uniform" == bdyclass ) { // uniform bdy conditions
      iret = GUpdateBdyFactory<Types>::bdy_block_conform_per(bdytree);
      if ( iret == 1 ) {
        bdytype [j] = GBDY_PERIODIC;
        igbdyft [j] = GBDY_PERIODIC;
        bperiodic   = bdytype[j] == GBDY_PERIODIC;
      }
      else if ( iret == 2 ) {
        cout << "GGridBox<Types>:: config_gbdy: Attempt to specify PERIODIC boundary failed" << endl;
        assert(FALSE);
      }

      // May have different uniform bdys for different state comps;
      // step through them in order to point to correct bdy indices:
      k = 0;
      while ( !bperiodic 
           && GUpdateBdyFactory<Types>::get_bdy_block(ptree, sbdy, k, bcblock) ) {
        bcblock.bdyid = j;
        base_ptr = GUpdateBdyFactory<Types>::build(ptree, sbdy, *this, bcblock, itmp, utmp, igbdy_start);
        igbdyft[j] = bcblock.tbdy;
        this->bdy_update_list_[j].push_back(base_ptr);
        k++;
      } 
    }
    else if ( "mixed" == bdyclass ) { // mixed bdy conditions
      cout << "GGridBox<Types>:: config_gbdy: 'mixed' bdy_class is not available"<< endl;
      assert(FALSE);
    }
    else {
      assert(FALSE && "Invalid bdy_class");
    }
    
    igbdy_start += itmp.size();

  } // end, global bdy face loop


  // With global list of domain boundaries, compute bdy data:
  do_gbdy_normals(this->dXdXi_, igbdy, degbdy, this->bdyNormals_, this->bdyTangents_); // all bdy nodes 
#if 0
GSIZET *ind=NULLPTR;
nind = 0;
igbdy.contains(0, ind, nind);
for ( auto j=0; j<nind; j++ )
cout << "GridBox::config_gbdy: igbdy[" << j << "]=" << igbdy[ind[j]] 
<< " n=(" << this->bdyNormals_[0][ind[j]] << "," << this->bdyNormals_[1][ind[j]] << ")" << endl;
delete [] ind;
#endif

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


} // end of method config_gbdy


//**********************************************************************************
//**********************************************************************************
// METHOD : is_global_vertex
// DESC   : Utilitiy method to determine of specified point is one of the 
//          global 2d boundary vertices
// ARGS   : pt : point to check
// RETURNS: TRUE on yes, else FALSE
//**********************************************************************************
template<typename Types>
GBOOL GGridBox<Types>::is_global_vertex(GTPoint<Ftype> &pt)
{
  GEOFLOW_TRACE();
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
template<typename Types>
GBOOL GGridBox<Types>::on_global_edge(GINT iface, GTPoint<Ftype> &pt)
{
  GEOFLOW_TRACE();
  assert( iface >=0 && iface <=5 && "Invalid face ID specification");

  GBOOL           bret = FALSE;
  GINT            nface=0;
  GTVector<GINT>  face(3); // at most 3 faces that point _can_ belong to
//GTPoint<Ftype> pt(ndim_);

  assert(this->minnodedist_ > 0.0);
  this->eps_ = 0.00125*this->minnodedist_;

  // Find faces point belongs to:
  for ( GINT j=0; j<ndim_; j++ ) {
    if     ( FUZZYEQ(pt.x1,P0_.x1,this->eps_) && !face.containsn(0,nface) ) 
      { face[nface] = 0; nface++; }
    else if( FUZZYEQ(pt.x2,P1_.x2,this->eps_) && !face.containsn(1,nface) ) 
      { face[nface] = 1; nface++; } 
    else if( FUZZYEQ(pt.x1,P1_.x1,this->eps_) && !face.containsn(2,nface) ) 
      { face[nface] = 2; nface++; }
    else if( FUZZYEQ(pt.x2,P0_.x2,this->eps_) && !face.containsn(3,nface) ) 
      { face[nface] = 3; nface++; }
    else if( FUZZYEQ(pt.x3,P0_.x3,this->eps_) && !face.containsn(4,nface) ) 
      { face[nface] = 4; nface++; }
    else if( FUZZYEQ(pt.x3,P1_.x3,this->eps_) && !face.containsn(5,nface) ) 
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
// METHOD : do_gbdy_normals
// DESC   : Compute normals to each domain bdy 
// ARGS   : 
//          dXdXi    : matrix of dX_i/dXi_j matrix elements, s.t.
//                     dXdX_i(i,j) = dx^i/dxi^j
//          igbdy    : vector of bdy indices into global volume fields 
//          debdy    : array of node 'descriptions'
//          normals  : vector of normal components
//          tangents : vector of tanget vectors components
// RETURNS: none
//**********************************************************************************
template<typename Types>
void GGridBox<Types>::do_gbdy_normals(const GTMatrix<GTVector<Ftype>> &dXdXi,
                               const GTVector<GSIZET>           &igbdy,
                               const GTVector<GUINT>            &debdy,
                               VVecFtype                        &normals,
                               GTVector<VVecFtype>              &tangents)
{
  GEOFLOW_TRACE();

  GSIZET nbdy;

  nbdy = igbdy.size();
  for ( auto j=0; j<normals.size(); j++ ) {
    normals[j].resize(nbdy);
    normals[j] = 0.0;
  }


  // Compute global boundary normals and associated data:

#if defined(_G_IS2D)
  do_gbdy_normals2d(dXdXi, igbdy, debdy, normals, tangents);
#elif defined(_G_IS3D)
  do_gbdy_normals3d(dXdXi, igbdy, debdy, normals, tangents);
#else
  #error Invalid problem dimensionality
#endif


} // end, method do_gbdy_normals


//**********************************************************************************
//**********************************************************************************
// METHOD : do_gbdy_normals2d
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
//          dXdXi    : matrix of dX_i/dXi_j matrix elements, s.t.
//                     dXdX_i(i,j) = dx^i/dxi^j
//          igbdy    : vector of bdy indices into global volume fields 
//          debdy    : array of node 'descriptions', with dimension of igbdy
//          normals  : vector of normal components, each of dim of igbdy
//          tangets  : vector of tangetn vector components, each of dim of igbdy
//          Note:
//          dXdXi   : matrix of dX_i/dXi_j matrix elements, s.t.
//                    dXdX_i(i,j) = dx^j/dxi^i,
//                    must be computed prior to entry.
// RETURNS: none
//**********************************************************************************
template<typename Types>
void GGridBox<Types>::do_gbdy_normals2d(const GTMatrix<GTVector<Ftype>> &dXdXi,
                                 const GTVector<GSIZET>           &igbdy,
                                 const GTVector<GUINT>            &debdy,
                                 VVecFtype                        &normals,
                                 GTVector<VVecFtype>              &tangents)
{
   GEOFLOW_TRACE();
   GINT           ib, ip;
   GUINT          id;
   Ftype          tiny;
   Ftype          xm;
   GTPoint<Ftype> kp(3), xp(3), p1(3), p2(3);

   tiny  = 100.0*std::numeric_limits<Ftype>::epsilon(); 
   kp    = 0.0;
   kp[2] = 1.0; // k-vector

   for ( auto j=0; j<normals.size(); j++ ) {
     normals[j].resize(igbdy.size());
     normals[j] = 0.0;
   }
   for ( auto j=0; j<tangents.size(); j++ ) { 
     for ( auto i=0; i<tangents[j].size(); i++ ) { 
       tangents[j][i].resize(igbdy.size());
       tangents[j][i] = 0.0;
     }
   }


   // Normals depend on element type:
   if ( this->gtype_ == GE_REGULAR ) {
     // All normal components are 0, except the one
     // perp. to edge:
     for ( auto j=0; j<igbdy.size(); j++ ) { // all points on iedge
       ib = igbdy[j];
       id = GET_NDHOST(debdy[j]); // host face id
       ip = (id+1) % 2;
       xm = id == 1 || id == 2 ? -1.0 : 1.0;
       normals[ip][j] = xm;
       tangents[0][id%2][j] = 1;
     }
   }
   else if ( this->gtype_ == GE_DEFORMED ) {
     // Bdy normal is hat{k} X dvec{X} / dxi_iedge,
     // for edge id:
     for ( auto j=0; j<igbdy.size(); j++ ) { // all points on global bdy 
       ib = igbdy[j];
       id = GET_NDHOST(debdy[j]); // host face id
       xm = id == 2 || id == 3 ? -1.0 : 1.0;
       p1 = 0.0;
       for ( auto i=0; i<dXdXi.size(2); i++ ) { // over _X_
         p1[i] = dXdXi(i,id%2)[ib]; 
       }
       kp.cross(p1, xp);   // xp = k X p1
       xp *= xm;
       xp.unit();
       for ( auto i=0; i<normals.size(); i++ ) normals[i][j] = xp[i];
       // k X tangent = n ==>
       tangents[0][0][j] =  normals[1][j]; 
       tangents[0][1][j] = -normals[0][j];
     }
   }
   else if ( this->gtype_ == GE_2DEMBEDDED ) {
     assert(FALSE && "Not available!");
     // Bdy normal is dvec{X} / dxi_xi X dvec{X} / dxi_eta
#if 0
     for ( auto j=0; j<igbdy.size(); j++ ) { // all points on global bdy
       ib = igbdy[j];
       id = GET_NDHOST(debdy[j]); // host face id
       xm = id == 1 || id == 2 ? -1.0 : 1.0;
       for ( auto i=0; i<dXdXi.size(2); i++ ) { // d_X_/dXi
         p1[i] = dXdXi(0,i)[ib]; // d_X_/dxi
         p2[i] = dXdXi(1,i)[ib]; // d_X_/deta
       }
       p1.cross(p2, xp);   // xp = p1 X p2
       xp *= xm;
       xp.unit();
       for ( ic=0; ic<xp.dim(); ic++ ) if ( fabs(xp[ic]) > tiny ) break;
       assert(ic >= GDIM); // no normal components > 0
       for ( auto i=0; i<normals.size(); i++ ) normals[i][j] = xp[i];
     }
#endif
   }

} // end, method do_gbdy_normals2d


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
//                     dXdX_i(i,j) = dx^i/dxi^j
//          normals  : vector of normal components, each of dim of igbdy
//          tangents : vector of tangent vector components, each of dim of igbdy
// RETURNS: none
//**********************************************************************************
template<typename Types>
void GGridBox<Types>::do_gbdy_normals3d(const GTMatrix<GTVector<Ftype>> &dXdXi,
                                 const GTVector<GSIZET>           &igbdy,
                                 const GTVector<GUINT>            &debdy,
                                 VVecFtype                        &normals,
                                 GTVector<VVecFtype>              &tangents)
{
   GEOFLOW_TRACE();
   GSIZET          ib, ip; 
   GUINT           id;
   GINT            ixi[6][2] = { {0,2}, {1,2}, {0,2}, 
                                 {1,2}, {0,1}, {0,1} };
   Ftype          tiny;
   Ftype          xm;
   GTPoint<Ftype> xp(3), p1(3), p2(3);
   tiny  = 100.0*std::numeric_limits<Ftype>::epsilon(); 

   // There must be 2 tangent vectors defining tangent plane:
   assert(tangents.size() == 2 ); 

   for ( auto j=0; j<normals.size(); j++ ) {
     normals[j].resize(igbdy.size());
     normals[j] = 0.0;
   }
   for ( auto j=0; j<tangents.size(); j++ ) { 
     for ( auto i=0; i<tangents[j].size(); i++ ) { 
       tangents[j][i].resize(igbdy.size());
       tangents[j][i] = 0.0;
     }
   }

   // Normals depend on element type:
   if ( this->gtype_ == GE_REGULAR ) {
     // All normal components are 0, except the one
     // perp. to face:
     for ( auto j=0; j<igbdy.size(); j++ ) { // all points on face
       id = GET_NDHOST(debdy[j]); // host face id
       ip = id < 4 ? (id+1)%2 : 2; // perp component
       xm = id == 1 || id == 2 || id == 5 ? 1.0 : -1.0;
       for ( auto i=0; i<normals.size(); i++ ) normals[i][j] = 0.0; 
       normals[ip][j] = xm;
       for ( auto i=0; i<2; i++ ) tangents[i][(GSIZET)ixi[i]][j] = 1.0;
     }
   }
   else if ( this->gtype_ == GE_DEFORMED ) {
     // Bdy normal is dvec{X} / dxi_xi X dvec{X} / dxi_eta
     for ( auto j=0; j<igbdy.size(); j++ ) { // all points on iedge
       ib = igbdy[j];
       xm = id == 1 || id == 2 || id == 5 ? 1.0 : -1.0;
       for ( auto i=0; i<dXdXi.size(2); i++ ) { // over _X_
         p1[i] = dXdXi(i,ixi[id][0])[ib]; // d_X_/dxi
         p2[i] = dXdXi(i,ixi[id][1])[ib]; // d_X_/deta
       }
       p1.cross(p2, xp);   // xp = p1 X p2
       xp.unit(); 
       for ( auto i=0; i<normals.size(); i++ ) normals[i][j] = xp[i];

       // Use Gram-Schmidt orthogonalization on p1 & p2, and use
       // the orthogonal vectors for tangent space. Let p1 be one 
       // vector; create the other by subtracting from p2 the component 
       // of p2 parallel to p1: xp = p2 - (p1.p2)/p2^2 p1:
       xp = p1; xp *= -p2.dot(p1)/pow(p1.mag(),2); xp += p2;
       p1.unit(); xp.unit();
       for ( auto i=0; i<GDIM; i++ ) tangents[0][i][j] = p1[i];
       for ( auto i=0; i<GDIM; i++ ) tangents[1][i][j] = xp[i];
     }
   }

} // end, method do_gbdy_normals3d



//**********************************************************************************
//**********************************************************************************
// METHOD : find_gbdy_ind
// DESC   : Find global bdy indices (indices into xNodes_ arrays) that
//          corresp to specified bdy. This is the main entry point.
// ARGS   : bdyid    : box-boundary id
//          bunique  : don't include repeated vertex/edge points on 
//                     bdy face that are shared with another bdy face
//          ikeep    : used if bunique=TRUE; will store current list of bdy
//                     node ids for all bdy nodes that don't already appear in 
//                     ikeep buffer
//          debdy    : array of descriptors for each index in ibd
//          ibdy     : array of indices into xNodes that comprise this boundary
// RETURNS: none.
//**********************************************************************************
template<typename Types>
void GGridBox<Types>::find_gbdy_ind(GINT bdyid, GBOOL bunique, 
                             GTVector<GSIZET> &ikeep,
                             GTVector<GSIZET> &ibdy,
                             GTVector<GUINT>  &debdy) 
{
  GEOFLOW_TRACE();

#if defined(_G_IS2D)
  find_gbdy_ind2d(bdyid, bunique, ikeep, ibdy, debdy);
#elif defined(_G_IS3D)
  find_gbdy_ind2d(bdyid, bunique, ikeep, ibdy, debdy);
#else
  #error "Dimensionality undefined"
#endif
  
} // end, method find_gbdy_ind 


//**********************************************************************************
//**********************************************************************************
// METHOD : find_gbdy_ind2d
// DESC   : Find global bdy indices (indices into xNodes_ arrays) that
//          corresp to specified bdy in 2d
// ARGS   : bdyid    : box-boundary id
//          bunique  : don't include repeated vertex/edge points on 
//                     bdy face that are shared with another bdy face
//          ikeep    : used if bunique=TRUE; will store current list of bdy
//                     node ids for all bdy nodes that don't already appear in 
//                     ikeep buffer
//          debdy    : array of descriptors for each index in ibd
//          ibdy     : array of indices into xNodes that comprise this boundary
// RETURNS: none.
//**********************************************************************************
template<typename Types>
void GGridBox<Types>::find_gbdy_ind2d(GINT bdyid, GBOOL bunique, 
                               GTVector<GSIZET> &ikeep,
                               GTVector<GSIZET> &ibdy,
                               GTVector<GUINT>  &debdy) 
{
  GEOFLOW_TRACE();
  assert(bdyid >=0 && bdyid < 2*GDIM);

  GBOOL             bexist, bglobale, bglobalv;
  GUINT             if2v[2*GDIM][2] = { {0,1}, {1,2}, {3,2}, {0,3} }; // vertex ids for each face
  GUINT             iv2f[2*GDIM][2] = { {0,3}, {0,1}, {1,2}, {2,3} }; // face ids for each vertex
  GUINT             ivert;
  GSIZET            ie, istart, nbdy, nind, nkeep, ntmp, nnodes;
  GTVector<GUINT>   utmp;
  GTVector<GSIZET>  ind, itmp, ktmp;
  GTPoint<Ftype>   xp;
  GTVector<GTPoint<Ftype>>
                    v(2);
  GTVector<GTVector<Ftype>>
                    *xlnodes;
  
  nkeep = ikeep.size();

  assert(this->minnodedist_ > 0.0);
  this->eps_ = 0.00125*this->minnodedist_;
  
  ntmp  = bunique ? this->xNodes_[0].size()
        : this->xNodes_[0].size() + pow(2,GDIM) + (GDIM > 2 ? 2*GDIM : 0);
  itmp.resize(ntmp);
  if ( bunique ) ktmp.resize(ntmp);
  utmp.resize(ntmp); utmp = 0;

  for ( auto j=0; j<nkeep; j++ ) ktmp[j] = ikeep[j];

  v[0].setBracket(this->eps_);
  v[1].setBracket(this->eps_);
  xp  .setBracket(this->eps_);
  switch ( bdyid ) { // get vertices defining bdy:
    // NOTE: the order of these vertices should
    //       corresp with those in if2v array:
    case 0: // lower horiz bdy:
      v[0][0] = P0_[0]; v[0][1] = P0_[1];
      v[1][0] = P1_[0]; v[1][1] = P0_[1];
      break;
    case 1: // east vert bdy:
      v[0][0] = P1_[0]; v[0][1] = P0_[1];
      v[1][0] = P1_[0]; v[1][1] = P1_[1];
      break;
    case 2: // top horiz bdy:
      v[0][0] = P0_[0]; v[0][1] = P1_[1];
      v[1][0] = P1_[0]; v[1][1] = P1_[1];
      break;
    case 3: // west vert bdy:
      v[0][0] = P0_[0]; v[0][1] = P0_[1];
      v[1][0] = P0_[0]; v[1][1] = P1_[1];
      break;
    default:
      assert(FALSE);
  }

   
  istart = 0;
  nbdy   = 0;
  for ( auto e=0; e<this->gelems_.size(); e++ ) { 
    xlnodes= &this->gelems_[e]->xNodes();
    nnodes = this->gelems_[e]->nnodes();
    nind   = geoflow::in_seg<Ftype>(v, *xlnodes, this->eps_, ind); // get indices on segment
    // For each index on bdy, set description:
    for ( auto i=0; i<nind; i++ ) { 
      xp.assign(*xlnodes, ind[i]);
      ie      = ind[i] + istart ;
      if ( bunique ) bexist = ktmp.containsn(ie,nkeep);
      if ( !bunique || !bexist ) {
        // Check if point is a global vertex; set node description:
        bglobalv = FALSE;
        if      ( xp == v[0] ) { ivert = if2v[bdyid][0]; bglobalv = TRUE; }
        else if ( xp == v[1] ) { ivert = if2v[bdyid][1]; bglobalv = TRUE; }
        if ( bglobalv ) { // elem vertex on global  bdy vertex
//cout << "config_gbdy: bdyid=" << bdyid << " ie=" << ie << " gblobalv=" << bglobalv << " xp=" << xp << " v0=" << v[0] << " `v1=" << v[1] << endl;
          SET_ND(utmp[nbdy], ivert, GElem_base::VERTEX, bdyid);
          itmp [nbdy] = ie; // 'global' index 
          if ( bunique ) ktmp[nkeep] = ie; // to store in keep
          nbdy++; 
          if ( bunique ) nkeep++;
        }
        else { // is on bdy face
          SET_ND(utmp[nbdy], bdyid, GElem_base::FACE, bdyid);
          itmp [nbdy] = ie; // 'global' index 
           if ( bunique ) ktmp[nkeep] = ie; // to store in keep
          nbdy++; 
          if ( bunique ) nkeep++;
        }
      }
    } // end, element's glob bdy
    istart += nnodes;
  } // end, elem loop


  // Fill return arrays:
  ibdy.resize(nbdy);
  debdy.resize(nbdy);
  if ( bunique ) ikeep.resize(nkeep);

  for( auto j=0; j<nbdy; j++ ) {
    ibdy [j] = itmp[j];  // bdy node in volume 
    debdy[j] = utmp[j];  // bdy node description
  }
  for( auto j=0; j<bunique && nkeep; j++ ) ikeep[j] = ktmp[j];

} // end, method find_gbdy_ind2d


//**********************************************************************************
//**********************************************************************************
// METHOD : find_gbdy_ind3d
// DESC   : Find global bdy indices (indices into xNodes_ arrays) that
//          corresp to specified bdy in 3d
// ARGS   : bdyid    : box-boundary id
//          bunique  : don't include repeated vertex/edge points on 
//                     bdy face that are shared with another bdy face
//          ikeep    : used if bunique=TRUE; will store current list of bdy
//                     node ids for all bdy nodes that don't already appear in 
//                     ikeep buffer
//          debdy    : array of descriptors for each index in ibd
//          ibdy     : array of indices into xNodes that comprise this boundary
// RETURNS: none.
//**********************************************************************************
template<typename Types>
void GGridBox<Types>::find_gbdy_ind3d(GINT bdyid, GBOOL bunique, 
                               GTVector<GSIZET> &ikeep,
                               GTVector<GSIZET> &ibdy,
                               GTVector<GUINT>  &debdy) 
{
  GEOFLOW_TRACE();
  assert(bdyid >=0 && bdyid < 2*GDIM);

  GBOOL             bexist, bglobale, bglobalv;
  GUINT             ev[4][2]; 
  GUINT             ivert;
  GSIZET            ie, istart, nbdy, nind, nind1, nkeep, ntmp, nnodes;
  GTVector<GUINT>   utmp;
  GTVector<GSIZET>  ind, ind1, itmp, ktmp;
  GTPoint<Ftype>   xp;
  GTVector<GTPoint<Ftype>>
                    e(2), v(4);
  GTVector<GTVector<Ftype>>
                   *xlnodes;
  
  nkeep = ikeep.size();

  assert(this->minnodedist_ > 0.0);
  this->eps_ = 0.00125*this->minnodedist_;
  
  ntmp  = bunique ? this->xNodes_[0].size()
        : this->xNodes_[0].size() + pow(2,GDIM) + (GDIM > 2 ? 2*GDIM : 0);
  itmp.resize(ntmp);
  if ( bunique ) ktmp.resize(ntmp);
  utmp.resize(ntmp); utmp = 0;

  for ( auto j=0; j<bunique && nkeep; j++ ) ktmp[j] = ikeep[j];

  v[0].setBracket(this->eps_);
  v[1].setBracket(this->eps_);
  xp  .setBracket(this->eps_);
  
  switch ( bdyid ) { // get vertices defining bdy:
    // NOTE: the order of these vertices should
    //       corresp with those in iv array:
    case 0: // face id
      v [0][0] = P0_[0]; v [0][1] = P0_[1]; v [0][2] = P0_[2]; // vertices
      v [1][0] = P1_[0]; v [1][1] = P0_[1]; v [1][2] = P0_[2];
      v [2][0] = P1_[0]; v [2][1] = P0_[1]; v [2][2] = P1_[2];
      v [3][0] = P0_[0]; v [3][1] = P0_[1]; v [3][2] = P1_[2];
      ev[0][0] = 0     ; ev[0][1] = 1; // edge vertices ids
      ev[1][0] = 1     ; ev[1][1] = 2;
      ev[2][0] = 3     ; ev[2][1] = 2;
      ev[3][0] = 0     ; ev[3][1] = 3;
      break;
    case 1: 
      v [0][0] = P1_[0]; v [0][1] = P0_[1]; v [0][2] = P0_[2]; 
      v [1][0] = P1_[0]; v [1][1] = P1_[1]; v [1][2] = P0_[2];
      v [2][0] = P1_[0]; v [2][1] = P1_[1]; v [2][2] = P1_[2];
      v [3][0] = P1_[0]; v [3][1] = P0_[1]; v [3][2] = P0_[2];
      ev[0][0] = 0; ev[0][1] = 1;
      ev[1][0] = 1; ev[1][1] = 2;
      ev[2][0] = 3; ev[2][1] = 2;
      ev[3][0] = 0; ev[3][1] = 3;
      break;
    case 2: 
      v [0][0] = P0_[0]; v [0][1] = P1_[1]; v [0][2] = P0_[2]; 
      v [1][0] = P1_[0]; v [1][1] = P1_[1]; v [1][2] = P0_[2];
      v [2][0] = P1_[0]; v [2][1] = P1_[1]; v [2][2] = P1_[2];
      v [3][0] = P0_[0]; v [3][1] = P1_[1]; v [3][2] = P1_[2];
      ev[0][0] = 0; ev[0][1] = 1;
      ev[1][0] = 1; ev[1][1] = 2;
      ev[2][0] = 3; ev[2][1] = 2; 
      ev[3][0] = 0; ev[3][1] = 3;
      break;
    case 3:
      v [0][0] = P0_[0]; v [0][1] = P0_[1]; v [0][2] = P0_[2]; 
      v [1][0] = P0_[0]; v [1][1] = P1_[1]; v [1][2] = P0_[2];
      v [2][0] = P0_[0]; v [2][1] = P1_[1]; v [2][2] = P1_[2];
      v [3][0] = P0_[0]; v [3][1] = P0_[1]; v [3][2] = P1_[2];
      ev[0][0] = 0; ev[0][1] = 1; 
      ev[1][0] = 1; ev[1][1] = 2; 
      ev[2][0] = 3; ev[2][1] = 2;
      ev[3][0] = 0; ev[3][1] = 3;
      break;
    case 4: 
      v [0][0] = P0_[0]; v [0][1] = P0_[1]; v [0][2] = P0_[2]; 
      v [1][0] = P1_[0]; v [1][1] = P0_[1]; v [1][2] = P0_[2];
      v [2][0] = P1_[0]; v [2][1] = P1_[1]; v [2][2] = P0_[2];
      v [3][0] = P0_[0]; v [3][1] = P1_[1]; v [3][2] = P0_[2];
      ev[0][0] = 0; ev[0][1] = 1;
      ev[1][0] = 1; ev[1][1] = 2;
      ev[2][0] = 3; ev[2][1] = 2;
      ev[3][0] = 0; ev[3][1] = 3;
      break;
    case 5: 
      v [0][0] = P0_[0]; v [0][1] = P0_[1]; v [0][2] = P1_[2]; 
      v [1][0] = P1_[0]; v [1][1] = P0_[1]; v [1][2] = P1_[2];
      v [2][0] = P1_[0]; v [2][1] = P1_[1]; v [2][2] = P1_[2];
      v [3][0] = P0_[0]; v [3][1] = P1_[1]; v [3][2] = P1_[2];
      ev[0][0] = 0; ev[0][1] = 1;
      ev[1][0] = 1; ev[1][1] = 2;
      ev[2][0] = 3; ev[2][1] = 2;
      ev[3][0] = 0; ev[3][1] = 3;
      break;
    default:
      assert(FALSE);
  }

   
  istart = 0;
  nbdy   = 0;
  for ( auto e=0; e<this->gelems_.size(); e++ ) { 
    xlnodes= &this->gelems_[e]->xNodes();
    nnodes = this->gelems_[e]->nnodes();
    nind   = geoflow::in_poly<Ftype>(v, *xlnodes, this->eps_, ind); // get indices on domain surface
    // For each index on bdy, set description:
    for ( auto i=0; i<nind; i++ ) { 
      xp.assign(*xlnodes, ind[i]);
      ie      = ind[i] + istart ;
      if ( bunique ) bexist = ktmp.containsn(ie,nkeep);
      if ( !bunique || !bexist ) {
        // Do not check for node identity, just uniquenes:
        SET_ND(utmp[nbdy], bdyid, GElem_base::FACE, bdyid);
        itmp [nbdy] = ie; // 'global' index 
        if ( bunique ) ktmp[nkeep] = ie; // to store in keep
        nbdy++; 
        if ( bunique ) nkeep++;
      } // end, unique conditional
    } // end, element's glob bdy
    istart += nnodes;
  } // end, elem loop


  // Fill return arrays:
  ibdy.resize(nbdy);
  debdy.resize(nbdy);
  if ( bunique ) ikeep.resize(nkeep);

  for( auto j=0; j<nbdy; j++ ) {
    ibdy [j] = itmp[j];  // bdy node in volume 
    debdy[j] = utmp[j];  // bdy node description
  }
  for( auto j=0; j<bunique && nkeep; j++ ) ikeep[j] = ktmp[j];

} // end, method find_gbdy_ind3d


//**********************************************************************************
//**********************************************************************************
// METHOD : elem_face_data
// DESC   : Main entry point into methods that compute ids, normals, 
//          and mass x Jacobian for each face of a 2d or 3d element
//          
// ARGS   : dXdXi     : matrix of dX_i/dXi_j matrix elements, s.t.
//                      dXdX_i(i,j) = dx^j/dxi^i, specified 'globally'
//          gieface   : index of elem bdys in global volume, set here
//          face_mass : mass*Jac at face nodes; should contain
//                      Gauss weights on element faces on entry, and
//                      contain the weights at each face node
//          normals   : vector of normal components
// RETURNS: none.
//**********************************************************************************
template<typename Types>
void GGridBox<Types>::elem_face_data(GTMatrix<GTVector<Ftype>> &dXdXi,
                              GTVector<GSIZET>                 &gieface,
                              GTVector<Ftype>                 &face_mass,
                              GTVector<GTVector<Ftype>>       &normals)
{
    
#if defined(_G_IS2D)
  elem_face_data2d(dXdXi, gieface, face_mass, normals);
#elif defined (_G_IS3D)
  elem_face_data3d(dXdXi, gieface, face_mass, normals);
#else
  #error "Dimensionality undefined"
#endif

} // end. elem_face_data main entry point


//**********************************************************************************
//**********************************************************************************
// METHOD : elem_face_data2d
// DESC   : Find 2d element boundary normals and mass
// ARGS   : dXdXi     : matrix of dX_i/dXi_j matrix elements, s.t.
//                      dXdX_i(i,j) = dx^j/dxi^i, specified 'globally'
//          gieface   : index of elem bdys in global volume, set here
//          face_mass : mass*Jac at face nodes; should contain
//                      Gauss weights on element faces on entry, and
//                      contain the weights at each face node
//          normals   : vector of normal components at each elem bdy node
// RETURNS: none.
//**********************************************************************************
template<typename Types>
void GGridBox<Types>::elem_face_data2d(GTMatrix<GTVector<Ftype>>       &dXdXi,
                                GTVector<GSIZET>                 &gieface,
                                GTVector<Ftype>                 &face_mass,
                                GTVector<GTVector<Ftype>>       &normals)
{
   GEOFLOW_TRACE();

   GINT              fi, ib, id, ip;
   GSIZET            istart, nbdy;
   Ftype            tiny;
   Ftype            jac, xm;
   GTPoint<Ftype>   kp(3), xp(3), p1(3), p2(3);
   GTVector<GINT>   *face_ind;
   GTVector<Ftype> *mass, *pdX;

   tiny  = 100.0*std::numeric_limits<Ftype>::epsilon();
   kp    = 0.0;
   kp[2] = 1.0; // k-vector

   
   // Get number of elem bdy nodes, and 
   // allocate arrays:
   nbdy = 0;
   for ( auto e=0; e<this->gelems_.size(); e++ ) {
     for ( auto j=0; j<this->gelems_[e]->nfaces(); j++ ) {
       nbdy += this->gelems_[e]->face_indices(j).size();
     }
   }
   gieface.resize(nbdy);
   face_mass  .resize(nbdy);
   for ( auto j=0; j<normals.size(); j++ ) {
     normals[j].resize(nbdy);
     normals[j] = 0.0;
   }
   

   // Normals, data depend on element type:
   nbdy = 0;
   if ( this->gtype_ == GE_REGULAR ) {
    
     istart = 0;
     for ( auto e=0; e<this->gelems_.size(); e++ ) { // over all elements
       mass      = &this->gelems_[e]->face_mass();
       for ( auto j=0; j<this->gelems_[e]->nfaces(); j++ ) { // loop over faces
         face_ind  = &this->gelems_[e]->face_indices(j);
         id = j % 2;
         ip = (j+1) % 2; pdX = &dXdXi(id,0);
         xm = j == 1 || j == 2 ? 1.0 : -1.0;
         for ( auto i=0; i<face_ind->size(); i++ ) { // over elem bdy points
           fi = (*face_ind)[i]; // index into elem data
           ib = fi + istart;    // index into global volume data
           gieface      [nbdy] = ib;
           face_mass    [nbdy] = (*mass)[fi] * (*pdX)[ib];
           normals[ip][nbdy++] = xm;
//cout << "GBox::elem_face_data2d: nbdy=" << nbdy << " e=" << e << " face=" << j << " nx=" << normals[0][nbdy-1] << " ny=" << normals[1][nbdy-1] << " faceMass=" << face_mass[nbdy-1] << " gie=" << ib << endl;
         }
       } // end, elem face loop
       istart += this->gelems_[e]->nnodes();
     } // end, elem loop

   }
   else if ( this->gtype_ == GE_DEFORMED   ) {

     // Bdy normal is hat{k} X dvec{X} / dxi_iedge,
     // for face:
     istart = 0;
     for ( auto e=0; e<this->gelems_.size(); e++ ) { // over all elements
       mass      = &this->gelems_[e]->face_mass();
       for ( auto j=0; j<this->gelems_[e]->nfaces(); j++ ) { 
         face_ind  = &this->gelems_[e]->face_indices(j);
         id = j % 2;
         xm = j == 1 || j == 2 ? 1.0 : -1.0;
         for ( auto i=0; i<face_ind->size(); i++ ) { // over elem bdy points
           fi = (*face_ind)[i]; // index into elem data
           ib = fi + istart;    // index into global volume data
           for ( auto k=0; k<dXdXi.size(2); k++ ) p1[k] = dXdXi(id,k)[ib];
           kp.cross(p1, xp);   // xp = k X p1 == elem face normal
           xp *= xm; xp.unit();
           gieface      [nbdy] = ib;
           for ( auto k=0; k<normals.size(); k++ ) normals[k][nbdy] = xp[k];
           jac = p1.mag();
           face_mass [nbdy++] = (*mass)[fi] * jac; 
         }
       } // end, elem face loop
       istart += this->gelems_[e]->nnodes();
     } // end, elem loop

   }
   else {
     assert(FALSE);
   }

} // end, method elem_face_data2d


//**********************************************************************************
//**********************************************************************************
// METHOD : elem_face_data3d
// DESC   : Find 3d element boundary normals and mass
// ARGS   : dXdXi     : matrix of dX_i/dXi_j matrix elements, s.t.
//                      dXdX_i(i,j) = dx^j/dxi^i, specified 'globally'
//          gieface   : index of elem bdys in global volume, set here
//          face_mass : mass*Jac at face nodes; should contain
//                      Gauss weights on element faces on entry, and
//                      contain the weights at each face node
//          normals   : vector of normal components at each elem bdy node
// RETURNS: none.
//**********************************************************************************
template<typename Types>
void GGridBox<Types>::elem_face_data3d(GTMatrix<GTVector<Ftype>>       &dXdXi,
                                GTVector<GSIZET>                 &gieface,
                                GTVector<Ftype>                 &face_mass,
                                GTVector<GTVector<Ftype>>       &normals)
{
   GEOFLOW_TRACE();

   // Reference coords for each face:
   GUINT             if2r[6][2] = { {0,2}, {1,2}, {0,2}, 
                                    {1,2}, {0,1}, {0,1} }; 
   // normal direction for each rergular face:
   GUINT             if2n[6]    = { 0, 1, 0, 1, 2, 2 };
   GINT              fi, ib, ip;
   GSIZET            istart, nbdy;
   Ftype            tiny;
   Ftype            jac, xm;
   GTPoint<Ftype>   kp(3), xp(3), p1(3), p2(3);
   GTVector<GINT>   *face_ind;
   GTVector<Ftype> *mass;

   tiny  = 100.0*std::numeric_limits<Ftype>::epsilon();
   kp    = 0.0;
   kp[2] = 1.0; // k-vector

   
   // Get number of elem bdy nodes, and 
   // allocate arrays:
   nbdy = 0;
   for ( auto e=0; e<this->gelems_.size(); e++ ) {
     for ( auto j=0; j<this->gelems_[e]->nfaces(); j++ ) {
       nbdy += this->gelems_[e]->face_indices(j).size();
     }
   }
   gieface  .resize(nbdy);
   face_mass.resize(nbdy);
   for ( auto j=0; j<normals.size(); j++ ) {
     normals[j].resize(nbdy);
     normals[j] = 0.0;
   }

   // Normals, data  depend on element type:
   nbdy = 0;
   if ( this->gtype_ == GE_REGULAR ) {
    
     istart = 0;
     for ( auto e=0; e<this->gelems_.size(); e++ ) { // over all elements
       mass      = &this->gelems_[e]->face_mass();
       for ( auto j=0; j<this->gelems_[e]->nfaces(); j++ ) { 
         face_ind  = &this->gelems_[e]->face_indices(j);
         ip = if2n[j];
         xm = j < 4 || j == 5 ? 1.0 : -1.0;
         for ( auto i=0; i<face_ind->size(); i++ ) { // over elem bdy points
           fi = (*face_ind)[i]; // index into elem data
           ib = fi + istart;    // index into global volume data
           gieface    [nbdy] = ib;
           jac               = dXdXi(if2r[j][0],0)[ib] * dXdXi(if2r[j][0],1)[ib];
           normals[ip][nbdy] = xm;
           face_mass[nbdy++] = (*mass)[fi] * jac;
         }
       } // end, elem face loop
       istart += this->gelems_[e]->nnodes();
     } // end, elem loop

   }
   else if ( this->gtype_ == GE_DEFORMED   ) {

     // Bdy normal is hat{k} X dvec{X} / dxi_iedge,
     // for face:
     istart = 0;
     for ( auto e=0; e<this->gelems_.size(); e++ ) { // over all elements
       mass      = &this->gelems_[e]->face_mass();
       for ( auto j=0; j<this->gelems_[e]->nfaces(); j++ ) { 
         face_ind  = &this->gelems_[e]->face_indices(j);
         xm = j < 4 || j == 5 ? 1.0 : -1.0;
         for ( auto i=0; i<face_ind->size(); i++ ) { // over elem bdy points
           fi = (*face_ind)[i]; // index into elem data
           ib = fi + istart;    // index into global volume data
           for ( auto k=0; k<dXdXi.size(2); k++ ) p1[k] = dXdXi(if2r[j][0],k)[ib];
           for ( auto k=0; k<dXdXi.size(2); k++ ) p2[k] = dXdXi(if2r[j][1],k)[ib];
           p1.cross(p2,xp); // normal = p1 X p2
           gieface    [nbdy] = ib;
           // faceJac = (dX/dxi X dX/deta) . dX/dzeta, with (xi, eta) in 
           // plane; zeta the last reference coord:
           jac               = xp.x1 * dXdXi(2,0)[ib] 
                             + xp.x2 * dXdXi(2,1)[ib]
                             + xp.x2 * dXdXi(2,2)[ib];
           xp *= xm; xp.unit();
           for ( auto k=0; k<normals.size(); k++ ) normals[k][nbdy] = xp[k];
           face_mass[nbdy++] = (*mass)[fi] * jac; 
         }
       } // end, elem face loop
       istart += this->gelems_[e]->nnodes();
     } // end, elem loop

   }
   else {
     assert(FALSE);
   }

} // end, method elem_face_data3d


