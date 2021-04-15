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
template<typename Types> 
GGridIcos<Types>::GGridIcos(const geoflow::tbox::PropertyTree &ptree, GTVector<GNBasis<GCTYPE,typename Types::Ftype>*> &b, GC_COMM &comm)
:          GGrid<Types>(ptree, b, comm),
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
  
  GString gname   = ptree.getValue<GString>("grid_type");
  assert(gname == "grid_icos" || gname == "grid_sphere");
  geoflow::tbox::PropertyTree gridptree = ptree.getPropertyTree(gname);

  gbasis_.resize(b.size());
  gbasis_ = b;
  lshapefcn_ = new GShapeFcn_linear<GTICOS>(2);
  ilevel_  = gridptree.getValue<GINT>("ilevel");
  sreftype_= gridptree.getValue<GString>("refine_type","GICOS_LAGRANGIAN");

  
  if ( ndim_ == 2 ) {
    assert(GDIM == 2 && "GDIM must be 2");
    radiusi_ = gridptree.getValue<Ftype>("radius");
    init2d();
  }
  else if ( ndim_ == 3 ) {
    assert(GDIM == 3 && "GDIM must be 3");
    std::vector<GINT> ne(3);
    radiusi_   = gridptree.getValue<Ftype>("radiusi");
    radiuso_   = gridptree.getValue<Ftype>("radiuso");
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
template<typename Types> 
GGridIcos<Types>::~GGridIcos()
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
std::ostream &operator<<(std::ostream &str, GGridIcos<Types> &e)
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
#endif


//**********************************************************************************
//**********************************************************************************
// METHOD : set_partitioner
// DESC   : Set domain decomposition object
// ARGS   : GDD_base pointer
// RETURNS: none
//**********************************************************************************
template<typename Types> 
void GGridIcos<Types>::set_partitioner(GDD_base<GTICOS> *gdd)
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
template<typename Types> 
void GGridIcos<Types>::init2d()
{
	GEOFLOW_TRACE();
  GString serr = "GridIcos::init2d: ";

  Ftype phi = (1.0+sqrt(5.0))/2.0;  // Golden ratio


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
  Ftype fact = 1.0/sqrt(phi*phi + 1.0);
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
template<typename Types> 
void GGridIcos<Types>::init3d()
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
template<typename Types> 
void GGridIcos<Types>::lagrefine()
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
template<typename Types> 
void GGridIcos<Types>::do_elems()
{
  GEOFLOW_TRACE();
  if ( ndim_ == 2 ) do_elems2d(this->irank_);
  if ( ndim_ == 3 ) do_elems3d(this->irank_);

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
void GGridIcos<Types>::do_elems(GTMatrix<GINT> &p,
                        GTVector<GTVector<Ftype>> &xnodes)
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
template<typename Types> 
void GGridIcos<Types>::do_elems2d(GINT irank)
{
	GEOFLOW_TRACE();
  GString           serr = "GridIcos::do_elems2d (1): ";
  Ftype             fact;
  GTVector<GTPoint<GTICOS>> cverts(4), gverts(4), tverts(4);
  GTPoint<GTICOS>    ct(3), v1(3), v2(3), v3(3); // 3d points
  GElem_base        *pelem;
  GTVector<GINT>    iind;
  GTVector<GINT>    I(1);
  GTVector<GTICOS>  Ni;
  GTVector<GTVector<Ftype>>   *xNodes;
  GTVector<GTVector<Ftype>*>  *xiNodes;
  GTVector<GTVector<GTICOS>>    xid;
  GTVector<GTVector<GTICOS>*>   pxid;
  GTVector<GTVector<GTICOS>>    xgtmp(3);

  // Do eveything on unit sphere, then project to radiusi_
  // as a final step.

  assert(gbasis_.size()>0 && "Must set basis first");

  if ( gdd_ == NULLPTR ) gdd_ = new GDD_base<GTICOS>(this->nprocs_);

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
      copycast<Ftype,GTICOS>(*xiNodes, pxid);
      
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

      this->gelems_.push_back(pelem);

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
template<typename Types> 
void GGridIcos<Types>::do_elems3d(GINT irank)
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
  GTVector<GTVector<Ftype>>   *xNodes;
  GTVector<GTVector<Ftype>>    xNodes2d(2);
  GTVector<GTVector<Ftype>*>   xiNodes2d(2);
  GTVector<Ftype>             *xiNodesr;
  GTVector<GTVector<GTICOS>>   xd, xd2d;
  GTVector<GTVector<GTICOS>>   xid, xid2d;
  GTVector<GTVector<GTICOS>*>  pxid;
  GTVector<GTVector<GTICOS>>   xgtmp(3);

  // Do eveything on unit sphere, then project to radiusi_
  // as a final step.

  assert(gbasis_.size()>0 && "Must set basis first");

  if ( gdd_ == NULLPTR ) gdd_ = new GDD_base<GTICOS>(this->nprocs_);

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
      copycast<Ftype,GTICOS>(xiNodes2d, pxid);

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
        spherical2xyz<Ftype>(*xNodes);

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

        this->gelems_.push_back(pelem);
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
template<typename Types> 
void GGridIcos<Types>::do_elems2d(GTMatrix<GINT> &p,
                           GTVector<GTVector<Ftype>> &gxnodes)
{
	GEOFLOW_TRACE();
  GString                     serr = "GridIcos::do_elems2d (2): ";
  GElem_base                  *pelem;
  GTVector<GTVector<Ftype>>  *xNodes;
  GTVector<GNBasis<GCTYPE,Ftype>*>
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
    this->gelems_.push_back(pelem);

    assert(nvnodes == this->gelems_[i]->nnodes() && "Incompatible node count");
    nfnodes = this->gelems_[i]->nfnodes();
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
template<typename Types> 
void GGridIcos<Types>::do_elems3d(GTMatrix<GINT> &p,
                           GTVector<GTVector<Ftype>> &gxnodes)
{
	GEOFLOW_TRACE();
  GString                      serr = "GridIcos::do_elems3d (2): ";
  GElem_base                  *pelem;
  GTVector<GTVector<Ftype>>   *xNodes;
  GTVector<GTVector<Ftype>*>  *xiNodes;
  GTVector<GNBasis<GCTYPE,Ftype>*>
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
    this->gelems_.push_back(pelem);

    assert(nvnodes == this->gelems_[i]->nnodes() && "Incompatible node count");
    nfnodes = this->gelems_[i]->nfnodes();
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
template<typename Types> 
void GGridIcos<Types>::print(const GString &filename, GCOORDSYST icoord)
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
void GGridIcos<Types>::config_gbdy(const geoflow::tbox::PropertyTree &ptree, 
                            GBOOL                         bterrain,
                            GTVector<GTVector<GSIZET>>   &igbdyf, 
                            GTVector<GTVector<GBdyType>> &igbdyft,
                            GTVector<GSIZET>             &igbdy,
                            GTVector<GUINT>              &degbdy)
{
	GEOFLOW_TRACE();
  // Cycle over all geometric boundaries, and configure:

  GBOOL              bret;
  GINT               iret, k;
  GSIZET             igbdy_start, nind;
  GTVector<GUINT>    utmp;
  GTVector<GSIZET>   ikeep, itmp;
  GTVector<Ftype>    rbdy(2);
  GTVector<GString>  bdynames(2);
  std::vector<GString>  svec;
  GString            gname, sbdy, bdyclass;
  PropertyTree       bdytree, gridptree;
  stBdyBlock         bcblock;
  UpdateBasePtr      base_ptr;


  this->bdyNormals_.resizem(GDIM);

  // If doing re-doing with terrain, assume that the incomming 
  // bdy spec (igbdy* and degbdy) is done, and just compute
  // the normals:
  // NOTE: Later, we'll have to check that PERIODIC bdy conditions
  //       & terrain are consistent, but for now, we just assume
  //       user has specified these properly
  if  ( bterrain ) {
    for ( auto j=0; j<2*GDIM; j++ ) { // over each canonical bdy
      do_gbdy_normals(this->dXdXi_, igbdy, degbdy, this->bdyNormals_, this->idepComp_); // all bdy nodes 
    }
    return;
  }




  // Clear input arrays:
  igbdyf .clear();
  igbdyft.clear();

  if ( ndim_ == 2 ) return; // no boundaries to configure
 
  bdynames[0] = "bdy_inner";
  bdynames[1] = "bdy_outer";

  if ( !ptree.isValue<GString>("grid_type") ) {
    cout << "GGridIcos<Types>::config_gbdy: grid_type not set" << endl;
    assert(FALSE);
  }
  gname     = ptree.getValue<GString>("grid_type");

  if ( !ptree.isPropertyTree(gname) ) {
    cout << "GGridIcos<Types>::config_gbdy: grid_type block " << gname << " not found" << endl;
    assert(FALSE);
  }
  gridptree = ptree.getPropertyTree(gname);

  rbdy[0] = radiusi_;
  rbdy[1] = radiuso_;

  igbdyf.resize(2); // 2 canonical bdys
  igbdyft.resize(2); // 2 canonical bdys

  this->bdy_update_list_.resize(2);
 

  // Get properties from the main prop tree. 
  // Note: bdys are configured by way of geometry's
  //       natural decomposition: here, by inner and
  //       outer spherical surfaces. But the bdy indices 
  //       and types returned on exit contain info for all bdys:
  igbdy_start = 0;
  for ( auto j=0; j<2; j++ ) { // cycle over 2 spherical surfaces
    sbdy         = gridptree.getValue<GString>(bdynames[j]);
    bdytree      = ptree.getPropertyTree(sbdy);
    bdyclass     = bdytree.getValue<GString>("bdy_class", "uniform");
    find_gbdy_ind3d(rbdy[j], itmp, utmp); // bdy node ids only

    igbdyf [j].resize(itmp.size()); igbdyf [j] = itmp;
    igbdyft[j].resize(itmp.size()); igbdyft[j] = GBDY_NONE;
    nind = 0;
    for ( auto k=0; k<igbdyf.size(); k++ ) nind += igbdyf[k].size();
    igbdy  .resize(nind); // vol indices of bdy nodes in base; bdy update needs this
    degbdy .resize(nind); 
    nind = 0;
    for ( auto i=0; i<igbdyf[j].size(); i++ ) {
      igbdy[nind] = igbdyf[j][i];
      degbdy[nind++] = igbdyf[j][i];
    }

    if ( "uniform" == bdyclass ) { // uniform bdy conditions
      iret = GUpdateBdyFactory<Types>::bdy_block_conform_per(bdytree);
      if ( iret == 1 || iret == 2 ) {
        cout << "GGridIcos<Types>:: config_gbdy: PERIODIC boundary conditions invalid for this grid" << endl;
        assert(FALSE);
      }
 
      // May have different uniform bdys for different state comps;
      // step through them in order to point to correct bdy indices:
      k = 0;
      while ( GUpdateBdyFactory<Types>::get_bdy_block(ptree, sbdy, k, bcblock) ) {;
        bcblock.bdyid = j;
        base_ptr = GUpdateBdyFactory<Types>::build(ptree, sbdy, *this, bcblock, itmp, igbdy_start);
        igbdyft[j] = bcblock.tbdy;
        this->bdy_update_list_[j].push_back(base_ptr);
        k++;
      }
    }
    else if ( "mixed" == bdyclass ) { // mixed bdy conditions
      cout << "GGridIcos<Types>:: config_gbdy: 'mixed' bdy_class is not available"<< endl;
      assert(FALSE);
    }
    else {
      assert(FALSE && "Invalid bdy_class");
    }

    igbdy_start += itmp.size();

  } // end, canonical bdy loop

  // With global list of domain boundaries, compute bdy data:
  do_gbdy_normals(this->dXdXi_, igbdy, degbdy, this->bdyNormals_, this->idepComp_); 

} // end of method config_gbdy



//**********************************************************************************
//**********************************************************************************
// METHOD : find_gbdy_ind3d
// DESC   : Find global bdy indices (indices into xNodes_ arrays) that
//          corresp to specified bdy in 3d
// ARGS   : radius   : booundary radius
//          debdy    : array of descriptors for each index in ibd
//          ibdy     : array of indices into xNodes that comprise this boundary
// RETURNS: none.
//**********************************************************************************
template<typename Types> 
void GGridIcos<Types>::find_gbdy_ind3d(Ftype radius,
                                GTVector<GSIZET> &ibdy,
                                GTVector<GUINT>  &debdy) 
{
  GEOFLOW_TRACE();

  GBOOL             bexist, bglobale, bglobalv;
  GSIZET            ie, istart, nbdy, nind, ntmp, nnodes;
  GTVector<GUINT>   utmp;
  GTVector<GINT>   *face_ind;
  GTVector<GSIZET>  ind, itmp;
  GTVector<GSIZET>  fi;
  GTVector<Ftype>   r;
  GTVector<GTVector<Ftype>>
                   *xlnodes;
  
  ntmp  = this->xNodes_[0].size() + pow(2,GDIM) + (GDIM > 2 ? 2*GDIM : 0);
  itmp.resize(ntmp);
  utmp.resize(ntmp); utmp = 0;

  this->eps_ = 0.125*this->minnodedist_;

  
  istart = 0;
  nbdy   = 0;
  for ( auto e=0; e<this->gelems_.size(); e++ ) { 
    xlnodes= &this->gelems_[e]->xNodes();
    nnodes = this->gelems_[e]->nnodes();
    r.resize(nnodes);
    for ( auto j=0; j<nnodes; j++ ) r[j] = sqrt(pow((*xlnodes)[0][j],2)
                                         +      pow((*xlnodes)[1][j],2)
                                         +      pow((*xlnodes)[2][j],2));
    for ( auto j=0; j<this->gelems_[e]->nfaces(); j++ ) {
       face_ind  = &this->gelems_[e]->face_indices(j);
       geoflow::convert<GINT,GSIZET>(*face_ind, fi);
       nind      = geoflow::fuzzyeq<Ftype>(r, fi, radius, this->eps_, ind); // get indices on domain surface
      // For each index on bdy, set description:
      for ( auto i=0; i<nind; i++ ) { 
        ie      = ind[i] + istart ;
        // Set node host face: 
        SET_ND(utmp[nbdy], 0, GElem_base::FACE, (GUINT)j );
        itmp [nbdy] = ie; // 'global' index 
        nbdy++; 
      } // end, element's glob bdy
    } // end, j-loop
    istart += nnodes;
  } // end, elem loop


  // Fill return arrays:
  ibdy.resize(nbdy);
  debdy.resize(nbdy);

  for( auto j=0; j<nbdy; j++ ) {
    ibdy [j] = itmp[j];  // bdy node in volume 
    debdy[j] = utmp[j];  // bdy node description
  }

} // end, method find_gbdy_ind3d



//**********************************************************************************
//**********************************************************************************
// METHOD : elem_face_data
// DESC   : Compute normals to each element face
// ARGS   : 
//          dXdXi     : matrix of dX_i/dXi_j matrix elements, s.t.
//                      dXdX_i(i,j) = dx^j/dxi^i
//          gieface   : vector of face indices into global volume fields 
//                      for all facase
//          gdeface   : description for each face node
//          face_mass : mass on eleme faces
//          normals   : vector of normal components
// RETURNS: none
//**********************************************************************************
template<typename Types> 
void GGridIcos<Types>::elem_face_data(GTMatrix<GTVector<Ftype>> &dXdXi,
                               GTVector<GSIZET>                 &gieface,
                               GTVector<Ftype>                  &face_mass,
                               GTVector<GTVector<Ftype>>        &normals)
{

	GEOFLOW_TRACE();
  #if defined(_G_IS2D)
    elem_face_data2d(dXdXi, gieface, face_mass, normals);
  #elif defined(_G_IS3D)
    elem_face_data3d(dXdXi, gieface, face_mass, normals);
  #else
    #error Invalid problem dimensionality
  #endif

} // end, method elem_face_data


//**********************************************************************************
//**********************************************************************************
// METHOD : elem_face_data2d
// DESC   : Compute normals to each element face in 2d
// ARGS   : 
//          dXdXi     : matrix of dX_i/dXi_j matrix elements, s.t.
//                      dXdX_i(i,j) = dx^j/dxi^i
//          gieface   : vector of face indices into global volume fields 
//                      for all facase
//          face_mass : mass on eleme faces
//          normals   : vector of normal components at each elem bdy node
// RETURNS: none
//**********************************************************************************
template<typename Types> 
void GGridIcos<Types>::elem_face_data2d(GTMatrix<GTVector<Ftype>> &dXdXi,
                                GTVector<GSIZET>                  &gieface,
                                GTVector<Ftype>                   &face_mass,
                                GTVector<GTVector<Ftype>>         &normals)
{
   GEOFLOW_TRACE();
   GINT              fi, ib, ip;
   GSIZET            istart, nbdy;
   Ftype             tiny;
   Ftype             jac, xm;
   GTPoint<Ftype>    kp(3), xp(3), p1(3), p2(3);
   GTVector<GINT>   *face_ind;
   GTVector<Ftype>  *mass;


   tiny  = 100.0*std::numeric_limits<Ftype>::epsilon();
   kp    = 1.0/sqrt(3.0); // unit vector in radial direction

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

   nbdy = 0;
   if ( this->gtype_ == GE_DEFORMED
    ||  this->gtype_ == GE_2DEMBEDDED ) {

     istart = 0;
     for ( auto e=0; e<this->gelems_.size(); e++ ) {
       for ( auto j=0; j<this->gelems_[e]->nfaces(); j++ ) { // over all faces
         face_ind  = &this->gelems_[e]->face_indices(j);
         mass      = &this->gelems_[e]->face_mass();
         xm = j == 2 || j == 3 ? 1.0 : -1.0;
         for ( auto i=0; i<face_ind->size(); i++ ) { // over elem bdy points
           fi = (*face_ind)[i]; // index into elem data
           ib = fi + istart;    // indir. index into global vol
           for ( auto i=0; i<dXdXi.size(2); i++ ) { // over _X_
             p1[i] = dXdXi(j%2,i)[ib];
           }
           kp.cross(p1, xp);   // xp =  p1 X k
           gieface[nbdy] = ib;
           xp *= xm; xp.unit();
           for ( auto k=0; k<normals.size(); k++ ) normals[k][nbdy] = xp[k];
           jac = p1.mag();
           face_mass [nbdy++] = (*mass)[fi] * jac;
         }
       } // end, elem face loop
       istart += this->gelems_[e]->nnodes();
     } // end, elem loop
   }
   else {
     assert(FALSE && "Invalid grid type");
   }

} // end, method elem_face_data2d


//**********************************************************************************
//**********************************************************************************
// METHOD : elem_face_data3d
// DESC   : Compute normals to each element face in 3d
// ARGS   : 
//          dXdXi     : matrix of dX_i/dXi_j matrix elements, s.t.
//                      dXdX_i(i,j) = dx^j/dxi^i
//          gieface   : vector of face indices into global volume fields 
//                      for all facase
//          face_mass : mass on eleme faces
//          normals   : vector of normal components at each elem bdy node
// RETURNS: none
//**********************************************************************************
template<typename Types> 
void GGridIcos<Types>::elem_face_data3d(GTMatrix<GTVector<Ftype>> &dXdXi,
                                 GTVector<GSIZET>                 &gieface,
                                 GTVector<Ftype>                  &face_mass,
                                 GTVector<GTVector<Ftype>>        &normals)
{
   GEOFLOW_TRACE();

   GINT             fi, ib, ic, id;
   GINT             ixi[6][2] = { {0,2}, {1,2}, {2,0},
                                  {2,1}, {1,0}, {0,1} };
   GSIZET           istart, nbdy;
   Ftype            jac, tiny, xm;
   GTPoint<Ftype>   xp(3), p1(3), p2(3);
   GTVector<GINT>   *face_ind;
   GTVector<Ftype>  *mass;

   tiny  = 100.0*std::numeric_limits<Ftype>::epsilon();

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

   nbdy = 0;
   if ( this->gtype_ == GE_DEFORMED
    ||  this->gtype_ == GE_2DEMBEDDED ) {

     istart = 0;
     for ( auto e=0; e<this->gelems_.size(); e++ ) {
       mass      = &this->gelems_[e]->face_mass();
       for ( auto j=0; j<this->gelems_[e]->nfaces(); j++ ) { // over all faces
         face_ind  = &this->gelems_[e]->face_indices(j);
         xm = j == 0 || j == 1 || j == 5 ? 1.0 : -1.0;
         // Bdy normal is dvec{X} / dxi_xi X dvec{X} / dxi_eta
         for ( auto i=0; i<face_ind->size(); i++ ) { // over elem bdy points
           fi = (*face_ind)[i]; // index into elem data
           ib = fi + istart;    // index into global volume data
           // Find derivs of _X_ wrt face's reference coords;
           // the cross prod of these vectors is the normal:
           for ( auto i=0; i<dXdXi.size(2); i++ ) { // d_X_/dXi
             p1[i] = dXdXi(ixi[j][0],i)[ib]; // d_X_/dxi
             p2[i] = dXdXi(ixi[j][1],i)[ib]; // d_X_/deta
           }
           p1.cross(p2, xp);   // xp = p1 X p2
           xp *= xm; xp.unit(); 
           gieface      [nbdy] = ib;
           for ( auto k=0; k<normals.size(); k++ ) normals[k][nbdy] = xp[k];
           jac = xp.mag();
           face_mass  [nbdy++] = (*mass)[fi] * jac;
         } // end, bdy points loop, i
       } // end, elem face loop, j
       istart += this->gelems_[e]->nnodes();
     } // end, elem loop, e
   }
   else {
     assert(FALSE && "Invalid grid type");
   }

} // end, method elem_face_data3d


//**********************************************************************************
//**********************************************************************************
// METHOD : do_gbdy_normals
// DESC   : Compute normals to each domain bdy 
// ARGS   : 
//          dXdXi    : matrix of dX_i/dXi_j matrix elements, s.t.
//                     dXdX_i(i,j) = dx^j/dxi^i
//          igbdy    : vector of bdy indices into global volume fields 
//          debdy    : array of node 'descriptions', with dimension of igbdy
//          normals  : vector of normal components, each of dim of igbdy
//          idepComp : vector index dependent on the other indices (first 
//                     component index whose normal component is nonzero)

// RETURNS: none
//**********************************************************************************
template<typename Types> 
void GGridIcos<Types>::do_gbdy_normals(const GTMatrix<GTVector<Ftype>> &dXdXi,
                                const GTVector<GSIZET>                 &igbdy,
                                const GTVector<GUINT>                  &debdy,
                                GTVector<GTVector<Ftype>>              &normals,
                                GTVector<GINT>                         &idepComp)
{
  GEOFLOW_TRACE();
  GSIZET nbdy;

  #if defined(_G_IS2D)
    return;
  #elif defined(_G_IS3D)
    nbdy = igbdy.size();
    idepComp.resize(nbdy);
    for ( auto j=0; j<normals.size(); j++ ) normals[j].resize(nbdy);

    do_gbdy_normals3d(dXdXi, igbdy, debdy, normals, idepComp);

  #else
    #error "Invalid problem dimensionality"
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
//          idepComp : vector index dependent on the other indices (first 
//                     component index whose normal component is nonzero)

// RETURNS: none
//**********************************************************************************
template<typename Types> 
void GGridIcos<Types>::do_gbdy_normals3d(const GTMatrix<GTVector<Ftype>> &dXdXi,
                                  const GTVector<GSIZET>                 &igbdy,
                                  const GTVector<GUINT>                  &debdy,
                                  GTVector<GTVector<Ftype>>              &normals,
                                  GTVector<GINT>                         &idepComp)
{
  GEOFLOW_TRACE();
  GSIZET         ib, ic, ip;
  GUINT          id;
  Ftype          tiny;
  Ftype          xm;
  GTPoint<Ftype> xp(3), p1(3), p2(3);

  tiny  = 100.0*std::numeric_limits<Ftype>::epsilon();

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
       idepComp[j] = ic;  // dependent component
     }
   }
   else {
     assert(FALSE && "Invalid grid type");
   }

} // end, method do_gbdy_normals3d


//**********************************************************************************
//**********************************************************************************
// METHOD : project2sphere (1)
// DESC   : Project Cartesian coords to sphere, specified by rad argument, and
//          express in Cartesian coords.  Necessary for 2d grids.
// ARGS   : tmesh: Vector of triangles (vertices), modified to contain 
//                 projected coordinates
// RETURNS: none.
//**********************************************************************************
template<typename Types>
template<typename T>
void GGridIcos<Types>::project2sphere(GTVector<GTriangle<T>> &tmesh, T rad)
{
	GEOFLOW_TRACE();
  GString serr = "GridIcos::project2sphere (1): ";

  T          r, xlat, xlong;
  GTPoint<T> v;

  for ( GSIZET i=0; i<tmesh_.size(); i++ ) { // loop over all triangles in tmesh
    for ( GSIZET j=0; j<3; j++ ) { // guaranteed to have 3 points each
      v = *tmesh[i].v[j];
      r = v.norm();
      xlat  = asin(v.x3/r);
      xlong = atan2(v.x2,v.x1);
      xlong = xlong < 0.0 ? 2.0*PI+xlong : xlong;
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
template<typename Types>
template<typename T>
void GGridIcos<Types>::project2sphere(GTVector<GTPoint<T>> &plist, T rad)
{
	GEOFLOW_TRACE();
  GString serr = "GridIcos::project2sphere (2): ";

  T          r, xlat, xlong;
  GTPoint<T> v;

  for ( GSIZET i=0; i<plist.size(); i++ ) { // loop over all points
    v = plist[i];
    r = v.norm();
    xlat  = asin(v.x3/r);
    xlong = atan2(v.x2,v.x1);
    xlong = xlong < 0.0 ? 2.0*PI+xlong : xlong;
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
template<typename Types>
template<typename T>
void GGridIcos<Types>::project2sphere(GTVector<GTVector<T>> &plist, T rad)
{
	GEOFLOW_TRACE();
  GString serr = "GridIcos::project2sphere (3): ";

  T r, xlat, xlong, x, y, z;

  for ( GSIZET i=0; i<plist[0].size(); i++ ) { // loop over all points
    x = plist[0][i]; y = plist[1][i]; z = plist[2][i];
    r = sqrt(x*x + y*y + z*z);
    xlat  = asin(z/r);
    xlong = atan2(y,x);
    xlong = xlong < 0.0 ? 2.0*PI+xlong : xlong;
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
template<typename Types>
template<typename T>
void GGridIcos<Types>::spherical2xyz(GTVector<GTPoint<T>*> &plist)
{
	GEOFLOW_TRACE();
  GString serr = "GridIcos::spherical2xyz(1): ";

  T r, xlat, xlong;

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
template<typename Types>
template<typename T>
void GGridIcos<Types>::spherical2xyz(GTVector<GTPoint<T>> &plist)
{
	GEOFLOW_TRACE();
  GString serr = "GridIcos::spherical2xyz(2): ";

  T r, xlat, xlong;

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
template<typename Types>
template<typename T>
void GGridIcos<Types>::spherical2xyz(GTVector<GTVector<T>> &plist)
{
	GEOFLOW_TRACE();
  GString serr = "GridIcos::spherical2xyz(3): ";
  assert(plist.size() >= 3 && "Invalid dimensionality on input array");

  T r, xlat, xlong;

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
template<typename Types>
template<typename T>
void GGridIcos<Types>::xyz2spherical(GTVector<GTPoint<T>*> &plist)
{
	GEOFLOW_TRACE();
  GString serr = "GridIcos::xyz2spherica(1): ";

  T r, x, y, z;

  for ( GSIZET i=0; i<plist.size(); i++ ) { // loop over all points
   x = plist[i]->x1; y = plist[i]->x2; z = plist[i]->x3;
   r = sqrt(x*x + y*y + z*z);
   plist[i]->x1 = r;
   plist[i]->x2 = asin(z/r);
   plist[i]->x3 = atan2(y,x);
   plist[i]->x3 = plist[i]->x3 < 0.0 ? 2.0*PI+plist[i]->x3 : plist[i]->x3;
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
template<typename Types>
template<typename T>
void GGridIcos<Types>::xyz2spherical(GTVector<GTVector<T>> &plist)
{
	GEOFLOW_TRACE();
  GString serr = "GridIcos::xyz2spherica(2): ";

  T r, x, y, z;

  for ( GSIZET i=0; i<plist.size(); i++ ) { // loop over all points
   x = plist[i][0]; y = plist[i][1]; z = plist[i][2];
   r = sqrt(x*x + y*y + z*z);
   plist[i][0] = r;
   plist[i][1] = asin(z/r);
   plist[i][2] = atan2(y,x);
   plist[i][2] = plist[i][2] < 0.0 ? 2.0*PI+plist[i][2] : plist[i][2];
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
template<typename Types>
template<typename T>
void GGridIcos<Types>::xyz2spherical(GTVector<GTPoint<T>> &plist)
{
	GEOFLOW_TRACE();
  GString serr = "GridIcos::xyz2spherica(3): ";

  T r, x, y, z;

  for ( GSIZET i=0; i<plist.size(); i++ ) { // loop over all points
   x = plist[i].x1; y = plist[i].x2; z = plist[i].x3;
   r = sqrt(x*x + y*y + z*z);
   plist[i].x1 = r;
   plist[i].x2 = asin(z/r);
   plist[i].x3 = atan2(y,x);
   plist[i].x3 = plist[i].x3 < 0.0 ? 2.0*PI+plist[i].x3 : plist[i].x3;
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
template<typename Types>
template<typename T>
void GGridIcos<Types>::cart2gnomonic(GTVector<GTPoint<T>> &clist, T rad, T xlatc, T xlongc, GTVector<GTPoint<T>> &glist)
{
	GEOFLOW_TRACE();
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

  T xlatp, xlongp;
  T den, r, xlat, xlong;
  T eps = 10.0*std::numeric_limits<GFTYPE>::epsilon();
  for ( GSIZET i=0; i<clist.size(); i++ ) { // loop over all points
    r      = clist[i].norm();
    xlat   = asin(clist[i].x3/r);
    xlong  = atan2(clist[i].x2,clist[i].x1);
    xlong  = xlong < 0.0 ? 2.0*PI+xlong : xlong;
#if 0
    den    = sin(xlatc)*sin(xlat) + cos(xlatc)*cos(xlat)*cos(xlong-xlongc);  
    xlatp  = asin( cos(xlatc)*sin(xlat) - sin(xlatc)*cos(xlat)*cos(xlong-xlongc) );
    xlongp = atan2( cos(xlat)*sin(xlong-xlongc),den );

    glist[i].x1 = rad*tan(xlongp);
    glist[i].x2 = rad*tan(xlatp)*sec(xlongp);
#else
    den    = sin(xlatc)*sin(xlat) + cos(xlatc)*cos(xlat)*cos(xlong-xlongc); 
    den    = fabs(den) < eps ? 0.0 : 1.0/den;
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
template<typename Types>
template<typename T>
void GGridIcos<Types>::gnomonic2cart(GTVector<GTVector<T>> &glist, T rad, T xlatc, T xlongc, GTVector<GTVector<T>> &clist)
{
	GEOFLOW_TRACE();
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

  T X, Y;
  T beta, rho, sign, x, xlat, xlong, y;
  T eps = 10.0*std::numeric_limits<GFTYPE>::epsilon();
  for ( GSIZET i=0; i<glist[0].size(); i++ ) { // loop over all points
    x      = glist[0][i];
    y      = glist[1][i];


    if ( fabs(x) < eps && fabs(y) < eps ) {
      xlat = xlatc;
      xlong = xlongc;
    } else {
      rho    = sqrt(x*x + y*y);
      beta   = atan(rho);
      Y      = cos(beta)*sin(xlatc) + (y*sin(beta)*cos(xlatc))/rho;
      sign   = copysign(1.0, Y);
      Y      = sign* MIN(fabs(Y),1.0);
      xlat   = asin(Y);
      Y      = x*sin(beta);
      X      = rho*cos(xlatc)*cos(beta)-y*sin(xlatc)*sin(beta);
      xlong  = xlongc + atan2(Y, X);
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
template<typename Types>
template<typename T>
void GGridIcos<Types>::gnomonic2cart(GTVector<GTPoint<T>> &glist, T rad, T xlatc, T xlongc, GTVector<GTPoint<T>> &clist)
{
	GEOFLOW_TRACE();
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

  T X, Y;
  T beta, rho, sign, x, xlat, xlong, y;
  T eps = 10.0*std::numeric_limits<GFTYPE>::epsilon();
  for ( GSIZET i=0; i<clist.size(); i++ ) { // loop over all points
    x      = glist[i][0];
    y      = glist[i][1];

    if ( fabs(x) < eps && fabs(y) < eps ) {
      xlat = xlatc;
      xlong = xlongc;
    } else {
      rho    = sqrt(x*x + y*y);
      beta   = atan(rho);
      Y      = cos(beta)*sin(xlatc) + (y*sin(beta)*cos(xlatc))/rho;
      sign   = copysign(1.0, Y);
      Y      = sign* MIN(fabs(Y),1.0);
      xlat   = asin(Y);
      Y      = x*sin(beta);
      X      = rho*cos(xlatc)*cos(beta)-y*sin(xlatc)*sin(beta);
      xlong  = xlongc + atan2(Y, X);
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
// DESC   : Reorder specified elem vertices to be consistent with
//          shape functions
// ARGS   : uverts : list of unordered 3-vertices
//          tverts : temp array of points equal in number to uverts
//          isort  : array of sorting indirection indices
//          overts : array of ordered 3-vertices, returned
// RETURNS: none.
//**********************************************************************************
template<typename Types>
template<typename T>
void GGridIcos<Types>::reorderverts2d(GTVector<GTPoint<T>> &uverts, GTVector<GTPoint<T>> &tverts, GTVector<GSIZET> &isort, 
                               GTVector<GTPoint<T>> &overts)
{
	GEOFLOW_TRACE();
  GString serr = "GridIcos::reorderverts2d: ";

  assert(uverts.size() == 4 
      && overts.size() == 4
      && "Incorrect number of vertices");

  T                xlatc, xlongc;
  GTPoint<T>       c(3);
  GTVector<T>      x(4);
  GTVector<GSIZET> Ixy(4);
  
  for ( auto j=0; j<tverts.size(); j++ ) tverts[j].resize(2);

  c = (uverts[0] + uverts[1] + uverts[2] + uverts[3])*0.25; // elem centroid
  xlatc  = asin(c.x3); // reference lat/long
  xlongc = atan2(c.x2,c.x1);
  xlongc = xlongc < 0.0 ? 2*PI+xlongc : xlongc;

  // Convert to gnomonic coords so that we
  // can examine truly 2d planar coords when
  // re-ordering:
  cart2gnomonic<T>(uverts, 1.0, xlatc, xlongc, tverts); // gnomonic vertices of quads

  isort.resize(4);
  for ( auto i=0; i<tverts.size(); i++ ) { 
    x[i] = tverts[i].x1;
  }

  x.sortincreasing(Ixy);

  // Do 'left' -hand vertices:
  if ( tverts[Ixy[0]].x2 < tverts[Ixy[1]].x2 ) {
    isort[0] = Ixy[0];
    isort[3] = Ixy[1];
  } else {
    isort[0] = Ixy[1];
    isort[3] = Ixy[0];
  }

  // Do 'right' -hand vertices:
  if ( tverts[Ixy[2]].x2 < tverts[Ixy[3]].x2 ) {
    isort[1] = Ixy[2];
    isort[2] = Ixy[3];
  } else {
    isort[1] = Ixy[3];
    isort[2] = Ixy[2];
  }

  for ( auto j=0; j<4; j++ ) overts[j] = uverts[isort[j]];
  
} // end of method reorderverts2d



//**********************************************************************************
//**********************************************************************************
// METHOD : copycast (1)
// DESC   : copy 'from' array to 'to' array, while casting. Also, make sure
//          sizes are correct.
// ARGS   : from : 'from' array
//          to   : 'to' array; may be of different type than 'from'
// RETURNS: none.
//**********************************************************************************
template<typename Types>
template<typename TF, typename TT>
void GGridIcos<Types>::copycast(GTVector<GTVector<TF>> &from, GTVector<GTVector<TT>> &to)
{
	GEOFLOW_TRACE();
  assert(to.size() == from.size() && "Incompatible dimensions");
  for ( auto j=0; j<to.size(); j++ ) {
    to[j].resize(from[j].size()); 
    for ( auto i=0; i<to[j].size(); i++ ) {
      to[j][i] = static_cast<TT>(from[j][i]);
    }
  }
  
} // end of method copycast (1)



//**********************************************************************************
//**********************************************************************************
// METHOD : copycast (2)
// DESC   : copy 'from' array to 'to' array, while casting. Also, make sure
//          sizes are correct.
// ARGS   : from : 'from' array
//          to   : 'to' array; may be of different type than 'from'
// RETURNS: none.
//**********************************************************************************
template<typename Types>
template<typename TF, typename TT>
void GGridIcos<Types>::copycast(GTVector<GTVector<TF>*> &from, GTVector<GTVector<TT>*> &to)
{
	GEOFLOW_TRACE();
  assert(to.size() == from.size() && "Incompatible dimensions");
  for ( auto j=0; j<to.size(); j++ ) {
    to[j]->resize(from[j]->size()); 
    for ( auto i=0; i<to[j]->size(); i++ ) {
      (*to[j])[i] = static_cast<TT>((*from[j])[i]);
    }
  }
  
} // end of method copycast (2)



//**********************************************************************************
//**********************************************************************************
// METHOD : copycast (3)
// DESC   : copy 'form' point to 'to' point, while casting
// ARGS   : from : 'from' point
//          to   : 'to' point; may be of different type than 'from'
// RETURNS: none.
//**********************************************************************************
template<typename Types>
template<typename TF, typename TT>
void GGridIcos<Types>::copycast(GTPoint<TF> &from, GTPoint<TT> &to)
{
	GEOFLOW_TRACE();
  for ( auto j=0; j<to.size(); j++ ) {
    to[j] = static_cast<TT>(from[j]);
  }
  
} // end of method copycast (3)



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
template<typename Types>
template<typename T>
void GGridIcos<Types>::interleave(GTVector<GTPoint<T>> &R0, GTVector<GTPoint<T>> &R1,
                   GINT I, GTVector<GTPoint<T>> &Rz)
{
	GEOFLOW_TRACE();
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
template<typename Types>
template<typename T>
void GGridIcos<Types>::lagvert(GTPoint<T>&a, GTPoint<T> &b, GTPoint<T> &c,
                   GINT I, GTVector<GTPoint<T>> &R)
{
	GEOFLOW_TRACE();
  GString serr = "GridIcos::lagvert: ";

  T   xI, xJ;

  GTPoint<T> rL(3);
  GTPoint<T> rR(3);

  R.resizem(I+1);

  T fact = 1.0/static_cast<T>(nrows_+1);

  // Build 'rail' points on L and R:
  xI = static_cast<T>(I);
  rL = a + ( (b - a) * (xI * fact) );
  rR = a + ( (c - a) * (xI * fact) );

  // Compute R vertices based on refinement indices:
  fact = I > 0 ? 1.0/static_cast<T>(I) : 1.0;
  for ( GSIZET j=0; j<I+1; j++ ) {
    xJ = static_cast<T>(j);
    R[j] = rL + (rR - rL)*(xJ*fact);
  }

} // end of method lagvert


//**********************************************************************************
//**********************************************************************************
// METHOD : order_latlong2d
// DESC   : Order 2d vertices on exit s.t. they roughly define a 'box' 
//          in spherical coords. Used only for quad elements.
// ARGS   : verts : Array of vertices, re-ordered on exit
// RETURNS: none.
//**********************************************************************************
template<typename Types>
template<typename T>
void GGridIcos<Types>::order_latlong2d(GTVector<GTPoint<T>> &verts)
{
	GEOFLOW_TRACE();
  assert(verts.size() == 4 && "4 vertices must be provided");

  GString              serr = "GGridIcos<Types>::order_latlong2d: ";
  GTVector<GSIZET>     isortlon(4);
  GTVector<T>  lon(4), lat(4);
  GTVector<GTPoint<T>> cverts(4);       // copy of input verts
  GTVector<GTPoint<T>> sverts(4);       // sverts in sph coords

  cverts = verts;
  sverts = verts;
  xyz2spherical<T>(sverts); // convert verts to latlon

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
template<typename Types>
template<typename T>
void GGridIcos<Types>::order_triangles(GTVector<GTriangle<T>> &tmesh)
{
	GEOFLOW_TRACE();
  GString              serr = "GGridIcos<Types>::order_triangles: ";
  GTVector<GSIZET>     isortlon(3);
  GTVector<T>          lon(3), lat(3);
  GTVector<GTPoint<T>> cverts(3);       // copy of input verts
  GTVector<GTPoint<T>> sverts(3);       // sverts in sph coords

  iup_.resize(tmesh.size());
  iup_ = 0;
  for ( GSIZET i=0; i<tmesh.size(); i++ ) {
    for ( GSIZET j=0; j<3; j++ ) cverts[j] = *tmesh[i].v[j];
    sverts = cverts;
    xyz2spherical<T>(sverts); // convert verts to latlon
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


