//==================================================================================
// Module       : ggrid_box
// Date         : 11/11/18 (DLR)
// Description  : Object defining a (global) 2d or 3d box grid.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include "ggrid_box.hpp"



//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Instantiate with box lengths (x, y, z), directional basis, and
//          domain decompoisition object. This generates a 2d or 3d grid
//          depending on whether b.size = 2 or 3..
// ARGS   : ptree  : proptery tree
//          b      : vector of basis pointers. Number of elements in b is
//                   what _must_ be in L, and ne.
//          comm   : communicator
// RETURNS: none
//**********************************************************************************
GGridBox::GGridBox(geoflow::tbox::PropertyTree &ptree, GTVector<GNBasis<GCTYPE,GFTYPE>*> &b, GC_COMM comm)
:
ndim_                     (GDIM),
comm_                     (comm),
nprocs_ (GComm::WorldSize(comm)),
gdd_                   (NULLPTR),
lshapefcn_             (NULLPTR),
bdycallback_           (NULLPTR)
{
  assert((b.size() == GDIM ) 
        && "Basis has incorrect dimensionalilty");

  gbasis_.resize(GDIM);
  gbasis_ = b;
  Lbox_.resize(GDIM);
  ne_.resize(GDIM);

  std::vector<GFTYPE> pt = ptree.getArray<GFTYPE>("xyz0");
  P0_ = pt;
  std::vector<GFTYPE> pt = ptree.getArray<GFTYPE>("delxyz");
  GTPoint<GFTYPE> dP(3);
  dP = pt;
  P1_ = P0_ + dP;

  GTVector<GBdyType> bdytype.resize(2*GDIM);
  bdytype[0] = geoflow::str2bdytype(ptree.getValue<GString>("bdy_y_0"));
  bdytype[1] = geoflow::str2bdytype(ptree.getValue<GString>("bdy_x_1"));
  bdytype[2] = geoflow::str2bdytype(ptree.getValue<GString>("bdy_y_1"));
  bdytype[3] = geoflow::str2bdytype(ptree.getValue<GString>("bdy_x_0"));
  if ( GDIM == 3 ) {
  bdytype[4] = geoflow::str2bdytype(ptree.getValue<GString>("bdy_z_0"));
  bdytype[5] = geoflow::str2bdytype(ptree.getValue<GString>("bdy_z_1"));
  }
  global_bdy_types_ = bdytype;

  ne         = ptree.getArray<GINT>("num_elems"));
  ne_        = ne;

  bPeriodic_.resize(3);
  bPeriodic_ = FALSE;
  assert(P0_.size() >= ndim_ && P1_.size() >= ndim_ && ne.size() >= b.size()
        && "Grid length and/or element count arrays of insufficient size");
  if ( global_bdy_types_[1] == GBDY_PERIODIC
    || global_bdy_types_[3] == GBDY_PERIODIC ) bPeriodic_[0] = TRUE;
  if ( global_bdy_types_[0] == GBDY_PERIODIC
    || global_bdy_types_[2] == GBDY_PERIODIC ) bPeriodic_[1] = TRUE;
  
  if ( GDIM==3 
    && global_bdy_types_[1] == GBDY_PERIODIC
    || global_bdy_types_[3] == GBDY_PERIODIC ) bPeriodic_[2] = TRUE;

  for ( GSIZET j=0; j<b.size(); j++ ) {
    Lbox_[j] = fabs(P1_[j] - P0_[j]);
    ne_  [j] = ne[j];
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
// METHOD : set_basis 
// DESC   : Set basis object for coord. direction, idir in (1, 2, 3)
// ARGS   :
// RETURNS: none
//**********************************************************************************
void GGridBox::set_basis(GTVector<GNBasis<GCTYPE,GFTYPE>*> &b)
{
  GString serr = "GGridBox::set_basis: ";

  gbasis_.resize(b.size());;
  gbasis_ = b;
} // end of method set_basis


//**********************************************************************************
//**********************************************************************************
// METHOD : init2d
// DESC   : Initialize base state/base icosahedron
// ARGS   : none.
// RETURNS: none.
//**********************************************************************************
void GGridBox::init2d()
{
  GString serr = "GGridBox::init2d: ";

  assert(gbasis_.size() == 2 && "Invalid dimension");

  GINT nelems = ne_[0]*ne_[1];
  qmesh_.resize(nelems);


  // Get uniform element sizes:
  GTPoint<GFTYPE> dx(ndim_);
  for ( GSIZET j=0; j<ndim_; j++ ) 
    dx[j] = Lbox_[j] / static_cast<GFTYPE>(ne_[j]);

  // Compute vertices of all hexes ('cubes'), starting 
  GINT   n=0;
  GTPoint<GFTYPE> v0(2);     // starting sph. point of elem
  GTPoint<GFTYPE> dv(2);     // delta-point
  for ( GSIZET j=0; j<ne_[1]; j++ ) { 
    for ( GSIZET i=0; i<ne_[0]; i++ ) { 
      v0.x1 = i*dx.x1; v0.x2 = j*dx.x2; 
      qmesh_[n].v1 = v0;
      dv.x1 = dx.x1 ; dv.x2 = 0.0  ;   qmesh_[n].v2 = v0 + dv;
      dv.x1 = dx.x1 ; dv.x2 = dx.x2;   qmesh_[n].v3 = v0 + dv;
      dv.x1 = 0.0   ; dv.x2 = dx.x2;   qmesh_[n].v4 = v0 + dv;
      n++;
    }
  }
  
  // Compute centroids of all hexes ('cubes'):
  GTPoint<GFTYPE> a(2);

  ftcentroids_.clear();
  for ( GSIZET j=0; j<qmesh_.size(); j++ ) { // for each cube
    for ( GSIZET k=0; k<qmesh_[j].v.size(); k++ ) a += *qmesh_[j].v[k];
    a *= (1.0/static_cast<GFTYPE>(qmesh_[j].v.size()));
    ftcentroids_.push_back(a);
  }

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
  GString serr = "GGridBox::init3d: ";

  assert(gbasis_.size() == 3 && "Invalid dimension");

  GINT nelems = ne_[0]*ne_[1]*ne_[2];
  hmesh_.resize(nelems);


  // Get uniform element sizes:
  GTPoint<GFTYPE> dx(ndim_);
  for ( GSIZET j=0; j<ndim_; j++ ) 
    dx[j] = Lbox_[j] / static_cast<GFTYPE>(ne_[j]);

  // Compute vertices of all hexes ('cubes'), starting 
  GINT   n=0;
  GTPoint<GFTYPE> v0(3);     // starting sph. point of elem
  GTPoint<GFTYPE> dv(3);     // delta-point
  for ( GSIZET k=0; k<ne_[2]; k++ ) { 
    for ( GSIZET j=0; j<ne_[1]; j++ ) { 
      for ( GSIZET i=0; i<ne_[0]; i++ ) { 
        v0.x1 = i*dx.x1; v0.x2 = j*dx.x2; v0.x3 = k*dx.x3; 
        hmesh_[n].v1 = v0;
        dv.x1 = dx.x1 ; dv.x2 = 0.0  ; dv.x3 = 0.0   ;  hmesh_[n].v2 = v0 + dv;
        dv.x1 = dx.x1 ; dv.x2 = dx.x2; dv.x3 = 0.0   ;  hmesh_[n].v3 = v0 + dv;
        dv.x1 = 0.0   ; dv.x2 = dx.x2; dv.x3 = 0.0   ;  hmesh_[n].v4 = v0 + dv;
        dv.x1 = v0.x1 ; dv.x2 = 0.0  ; dv.x3 = dx.x3 ;  hmesh_[n].v5 = v0 + dv;
        dv.x1 = v0.x1 ; dv.x2 = dx.x2; dv.x3 = dx.x3 ;  hmesh_[n].v6 = v0 + dv;
        dv.x1 = v0.x1 ; dv.x2 = dx.x2; dv.x3 = dx.x3 ;  hmesh_[n].v7 = v0 + dv;
        dv.x1 = v0.x1 ; dv.x2 = 0.0  ; dv.x3 = dx.x3 ;  hmesh_[n].v8 = v0 + dv;
        n++;
      }
    }
  }
  
  // Compute centroids of all hexes ('cubes'):
  GTPoint<GFTYPE> a(3);

  ftcentroids_.clear();
  for ( GSIZET j=0; j<hmesh_.size(); j++ ) { // for each cube
    for ( GSIZET k=0; k<hmesh_[j].v.size(); k++ ) a += *hmesh_[j].v[k];
    a *= (1.0/static_cast<GFTYPE>(hmesh_[j].v.size()));
    ftcentroids_.push_back(a);
  }

} // end, method init3d



//**********************************************************************************
//**********************************************************************************
// METHOD : do_grid
// DESC   : Public entry point for grid computation
// ARGS   : grid: GGrid object
//          rank: MPI rank or partition id
// RETURNS: none.
//**********************************************************************************
void GGridBox::do_grid(GGrid &grid, GINT irank)
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
// DESC   : Build 2d element list. It's assumed that init3d has been called
//          prior to entry, and that the qmesh_ has beed set. This arrays
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
void GGridBox::do_grid2d(GGrid &grid, GINT irank)
{
  assert(gbasis_.size()>0 && "Must set basis first");
  assert(ndim_ == 2 && "Dimension must be 2");

  GString                      serr = "GGridBox::do_grid2d: ";
  GTPoint<GFTYPE>              cent;
  GTVector<GINT>               iind;
  GTVector<GINT>               I(3);
  GTVector<GINT>              *bdy_ind;
  GTVector<GBdyType>          *bdy_typ;
  GTVector<GINT>              *face_ind;
  GTVector<GFTYPE>             Ni;
  GElemList                   *gelems = &grid.elems();
  GElem_base                  *pelem;
  GTVector<GTVector<GFTYPE>>  *xNodes;
  GTVector<GTVector<GFTYPE>*> *xiNodes;
  GTVector<GTVector<GFTYPE>>   xgtmp(3);


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
        (*xNodes)[l] += Ni * (*(qmesh_[i].v[m]))[l];
      }
    }

    pelem->init(*xNodes);

    // With edge/face centroids set, compute bdy_nodes:
    for ( GSIZET j=0; j<2*ndim_; j++ ) { // cycle over all edges
      cent = pelem->edgeCentroid(j);
      if ( cent.x2 == P0_.x2 ) face_ind = &pelem->edge_indices(0);
      if ( cent.x1 == P1_.x1 ) face_ind = &pelem->edge_indices(1);
      if ( cent.x2 == P1_.x2 ) face_ind = &pelem->edge_indices(2);
      if ( cent.x1 == P0_.x1 ) face_ind = &pelem->edge_indices(3);
      for ( GSIZET k=0; k<face_ind->size(); k++ ) bdy_ind->push_back((*face_ind)[k]); 
    }


    gelems->push_back(pelem);

    // Set global bdy types at each bdy_ind (this is a coarse 
    // application; finer control may be exercised in callback):
    set_global_bdytypes_2d(*pelem);

    // Find global global interior and face start & stop indices represented 
    // locally within element:
    nfnodes = 0;
    for ( GSIZET j=0; j<(*gelems)[i]->nfaces(); j++ )  // get # face nodes
      nfnodes += (*gelems)[i]->face_indices(j).size();
    pelem->igbeg() = icurr;      // beginning global index
    pelem->igend() = icurr + pelem->nnodes()-1; // end global index
    pelem->ifbeg() = fcurr;
    pelem->ifend() = fcurr+nfnodes-1;// end global face index
    icurr += pelem->nnodes();
    fcurr += nfnodes;
  } // end of quad mesh loop

  // Can set individual nodes and internal bdy conditions
  // with callback here:
  if ( bdycallback_ != NULLPTR ) {
    (*bdycallback_)(*gelems);
  }

} // end of method do_grid2d



//**********************************************************************************
//**********************************************************************************
// METHOD : do_grid3d
// DESC   : Build 3d element list. It's assumed that init3d has been called
//          prior to entry, and that the bmesh_ has beed set. This arrays
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
void GGridBox::do_grid3d(GGrid &grid, GINT irank)
{

  assert(gbasis_.size()>0 && "Must set basis first");
  assert(ndim_ == 3 && "Dimension must be 3");

  GString                      serr = "GGridBox::do_grid3d: ";
  GTPoint<GFTYPE>              cent;
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



  if ( gdd_       == NULLPTR ) gdd_ = new GDD_base(nprocs_);

  // Get indirection indices to create elements
  // only for this task:
  gdd_->doDD(ftcentroids_, irank, iind);

  GSIZET i, n;
  GSIZET nfnodes;   // no. face indices
  GSIZET icurr = 0; // current global index
  GSIZET fcurr = 0; // current global face index
  for ( GSIZET n=0; n<iind.size(); n++ ) { // for each hex in irank's mesh
    i = iind[n];

    pelem = new GElem_base(GE_REGULAR, gbasis_);
    xNodes  = &pelem->xNodes();  // node spatial data
    xiNodes = &pelem->xiNodes(); // node ref interval data
    Ni.resize(pelem->nnodes()); // tensor product shape function
    bdy_ind = &pelem->bdy_indices(); // get bdy indices data member
    bdy_ind->clear();
    pelem->igbeg() = icurr;      // beginning global index
    pelem->igend() = icurr + pelem->nnodes()-1; // end global index
    for ( GSIZET l=0; l<ndim_; l++ ) { // loop over element Cart coords
      (*xNodes)[l] = 0.0;
      for ( GSIZET m=0; m<pow(2,ndim_); m++ ) { // loop over verts given in Cart coords
        I[0] = m;
        lshapefcn_->Ni(I, *xiNodes, Ni);
        (*xNodes)[l] += Ni * (*(hmesh_[i].v[m]))[l];
      }
    }

    pelem->init(*xNodes);

    // With edge/face centroids set, compute bdy_nodes:
    for ( GSIZET j=0; j<2*ndim_; j++ ) { // cycle over faces
      cent = pelem->faceCentroid(j);
      if ( cent.x2 == P0_.x2 ) face_ind = &pelem->edge_indices(0);
      if ( cent.x1 == P1_.x1 ) face_ind = &pelem->edge_indices(1);
      if ( cent.x2 == P1_.x2 ) face_ind = &pelem->edge_indices(2);
      if ( cent.x1 == P0_.x1 ) face_ind = &pelem->edge_indices(3);
      if ( cent.x3 == P0_.x3 ) face_ind = &pelem->edge_indices(4);
      if ( cent.x3 == P1_.x3 ) face_ind = &pelem->edge_indices(5);
      face_ind = &pelem->face_indices(0);
      for ( GSIZET k=0; k<face_ind->size(); k++ ) bdy_ind->push_back((*face_ind)[k]); 
    }

    gelems->push_back(pelem);
    // Set global bdy types at each bdy_ind (this is a coarse 
    // application; finer control may be exercised in callback):
    set_global_bdytypes_3d(*pelem);

    nfnodes = 0;
    for ( GSIZET j=0; j<(*gelems)[i]->nfaces(); j++ )  // get # face nodes
      nfnodes += (*gelems)[i]->face_indices(j).size();
    pelem->igbeg() = icurr;      // beginning global index
    pelem->igend() = icurr + pelem->nnodes()-1; // end global index
    pelem->ifbeg() = fcurr;
    pelem->ifend() = fcurr+nfnodes-1; // end global face index
    icurr += pelem->nnodes();
    fcurr += nfnodes;
  } // end of hex mesh loop

  // Can set individual nodes and internal bdy conditions
  // with callback here:
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
void GGridBox::set_bdy_callback(std::function<void(GElemList &)> &callback)
{
  bdycallback_  = &callback;
} // end of method set_bdy_callback



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
// METHOD : set_global_bdytypes_2d
// DESC   : Set global boundary condition types in 2d
//            NOTE: element node coordinate may change in this call! Also,
//                  bdy indices must already have been set prior to 
//                  entry.
// ARGS   : pelem : Element under consideration
// RETURNS: none.
//**********************************************************************************
void GGridBox::set_global_bdytypes_2d(GElem_base &pelem)
{

  GTVector<GINT>              *bdy_ind=&pelem.bdy_indices();
  GTVector<GBdyType>          *bdy_typ=&pelem.bdy_types  ();
  GTVector<GTVector<GFTYPE>>  *xNodes =&pelem.xNodes();
  GTVector<GTVector<GINT>>    *face_ind;

  GSIZET ib;
  if ( bPeriodic_[0] == GBDY_PERIODIC ) { // check for periodicity in x
    for ( GSIZET m=0; m<2; m++ ) { // for each x edge
      for ( GSIZET k=0; k<bdy_ind->size(); k++ ) { // for each face index
        ib = (*bdy_ind)[k];
        if ( (*xNodes)[0][ib] == P0_.x1 || (*xNodes)[0][ib] == P1_.x1 ) 
          bdy_typ->push_back( GBDY_PERIODIC );
          // Set right x-coord equal to that on left-most bdy:
          if ( (*xNodes)[0][ib] == P1_.x1 ) (*xNodes)[0][ib] = P0_.x1;
      }
    }
  }
  else {                                         // not x-periodic
    for ( GSIZET m=0; m<2; m++ ) { // for each x edge
      for ( GSIZET k=0; k<bdy_ind->size(); k++ ) { // for each bdy index
        ib = (*bdy_ind)[k];
        if ( (*xNodes)[0][ib] == P0_.x1 || (*xNodes)[0][ib] == P1_.x1 ) 
          bdy_typ->push_back( global_bdy_types_[2*m+1] );
      }
    }
  }

  if ( bPeriodic_[1] == GBDY_PERIODIC ) { // check for periodicity in y
    for ( GSIZET m=0; m<2; m++ ) { // for each y edge
      for ( GSIZET k=0; k<bdy_ind->size(); k++ ) { // for each bdy index
        ib = (*bdy_ind)[k];
        if ( (*xNodes)[1][ib] == P0_.x2 || (*xNodes)[1][ib] == P1_.x2 ) 
          bdy_typ->push_back( GBDY_PERIODIC );
          // Set top y-coord equal to that on bottom-most bdy:
          if ( (*xNodes)[1][ib] == P1_.x2 ) (*xNodes)[1][ib] = P0_.x2;
      }
    }
  }
  else {                                         // not y-periodic
    for ( GSIZET m=0; m<2; m++ ) { // for each y edge
      for ( GSIZET k=0; k<bdy_ind->size(); k++ ) { // for each bdy index
        ib = (*bdy_ind)[k];
        if ( (*xNodes)[0][ib] == P0_.x2 || (*xNodes)[0][ib] == P1_.x2 ) 
          bdy_typ->push_back( global_bdy_types_[2*m] );
      }
    }
  }
} // end, method set_global_bdytypes_2d


//**********************************************************************************
//**********************************************************************************
// METHOD : set_global_bdytypes_3d
// DESC   : Set global boundary conditions in 3d
//            NOTE: element node coordinate may change in this call! Also,
//                  bdy indices must already have been set prior to 
//                  entry.
// ARGS   : pelem : Element under consideration
// RETURNS: none.
//**********************************************************************************
void GGridBox::set_global_bdytypes_3d(GElem_base &pelem)
{
  GTVector<GINT>              *bdy_ind=&pelem.bdy_indices();
  GTVector<GINT>              *face_ind=&pelem.bdy_indices();
  GTVector<GBdyType>          *bdy_typ=&pelem.bdy_types  ();
  GTVector<GTVector<GFTYPE>>  *xNodes =&pelem.xNodes();

  GSIZET ib;
  if ( bPeriodic_[0] == GBDY_PERIODIC ) { // check for periodicity in x
    for ( GSIZET m=0; m<2; m++ ) { // for each x face
      face_ind = &pelem.face_indices(2*m+1);
      for ( GSIZET k=0; k<face_ind->size(); k++ ) { // for each bdy index
        ib = (*bdy_ind)[k];
        if ( (*xNodes)[0][ib] == P0_.x1 || (*xNodes)[0][ib] == P1_.x1 ) 
          bdy_typ->push_back( GBDY_PERIODIC );
          // Set right x-coord equal to that on left-most bdy:
          if ( (*xNodes)[0][ib] == P1_.x1 ) (*xNodes)[0][ib] = P0_.x1;
      }
    }
  }
  else {                                         // not x-periodic
    for ( GSIZET m=0; m<2; m++ ) { // for each x face
      face_ind = &pelem.face_indices(2*m+1);
      for ( GSIZET k=0; k<face_ind->size(); k++ ) { // for each bdy index
        ib = (*bdy_ind)[k];
        if ( (*xNodes)[0][ib] == P0_.x1 || (*xNodes)[0][ib] == P1_.x1 ) 
          bdy_typ->push_back( global_bdy_types_[2*m+1] );
      }
    }
  }

  if ( bPeriodic_[1] == GBDY_PERIODIC ) { // check for periodicity in y
    for ( GSIZET m=0; m<2; m++ ) { // for each y face
      face_ind = &pelem.face_indices(2*m);
      for ( GSIZET k=0; k<face_ind->size(); k++ ) { // for each bdy index
        ib = (*bdy_ind)[k];
        if ( (*xNodes)[1][ib] == P0_.x2 || (*xNodes)[1][ib] == P1_.x2 ) 
          bdy_typ->push_back( GBDY_PERIODIC );
          // Set top y-coord equal to that on bottom-most bdy:
          if ( (*xNodes)[1][ib] == P1_.x2 ) (*xNodes)[1][ib] = P0_.x2;
      }
    }
  }
  else {                                         // not y-periodic
    for ( GSIZET m=0; m<2; m++ ) { // for each y face
      face_ind = &pelem.face_indices(2*m);
      for ( GSIZET k=0; k<face_ind->size(); k++ ) { // for each bdy index
        ib = (*bdy_ind)[k];
        if ( (*xNodes)[0][ib] == P0_.x2 || (*xNodes)[0][ib] == P1_.x2 ) 
          bdy_typ->push_back( global_bdy_types_[2*m] );
      }
    }
  }

  if ( bPeriodic_[1] == GBDY_PERIODIC ) { // check for periodicity in z
    for ( GSIZET m=0; m<2; m++ ) { // for each z face
      face_ind = &pelem.face_indices(m+4);
      for ( GSIZET k=0; k<face_ind->size(); k++ ) { // for each bdy index
        ib = (*bdy_ind)[k];
        if ( (*xNodes)[1][ib] == P0_.x3 || (*xNodes)[1][ib] == P1_.x3 ) 
          bdy_typ->push_back( GBDY_PERIODIC );
          // Set top z-coord equal to that on bottom-most bdy:
          if ( (*xNodes)[2][ib] == P1_.x3 ) (*xNodes)[2][ib] = P0_.x3;
      }
    }
  }
  else {                                         // not z-periodic
    for ( GSIZET m=0; m<2; m++ ) { // for each z face
      face_ind = &pelem.face_indices(m+4);
      for ( GSIZET k=0; k<face_ind->size(); k++ ) { // for each bdy index
        ib = (*bdy_ind)[k];
        if ( (*xNodes)[0][ib] == P0_.x3 || (*xNodes)[0][ib] == P1_.x3 ) 
          bdy_typ->push_back( global_bdy_types_[m+4] );
      }
    }
  }
} // end, method set_global_bdytypes_3d
