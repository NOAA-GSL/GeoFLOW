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
// ARGS   : P0, P1 : diag opposite corners of box
//          ne     : no. elements in each Cartesian direction
//          b      : vector of basis pointers. Number of elements in b is
//                   what _must_ be in L, and ne.
//          nprocs : no. MPI tasks
// RETURNS: none
//**********************************************************************************
GGridBox::GGridBox(GTPoint<GFTYPE> &P0, GTPoint<GFTYPE> &P1, GTVector<GINT> &ne, GTVector<GNBasis<GCTYPE,GFTYPE>*> &b, GINT nprocs) 
:
ndim_     (b.size()),
nprocs_     (nprocs),
gdd_       (NULLPTR),
lshapefcn_ (NULLPTR),
bdycallback_(NULLPTR)
{
  assert((b.size() == 2 || b.size() == 3) 
        && "Basis has incorrect dimensionalilty");
  assert(P0.size() >= ndim_ && P1.size() >= ndim_ && ne.size() >= b.size()
        && "Grid length and/or element count arrays of insufficient size");

  gbasis_.resize(b.size());
  gbasis_ = b;
  Lbox_.resize(b.size());
  ne_.resize(b.size());
  P0_ = P0;
  P1_ = P1;
  for ( GSIZET j=0; j<b.size(); j++ ) {
    Lbox_[j] = fabs(P1[j] - P0[j]);
    ne_  [j] = ne[j];
  }

  lshapefcn_ = new GShapeFcn_linear();
  if ( gbasis_.size() == 2 ) {
    init2d();
  }
  else if ( Lbox_.size() == 3 ) {
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

} // end, method do_grid


//**********************************************************************************
//**********************************************************************************
// METHOD : do_grid2d
// DESC   : Build 2d elemental grid. It's assumed that init3d has been called
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
  GSIZET icurr = 0; // current global index
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
    icurr += pelem->nnodes();
  } // end of hex mesh loop

  // If we have a callback function, set the boundary conditions here:
  if ( bdycallback_ != NULLPTR ) {
    (*bdycallback_)(grid);
  }

} // end of method do_grid2d


//**********************************************************************************
//**********************************************************************************
// METHOD : do_grid3d
// DESC   : Build d elemental grid. It's assumed that init3d has been called
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
  GSIZET icurr = 0; // current global index
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
    icurr += pelem->nnodes();
  } // end of hex mesh loop

  // If we have a callback function, set the boundary conditions here:
  if ( bdycallback_ != NULLPTR ) {
    (*bdycallback_)(grid);
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
void GGridBox::set_bdy_callback(std::function<void(GGrid &)> &callback)
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
// METHOD : reorderverts2d
// DESC   : Reorder specified vertices to be consistent with
//          shape functions
// ARGS   : uverts : list of unordered vertices
//          overts : array of ordered vertices, returned
// RETURNS: none.
//**********************************************************************************
void GGridBox::reorderverts2d(GTVector<GTPoint<GFTYPE>> &uverts, GTVector<GSIZET> &isort, 
                               GTVector<GTPoint<GFTYPE>> &overts)
{
  GString serr = "GGridBox::reorderverts2d: ";

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


