//==================================================================================
// Module       : ggrid
// Date         : 8/31/18 (DLR)
// Description  : GeoFLOW grid object. Is primarily a container for
//                a list of elements, as an indirection buffer
//                that breaks down elements by type.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#include <memory>
#include <cmath>
#include <limits>
#include <typeinfo>
#include "gtypes.h"
#include "gelem_base.hpp"
#include "ggrid.hpp"
#include "gmass.hpp"
#include "gcomm.hpp"
#include "ggfx.hpp"
#include "gutils.hpp"
#include "tbox/error_handler.hpp"

using namespace std;

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Default constructor
// ARGS   : ptree: property tree
//          b    : GNBasis
//          comm : communicator
// RETURNS: none
//**********************************************************************************
GGrid::GGrid(const geoflow::tbox::PropertyTree &ptree, GTVector<GNBasis<GCTYPE,GFTYPE>*> &b, GC_COMM &comm)
:
bInitialized_                   (FALSE),
bapplybc_                       (FALSE),
bupdatebc_                      (FALSE),
do_face_normals_                (FALSE),
nprocs_        (GComm::WorldSize(comm)),
ngelems_                            (0),
irank_         (GComm::WorldRank(comm)),
minnodedist_   (std::numeric_limits<GFTYPE>::max()),
volume_                           (0.0),
ivolume_                          (0.0),
comm_                            (comm),
mass_                         (NULLPTR),
imass_                        (NULLPTR),
ggfx_                         (NULLPTR),
ptree_                          (ptree),
bdy_apply_callback_           (NULLPTR),
bdy_update_callback_          (NULLPTR)
{
} // end of constructor method (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
GGrid::~GGrid()
{
  if ( mass_ != NULLPTR ) delete mass_;
  if ( imass_ != NULLPTR ) delete imass_;
  for ( auto j=0; j<gelems_.size(); j++ ) {
    if ( gelems_[j] != NULLPTR ) delete gelems_[j];
  }
} // end, destructor


//**********************************************************************************
//**********************************************************************************
// METHOD : do_typing
// DESC   : Do classification of element list into its
//          types, and store in member data structure. Each
//          type element in structure is an index into gelems_
//          corresp to that type
//          array.
// ARGS   : none.
// RETURNS: none.
//**********************************************************************************
void GGrid::do_typing()
{
  GTVector<GElemType> itmp(gelems_.size());

  GSIZET *ind=NULLPTR;
  GSIZET  nd=0;
  GSIZET  nfound;
  for ( auto j=0; j<gelems_.size(); j++ ) itmp[j] = gelems_[j]->elemtype();

  itype_.resize(GE_MAX);
  ntype_.resize(GE_MAX);
  ntype_ = 0;
  for ( auto j=0; j<GE_MAX; j++ ) {
    nfound = itmp.contains(static_cast<GElemType>(j),ind,nd);
    for ( auto i=0; i<nfound; i++ ) itype_[j].push_back(ind[i]);
    ntype_[j] = nfound;
  }
  if ( ind != NULLPTR ) delete [] ind;

  // Do sanity check:
  assert(ntype_.sum() == gelems_.size() && "Element typing failed");
} // end of method do_typing


//**********************************************************************************
//**********************************************************************************
// METHOD : print
// DESC   : Write basic grid info to specified file for plotting
// ARGS   : filename : string containing file name 
//          bdof     : print internal dof too?
// RETURNS: none.
//**********************************************************************************
void GGrid::print(GString filename, GBOOL bdof)
{
  GString serr = "GridIcos::print: ";
  std::ofstream ios;

  GSIZET n;
  GTVector<GINT> N;
  GTPoint<GFTYPE> *vert;
  GTVector<GTVector<GFTYPE>> *xnodes;

  ios.open(filename);

  if ( !bdof ) { // print only wire mesh (vertices)
    for ( auto i=0; i<gelems_.size(); i++ ) { // for each element
       
       for ( auto j=0; j<gelems_[i]->nvertices(); j++ ) { // for each vertex of element
        vert = &gelems_[i]->xVertices(j);
        for ( auto k=0; k<vert->size()-1; k++ ) {
            ios << (*vert)[k] << " ";
        }
        ios << (*vert)[vert->size()-1] << std::endl;
         
      }
    }
    return;
  }

  // Print internal dofs too:
  for ( auto i=0; i<gelems_.size(); i++ ) { // for each element
    xnodes  = &gelems_[i]->xNodes();  
    N       = gelems_[i]->size();
    n       = (*xnodes)[0].size();
    if ( GDIM == 1 ) {
      for ( auto k=0; k<n; k++ ) {
        ios << (*xnodes)[0][k] << " " << std::endl;
      }
    }
    else if ( GDIM == 2 && gelems_[i]->elemtype() != GE_2DEMBEDDED ) {
      xnodes  = &gelems_[i]->xNodes();  
      for ( auto k=0; k<n; k++ ) {
        ios << (*xnodes)[0][k] << " "
            << (*xnodes)[1][k] << std::endl;
      }
    }
    else if (GDIM==2 && gelems_[i]->elemtype() == GE_2DEMBEDDED ) {
      // Lay down in separate 'sub-quads':
      for ( auto k=0; k<N[1]-1; k++ ) {
        for ( auto j=0; j<N[0]-1; j++ ) {
          // For each sub-quad, print its vertices:
          ios << (*xnodes)[0][j+k*N[0]] << " "
              << (*xnodes)[1][j+k*N[0]] << " "
              << (*xnodes)[2][j+k*N[0]] << std::endl;
          ios << (*xnodes)[0][j+1+k*N[0]] << " "
              << (*xnodes)[1][j+1+k*N[0]] << " "
              << (*xnodes)[2][j+1+k*N[0]] << std::endl;
          ios << (*xnodes)[0][j+1+(k+1)*N[0]] << " "
              << (*xnodes)[1][j+1+(k+1)*N[0]] << " "
              << (*xnodes)[2][j+1+(k+1)*N[0]] << std::endl;
          ios << (*xnodes)[0][j+(k+1)*N[0]] << " "
              << (*xnodes)[1][j+(k+1)*N[0]] << " "
              << (*xnodes)[2][j+(k+1)*N[0]] << std::endl;
        }
      }
    }
    else if ( GDIM == 3 ) {
      for ( auto k=0; k<n-1; k++ ) {
        ios << (*xnodes)[0][k] << " "
            << (*xnodes)[1][k] << " "
            << (*xnodes)[2][k] << std::endl;
      }
    }
  }


  ios.close();

} // end of method print


//**********************************************************************************
//**********************************************************************************
// METHOD :  << operator method (1)
// DESC   : output stream operator
// ARGS   :
// RETURNS: ostream &
//**********************************************************************************
std::ostream &operator<<(std::ostream &str, GGrid &e)
{
  return str;
} // end of operator <<


//**********************************************************************************
//**********************************************************************************
// METHOD : ndof
// DESC   : Find number of dof in grid
// ARGS   : none
// RETURNS: GSIZET number local dof
//**********************************************************************************
GSIZET GGrid::ndof()
{
   assert(gelems_.size() > 0 && "Elements not set");

   GSIZET Ntot=0;
   for ( auto i=0; i<gelems_.size(); i++ ) Ntot += gelems_[i]->nnodes();

   return Ntot;
} // end of method ndof


//**********************************************************************************
//**********************************************************************************
// METHOD : nfacedof
// DESC   : Find number of elem face (not bdy) dof in grid. 
// ARGS   : none
// RETURNS: GSIZET number surface/face  dof
//**********************************************************************************
GSIZET GGrid::nfacedof()
{
   return gieface_.size();
} // end of method nfacedof


//**********************************************************************************
//**********************************************************************************
// METHOD : nbdydof
// DESC   : Find number of bdy dof in grid. These will include
//          global domain boundary _and_ may also include embedded
//          boundary surface nodes.
// ARGS   : none
// RETURNS: GSIZET number surface dof
//**********************************************************************************
GSIZET GGrid::nbdydof()
{
   return igbdt_.size();;
} // end of method nbdydof


//**********************************************************************************
//**********************************************************************************
// METHOD : minlength
// DESC   : Find elem length separation
// ARGS   : none
// RETURNS: GFTYPE separation
//**********************************************************************************
GFTYPE GGrid::minlength()
{
   assert(gelems_.size() > 0 && "Elements not set");

   GFTYPE lmin, gmin;
   GTPoint<GFTYPE> dr;
   GTVector<GTPoint<GFTYPE>> *xverts;

   lmin = std::numeric_limits<GFTYPE>::max();
   for ( auto i=0; i<gelems_.size(); i++ ) {
     xverts = &gelems_[i]->xVertices();
     #if defined(_G_IS2D)
     for ( auto j=0; j<xverts->size(); j++ ) {
       dr = (*xverts)[(j+1)%xverts->size()] - (*xverts)[j];
       lmin = MIN(lmin,dr.norm());
     }
     #elif defined(_G_IS3D)
     for ( auto j=0; j<4; j++ ) { // bottom
       dr = (*xverts)[(j+1)%xverts->size()] - (*xverts)[j];
       lmin = MIN(lmin,dr.norm());
     }
     for ( auto j=4; j<8; j++ ) { // top
       dr = (*xverts)[(j+1)%xverts->size()] - (*xverts)[j];
       lmin = MIN(lmin,dr.norm());
     }
     for ( auto j=0; j<4; j++ ) { // vertical edges
       dr = (*xverts)[j+4] - (*xverts)[j];
       lmin = MIN(lmin,dr.norm());
     }
     #endif
   }

   GComm::Allreduce(&lmin, &gmin, 1, T2GCDatatype<GFTYPE>() , GC_OP_MIN, comm_);

   return gmin;
} // end of method minlength


//**********************************************************************************
//**********************************************************************************
// METHOD : maxlength
// DESC   : Find max elem length 
// ARGS   : none
// RETURNS: GFTYPE length
//**********************************************************************************
GFTYPE GGrid::maxlength()
{
   assert(gelems_.size() > 0 && "Elements not set");

   GFTYPE lmax, gmax;
   GTPoint<GFTYPE> dr;
   GTVector<GTPoint<GFTYPE>> *xverts;

   lmax = 0.0;
   for ( auto i=0; i<gelems_.size(); i++ ) {
     xverts = &gelems_[i]->xVertices();
     #if defined(_G_IS2D)
     for ( auto j=0; j<xverts->size(); j++ ) {
       dr = (*xverts)[(j+1)%xverts->size()] - (*xverts)[j];
       lmax = MAX(lmax,dr.norm());
     }
     #elif defined(_G_IS3D)
     for ( auto j=0; j<4; j++ ) { // bottom
       dr = (*xverts)[(j+1)%xverts->size()] - (*xverts)[j];
       lmax = MAX(lmax,dr.norm());
     }
     for ( auto j=4; j<8; j++ ) { // top
       dr = (*xverts)[(j+1)%xverts->size()] - (*xverts)[j];
       lmax = MAX(lmax,dr.norm());
     }
     for ( auto j=0; j<4; j++ ) { // vertical edges
       dr = (*xverts)[j+4] - (*xverts)[j];
       lmax = MAX(lmax,dr.norm());
     }
     #endif
   }

   GComm::Allreduce(&lmax, &gmax, 1, T2GCDatatype<GFTYPE>() , GC_OP_MAX, comm_);

   return gmax;
} // end of method maxlength


//**********************************************************************************
//**********************************************************************************
// METHOD : avglength
// DESC   : Find avg elem length 
// ARGS   : none
// RETURNS: GFTYPE length
//**********************************************************************************
GFTYPE GGrid::avglength()
{
   assert(gelems_.size() > 0 && "Elements not set");

   GFTYPE gavg, lavg, navg, lv[2], gv[2];;
   GTPoint<GFTYPE> dr;
   GTVector<GTPoint<GFTYPE>> *xverts;

   navg = 0.0;
   lavg = 0.0;
   for ( auto i=0; i<gelems_.size(); i++ ) {
     xverts = &gelems_[i]->xVertices();
     #if defined(_G_IS2D)
     for ( auto j=0; j<xverts->size(); j++ ) {
       dr = (*xverts)[(j+1)%xverts->size()] - (*xverts)[j];
       lavg += dr.norm();
       navg += 1.0;
     }
     #elif defined(_G_IS3D)
     for ( auto j=0; j<4; j++ ) { // bottom
       dr = (*xverts)[(j+1)%xverts->size()] - (*xverts)[j];
       lavg += dr.norm();
       navg += 1.0;
     }
     for ( auto j=4; j<8; j++ ) { // top
       dr = (*xverts)[(j+1)%xverts->size()] - (*xverts)[j];
       lavg += dr.norm();
       navg += 1.0;
     }
     for ( auto j=0; j<4; j++ ) { // vertical edges
       dr = (*xverts)[j+4] - (*xverts)[j];
       lavg += dr.norm();
       navg += 1.0;
     }
     #endif
   }

   lv[0] = lavg;
   lv[1] = navg;
   GComm::Allreduce(&lv, &gv, 2, T2GCDatatype<GFTYPE>() , GC_OP_SUM, comm_);
   gavg = gv[0] / gv[1];

   return gavg;
} // end of method avglength


//**********************************************************************************
//**********************************************************************************
// METHOD : grid_init (1)
// DESC   : Initialize global (metric) variables. All elements are assumed to be
//          of the same type.
// ARGS   : none
// RETURNS: none
//**********************************************************************************
void GGrid::grid_init()
{

  GTimerStart("GGrid::grid_init: do_elems");
  do_elems(); // generate element list from derived class
  GTimerStop("GGrid::grid_init: do_elems");

  GComm::Synch(comm_);

  GTimerStart("GGrid::grid_init: do_typing");
  do_typing(); // do element-typing check
  GTimerStop("GGrid::grid_init: do_typing");


  // Have elements been set yet?
  assert(gelems_.size() > 0 && "Elements not set");

  // Restrict grid to a single element type:
  assert(ntype_.multiplicity(0) == GE_MAX-1
        && "Only a single element type allowed on grid");

  if      ( itype_[GE_2DEMBEDDED].size() > 0 ) gtype_ = GE_2DEMBEDDED;
  else if ( itype_[GE_DEFORMED]  .size() > 0 ) gtype_ = GE_DEFORMED;
  else if ( itype_[GE_REGULAR]   .size() > 0 ) gtype_ = GE_REGULAR;

  globalize_coords    (); // set glob vec of node coords
  init_local_face_info(); // find glob vec of face indices


  GTimerStart("GGrid::grid_init: init_bc_info");
  // All element bdy/face data should have been set by now:


  init_bc_info();
  GTimerStop("GGrid::grid_init: init_bc_info");


  GTimerStart("GGrid::grid_init: def_geom_init");
  if ( gtype_ == GE_2DEMBEDDED || gtype_  == GE_DEFORMED ) {
    def_geom_init();
  }
  GTimerStop("GGrid::grid_init: def_geom_init");

  GTimerStart("GGrid::grid_init: reg_geom_init");
  if ( gtype_ == GE_REGULAR ) {
    reg_geom_init();
  }
  GTimerStop("GGrid::grid_init: reg_geom_init");

  // Get global number of elements:
  GSIZET nelems = gelems_.size();
  GComm::Allreduce(&nelems, &ngelems_, 1, T2GCDatatype<GSIZET>() , GC_OP_SUM, comm_);

  bInitialized_ = TRUE;

  mass_ = new GMass(*this);
  
  GTimerStart("GGrid::grid_init: find_min_dist");
  minnodedist_ = find_min_dist();
  GTimerStop("GGrid::grid_init: find_min_dist");

  // Compute (global) grid volume:
  GTVector<GFTYPE> tmp0(ndof()), tmp1(ndof());
  tmp0 = 1.0;
  volume_  = integrate(tmp0, tmp1);
  ivolume_ = 1.0 / volume_;
} // end of method grid_init (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : grid_init (2)
// DESC   : Initialize global (metric) variables. All elements are assumed to be
//          of the same type. Called for restart.
// ARGS   : none
// RETURNS: none
//**********************************************************************************
void GGrid::grid_init(GTMatrix<GINT> &p,
                      GTVector<GTVector<GFTYPE>> &xnodes)
{

  GTimerStart("GGrid::grid_init: do_elems");
  do_elems(p, xnodes); // generate element list from derived class
  GTimerStop("GGrid::grid_init: do_elems");

  GComm::Synch(comm_);

  GTimerStart("GGrid::grid_init: do_typing");
  do_typing(); // do element-typing check
  GTimerStop("GGrid::grid_init: do_typing");

  // Have elements been set yet?
  assert(gelems_.size() > 0 && "Elements not set");


  // Restrict grid to a single element type:
  assert(ntype_.multiplicity(0) == GE_MAX-1
        && "Only a single element type allowed on grid");


  if      ( itype_[GE_2DEMBEDDED].size() > 0 ) gtype_ = GE_2DEMBEDDED;
  else if ( itype_[GE_DEFORMED]  .size() > 0 ) gtype_ = GE_DEFORMED;
  else if ( itype_[GE_REGULAR]   .size() > 0 ) gtype_ = GE_REGULAR;

  globalize_coords    (); // set glob vec of node coords
  init_local_face_info(); // find glob vec of face indices

  GTimerStart("GGrid::grid_init: init_bc_info");
  // All element bdy/face data should have been set by now:
  init_bc_info();
  GTimerStop("GGrid::grid_init: init_bc_info");

  GTimerStart("GGrid::grid_init: def_geom_init");
  if ( itype_[GE_2DEMBEDDED].size() > 0
    || itype_  [GE_DEFORMED].size() > 0 ) {
    def_geom_init();
  }
  GTimerStop("GGrid::grid_init: def_geom_init");

  GTimerStart("GGrid::grid_init: reg_geom_init");
  if ( itype_[GE_REGULAR].size() > 0 ) {
    reg_geom_init();
  }
  GTimerStop("GGrid::grid_init: reg_geom_init");


  GTimerStart("GGrid::grid_init: find_min_dist");
  minnodedist_ = find_min_dist();
  GTimerStop("GGrid::grid_init: find_min_dist");

  bInitialized_ = TRUE;

  // Compute grid volume:
  GTVector<GFTYPE> tmp0(ndof()), tmp1(ndof());
  tmp0 = 1.0;
  volume_  = integrate(tmp0, tmp1);
  ivolume_ = 1.0 / volume_;

} // end of method grid_init


//**********************************************************************************
//**********************************************************************************
// METHOD : def_geom_init
// DESC   : Initialize global (metric) variables for deformed elems. 
//          All elements are assumed to be of the same type. Note:
//          do_elems and globalize_coords must be called prior to entry.
// ARGS   : none
// RETURNS: none
//**********************************************************************************
void GGrid::def_geom_init()
{
   assert(gelems_.size() > 0 && "Elements not set");
   assert(gtype_ == GE_2DEMBEDDED
       || gtype_ == GE_DEFORMED && "Invalid element type");


   GString serr = "GGrid::def_geom_init: ";
   GSIZET nxy = gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;
   GTVector<GTVector<GFTYPE>> *xe;

   // Resize geometric quantities to global size:
   dXidX_.resize(nxy,nxy);
   dXdXi_.resize(nxy,nxy);
   for ( auto j=0; j<nxy; j++ ) {
     for ( auto i=0; i<nxy; i++ )  {
       dXidX_(i,j).resize(ndof());
       dXdXi_(i,j).resize(ndof());
     }
   }
   Jac_.resize(ndof());
   faceJac_.resize(nfacedof());

   // Resize surface-point-wise normals:
   faceNormals_.resize(nxy); // no. coords for each normal at each face point
   bdyNormals_.resize(nxy); // no. coords for each normal at each domain bdy point
   for ( auto i=0; i<bdyNormals_.size(); i++ ) {
     faceNormals_[i].resize(nfacedof());
     bdyNormals_ [i].resize(nbdydof());
   }

   // Now, set the geometry/metric quanties from the elements:
   GSIZET ibeg, iend; // beg, end indices for global arrays
   GSIZET ibbeg, ibend; // beg, end indices for global arrays for bdy quantities
   GSIZET ifbeg, ifend; // beg, end indices for global arrays for face quantities
   for ( auto e=0; e<gelems_.size(); e++ ) {
     ibeg  = gelems_[e]->igbeg(); iend  = gelems_[e]->igend();
     ifbeg = gelems_[e]->ifbeg(); ifend = gelems_[e]->ifend();
//   ibbeg = gelems_[e]->ibbeg(); ibend = gelems_[e]->ibend();

     xe    = &gelems_[e]->xNodes();

     // Restrict global arrays to local scope:
     for ( auto j=0; j<nxy; j++ ) {
       for ( auto i=0; i<nxy; i++ )  {
         dXidX_(i,j).range(ibeg, iend);
         dXdXi_(i,j).range(ibeg, iend);
       }
     }
     Jac_.range(ibeg, iend);
     faceJac_.range(ifbeg, ifend);

     // Set the geom/metric quantities using element data:
     if ( GDIM == 2 ) {
       gelems_[e]->dogeom2d(dXdXi_, dXidX_, Jac_, faceJac_);
     }
     else if ( GDIM == 3 ) {
       gelems_[e]->dogeom3d(dXdXi_, dXidX_, Jac_, faceJac_);
     }

     // Zero-out local xe; only global allowed now:
//   for ( auto j=0; j<nxy; j++ ) (*xe)[j].clear(); 
     
   } // end, element loop

   for ( auto j=0; j<nxy; j++ )  {
     for ( auto i=0; i<nxy; i++ )  {
       dXidX_(i,j).range_reset();
       dXdXi_(i,j).range_reset();
     }
   }
   Jac_.range_reset();
   faceJac_.range_reset();

   do_normals();

   GComm::Synch(comm_);
   
} // end of method def_geom_init


//**********************************************************************************
//**********************************************************************************
// METHOD : reg_geom_init
// DESC   : Initialize global (metric) variables for regular elemes. 
//          All elements are assumed to be of the same type. Note:
//          do_elems and globalize_coords must be called prior to entry.
// ARGS   : none
// RETURNS: none
//**********************************************************************************
void GGrid::reg_geom_init()
{
   assert(gelems_.size() > 0 && "Elements not set");
   assert( gtype_ == GE_REGULAR && "Invalid element type");

   GString serr = "GridIcos::reg_geom_init: ";
   GSIZET nxy = GDIM;
   GTMatrix<GTVector<GFTYPE>>  dXdXi_;
   GTVector<GTVector<GFTYPE>> *xe;

   // Resize geometric quantities to global size:
   dXidX_.resize(nxy,1);
   dXdXi_.resize(nxy,1);
   for ( auto i=0; i<nxy; i++ ) {
     dXidX_(i,0).resize(ndof());
     dXdXi_(i,0).resize(ndof());
   }
   Jac_.resize(ndof());
   faceJac_.resize(nfacedof());


   // Resize surface-point-wise normals:
   faceNormals_.resize(nxy); // no. coords for each normal at each face point
   bdyNormals_ .resize(nxy); // no. coords for each normal at each bdy point
   for ( auto i=0; i<nxy; i++ ) {
     faceNormals_[i].resize(nfacedof());
     bdyNormals_ [i].resize(nbdydof());
   }


   // Now, set the geometry/metric quanties from the elements:
   GSIZET ibeg, iend; // beg, end indices for global arrays
   GSIZET ibbeg, ibend; // beg, end indices for global arrays for bdy quantities
   GSIZET ifbeg, ifend; // beg, end indices for global arrays for face quantities
   for ( auto e=0; e<gelems_.size(); e++ ) {
     ibeg  = gelems_[e]->igbeg(); iend  = gelems_[e]->igend();
     ifbeg = gelems_[e]->ifbeg(); ifend = gelems_[e]->ifend();
//   ibbeg = gelems_[e]->ibbeg(); ibend = gelems_[e]->ibend();
  
     xe    = &gelems_[e]->xNodes();

     // Restrict global data to local scope:
     for ( auto j=0; j<nxy; j++ ) {
       faceNormals_[j].range(ifbeg, ifend); 
//     bdyNormals_ [j].range(ibbeg, ibend); 
     }
     for ( auto j=0; j<dXidX_.size(2); j++ ) {
       for ( auto i=0; i<dXidX_.size(1); i++ )  {
         dXidX_(i,j).range(ibeg, iend);
       }
     }
     Jac_.range(ibeg, iend);
     faceJac_.range(ifbeg, ifend);

     // Set the geom/metric quantities using element data:
     if ( GDIM == 2 ) {
       gelems_[e]->dogeom2d(dXdXi_, dXidX_, Jac_, faceJac_); 
     } 
     else if ( GDIM == 3 ) {
       gelems_[e]->dogeom3d(dXdXi_, dXidX_, Jac_, faceJac_);
     }
      
     // Zero-out local xe; only global allowed now:
//   for ( auto j=0; j<nxy; j++ ) (*xe)[j].clear(); 

   } // end, element loop

   // Reset global scope:
   for ( auto j=0; j<nxy; j++ ) {
     faceNormals_[j].range_reset();
     bdyNormals_ [j].range_reset();
   }
   for ( auto j=0; j<dXidX_.size(2); j++ )  {
     for ( auto i=0; i<dXidX_.size(1); i++ )  {
       dXidX_(i,j).range_reset();
     }
   }
   Jac_.range_reset();
   faceJac_.range_reset();

   GComm::Synch(comm_);
   
} // end of method reg_geom_init



//**********************************************************************************
//**********************************************************************************
// METHOD : do_normals
// DESC   : Compute normals to elem faces, and to domain boundary 
//          nodes. Note: init_bc_info must be called prior to entry.
// ARGS   : none
// RETURNS: none
//**********************************************************************************
void GGrid::do_normals()
{
  assert(gelems_.size() > 0 && "Elements not set");

  GString         serr = "GridIcos::do_normals: ";
  GSIZET          nxy = gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;
  GSIZET          n;
  GTPoint<GFTYPE> pt;

  // Note: at a given node, the (Cartesian) normals are
  // computed as n_j = (dX_/dxi  X  dX_/deta)_j, where
  // (xi, eta) define the surface, and X_ = (x, y, z)
  // are the physical coordinates. Knoweldge of the 

  // Set element face normals. Note: arrays for 
  // normals are allocated in these calls:
  if ( do_face_normals_ ) {
    do_face_normals();
  }
  
  // Set domain boundary node normals:
  do_bdy_normals(dXdXi_, igbdy_bdyface__, bdyNormals_, idepComp_);
   
} // end of method do_normals


//**********************************************************************************
//**********************************************************************************
// METHOD : globalize_coords
// DESC   : Create 'global' vectors from element-based coordinates
// ARGS   : none
// RETURNS: none
//**********************************************************************************
void GGrid::globalize_coords()
{
   GString serr = "GridIcos::globalize_coords: ";
   GSIZET  nxy = gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;
   GTVector<GTVector<GFTYPE>> *xe;

   xNodes_.resize(nxy);
   for ( auto j=0; j<nxy; j++ ) xNodes_[j].resize(ndof());

   GSIZET ibeg, iend; // beg, end indices for global arrays
   for ( auto e=0; e<gelems_.size(); e++ ) {
     ibeg  = gelems_[e]->igbeg(); iend  = gelems_[e]->igend();
     xe    = &gelems_[e]->xNodes();

     // Set global nodal Cart coords from element coords:
     for ( auto j=0; j<nxy; j++ ) {
       xNodes_[j].range(ibeg, iend);
       xNodes_[j] = (*xe)[j];
     }

   } // end, element loop

   // Reset global scope:
   for ( auto j=0; j<nxy; j++ ) xNodes_[j].range_reset();


} // end, method globalize_coords


//**********************************************************************************
//**********************************************************************************
// METHOD : dXidX (1)
// DESC   : return global metric terms
// ARGS   : none
// RETURNS: GTMatrix<GTVector<GFTYPE>> &
//**********************************************************************************
GTMatrix<GTVector<GFTYPE>> &GGrid::dXidX()
{
   assert(bInitialized_ && "Object not inititaized");
   return dXidX_;

} // end of method dXidX


//**********************************************************************************
//**********************************************************************************
// METHOD : dXidX (2)
// DESC   : return global metric element
// ARGS   : i,j : matrix element indices
// RETURNS: GTVector<GFTYPE> &
//**********************************************************************************
GTVector<GFTYPE> &GGrid::dXidX(GSIZET i, GSIZET j)
{
   assert(bInitialized_ && "Object not inititaized");
   return dXidX_(i,j);

} // end of method dXidX


//**********************************************************************************
//**********************************************************************************
// METHOD : massop
// DESC   : return global diagonal mass operator
// ARGS   : none
// RETURNS: GTVector<GFTYPE> &
//**********************************************************************************
GMass &GGrid::massop()
{
   assert(bInitialized_ && "Object not inititaized");
   return *mass_;

} // end of method massop


//**********************************************************************************
//**********************************************************************************
// METHOD : imassop
// DESC   : return inverse of global diagonal mass operator
// ARGS   : none
// RETURNS: GTVector<GFTYPE> &
//**********************************************************************************
GMass &GGrid::imassop()
{
   assert(bInitialized_ && "Object not inititaized");
   if ( imass_ == NULLPTR ) imass_ = new GMass(*this, TRUE);
   return *imass_;

} // end of method imassop


//**********************************************************************************
//**********************************************************************************
// METHOD : Jac
// DESC   : return global coord transform Jacobian
// ARGS   : none
// RETURNS: GTVector<GFTYPE> &
//**********************************************************************************
GTVector<GFTYPE> &GGrid::Jac()
{
   assert(bInitialized_ && "Object not inititaized");
   return Jac_;

} // end of method Jac


//**********************************************************************************
//**********************************************************************************
// METHOD : faceJac
// DESC   : Return global coord transform Jacobian for faces
// ARGS   : none
// RETURNS: GTVector<GFTYPE> &
//**********************************************************************************
GTVector<GFTYPE> &GGrid::faceJac()
{
   assert(bInitialized_ && "Object not inititaized");
   return faceJac_;

} // end of method faceJac


//**********************************************************************************
//**********************************************************************************
// METHOD : find_min_dist
// DESC   : Compute min inter-node distance for grid
//          NOTE: this method isn't recommended for regular
//                use, as it's slow.
// ARGS   : none
// RETURNS: none
//**********************************************************************************
GFTYPE GGrid::find_min_dist()
{
  assert(gelems_.size() > 0 && "Elements not set");

 
  GFTYPE           tiny = 1000.0*std::numeric_limits<GFTYPE>::epsilon();
  GTPoint<GFTYPE>  dx(3), p0(3), p1(3);

  // Find min node distance on grid:
  dx[2] = 0.0;
  p0[2] = 0.0;
  p1[2] = 0.0;

  GFTYPE xmin = std::numeric_limits<GFTYPE>::max();
  GFTYPE xgmin;

  GComm::Synch(comm_);

  // Get grid dimension:
  GINT nxy = gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;

  // Cycle over all node points and find min distance:
  // NOTE: this isn't fast....
  for ( auto j=0; j<xNodes_[0].size()-1; j++ ) {
    for ( auto i=0; i<nxy; i++ ) {
      p0[i] = xNodes_[i][j];
      p1[i] = xNodes_[i][j+1];
    }
    dx = p1 - p0; 
    if ( dx.norm() > tiny ) xmin = MIN(xmin, dx.norm()); 
  }
  GComm::Allreduce(&xmin, &xgmin, 1, T2GCDatatype<GFTYPE>() , GC_OP_MIN, comm_);
  return xgmin;

} // end of method find_min_dist


//**********************************************************************************
//**********************************************************************************
// METHOD : integrate
// DESC   : Compute spatial integral of global input vector. Result 
//          is a sum over all MPI tasks. NOTE: This method extracts weights
//          from the elements, so don't use this method very often.
// ARGS   : u     : 'global' integral argument
//          tmp   :  tmp vector, same size as u
//          blocal:  do global reduction?
// RETURNS: GFTYPE integral
//**********************************************************************************
GFTYPE GGrid::integrate(GTVector<GFTYPE> &u, GTVector<GFTYPE> &tmp, GBOOL bglobal)
{
  assert(bInitialized_ && "Object not inititaized");

  GSIZET                       ibeg, iend; // beg, end indices for global array
  GSIZET                       n;
  GFTYPE                       xint, xgint;
  GTVector<GINT>               N(GDIM);    // coord node sizes
  GTVector<GTVector<GFTYPE>*>  W(GDIM);    // element weights

  // NOTE: We could instantiate a mass matrix, but this would
  //       require computing a new GDIM-operator each time
  //       integration is required or setting it.
#if defined(_G_IS2D)
  for ( auto e=0; e<gelems_.size(); e++ ) {
    ibeg  = gelems_[e]->igbeg(); iend  = gelems_[e]->igend();

    // Restrict global data to local scope:
    tmp.range(ibeg, iend);
    u.range(ibeg, iend);

    for ( auto k=0; k<GDIM; k++ ) {
      W[k] = gelems_[e]->gbasis(k)->getWeights();
      N[k] = gelems_[e]->size(k);
    }
    n = 0;
    for ( auto k=0; k<N[1]; k++ ) {
      for ( auto j=0; j<N[0]; j++ ) {
        tmp[n] = (*W[1])[k]*(*W[0])[j] * u[n];
        n++;
      }
    }
  } // end, element loop
#elif defined(_G_IS3D)
  for ( auto e=0; e<gelems_.size(); e++ ) {
    ibeg  = gelems_[e]->igbeg(); iend  = gelems_[e]->igend();

    // Restrict global data to local scope:
    tmp.range(ibeg, iend);
    u.range(ibeg, iend);

    for ( auto k=0; k<GDIM; k++ ) {
      W[k] = gelems_[e]->gbasis(k)->getWeights();
      N[k] = gelems_[e]->size(k);
    }
    n = 0;
    for ( auto k=0; k<N[2]; k++ ) {
      for ( auto j=0; j<N[1]; j++ ) {
        for ( auto i=0; i<N[0]; i++ ) {
          tmp[n] = (*W[2])[k]*(*W[1])[j]*(*W[0])[i] * u[n];
          n++;
        }
      }
    }
  } // end, element loop
#endif
  tmp.range_reset(); 
  u.range_reset(); 

  // Multiply by Jacobian:
  tmp.pointProd(Jac_);
  xint = tmp.sum();  

  xgint = xint;
  if ( bglobal ) {
    GComm::Allreduce(&xint, &xgint, 1, T2GCDatatype<GFTYPE>() , GC_OP_SUM, comm_);
  }

  return xgint;

} // end of method integrate


//**********************************************************************************
//**********************************************************************************
// METHOD : deriv
// DESC   : Compute spatial derivative of u in direction idir, and
//          return in du.
// ARGS   : u   : 'global' integral argument
//          idir: coord wrt which to take derivative (1, 2, or 3)
//          utmp: tmp vector of same size as u
//          du  : derivtive, returned
// RETURNS: none.
//**********************************************************************************
void GGrid::deriv(GTVector<GFTYPE> &u, GINT idir, GTVector<GFTYPE> &utmp, 
                  GTVector<GFTYPE> &du)
{
  assert(bInitialized_ && "Object not inititialized");


  GTMatrix<GTVector<GFTYPE>> *dXidX = &this->dXidX();


  // du/dx_idir = Sum_j=[1:N] dxi_j/dx_idir D_j u:
  if ( this->gtype() == GE_REGULAR ) {
    assert(idir > 0 && idir <= GDIM && "Invalid derivative");
    GMTK::compute_grefderiv(*this, u, etmp_, idir, FALSE, du); // D_idir u
    du.pointProd((*dXidX)(idir-1, 0));
  }
  else {  // compute dXi_j/dX_idir D^j u:
    assert(idir > 0 && idir <= GDIM+1 && "Invalid derivative");
    GMTK::compute_grefderiv(*this, u, etmp_, 1, FALSE, du); // D_xi u
    du.pointProd((*dXidX)(0,idir-1));
    for ( auto j=1; j<GDIM; j++ ) {
      GMTK::compute_grefderiv(*this, u, etmp_, j+1, FALSE, utmp); // D_xi^j u
      utmp.pointProd((*dXidX)(j,idir-1));
      du += utmp; 

    }
  }
    
} // end of method deriv


//**********************************************************************************
//**********************************************************************************
// METHOD : init_local_face_info
// DESC   : Set local face info from element face data,
//          to be called after elements have been set.
// ARGS   : none
// RETURNS: none.
//**********************************************************************************
void GGrid::init_local_face_info()
{
  GBOOL                        bret;
  GSIZET                       ibeg, iend; // beg, end indices for global array
  GTVector<GINT>              *iebdy;  // domain bdy indices
  GTVector<GTVector<GINT>>    *ieface; // domain face indices

  GSIZET  m, n, nn; 
  GSIZET        ig; // index into global array

  n = 0;
  for ( auto e=0; e<gelems_.size(); e++ ) { // get global # face nodes
#if 0
    ieface = &gelems_[e]->face_indices(); // set in child class
    for ( auto j=0; j<ieface->size(); j++ ) { // count elem face nodes
      for ( auto k=0; k<(*ieface)[j].size(); k++) n++; 
    }
#endif
    n += gelems_[e]->nfnodes();
  }
  gieface_.resize(n);

  // Find list over all elemes of element face nodes:
  nn = 0; // global reference index
  m  = 0; // current num of global face indices
  for ( auto e=0; e<gelems_.size(); e++ ) { // cycle over all elems
    ibeg   = gelems_[e]->igbeg(); iend  = gelems_[e]->igend();
    ieface = &gelems_[e]->face_indices(); // set in child class
    for ( auto j=0; j<ieface->size(); j++ ) { // cycle over all elem faces
      for ( auto k=0; k<(*ieface)[j].size(); k++ ) {
        ig = nn + (*ieface)[j][k];
        if ( !gieface_.containsn(ig, m) ) { // don't include repeated face ind
          gieface_[m] = ig;
          m++;
        }
      }
    }
    nn += gelems_[e]->nnodes();
  } // end, element loop


} // end, init_local_face_info


//**********************************************************************************
//**********************************************************************************
// METHOD : init_bc_info
// DESC   : Set global bdy condition data from the element bdy data,
//          to be called after elements have been set.
// ARGS   : none
// RETURNS: none.
//**********************************************************************************
void GGrid::init_bc_info()
{
  GBOOL                        bret;
  GSIZET                       ibeg, iend; // beg, end indices for global array
  GBdyType                     btype;
  GTVector<GTVector<GBdyType>> igbdyt_bdyface;


  // Find boundary indices & types from config file 
  // specification, for _each_ natural/canonical domain face:
  config_bdy(ptree_, igbdy_bdyface_, igbdyt_bdyface);

  // Flatten bdy index indirection array:
  GSIZET      nind=0, nw=0;
  for ( auto j=0; j<igbdy_bdyface_.size(); j++ ) { // over dom can. bdy faces
    nind += igbdy_bdyface_[j].size(); // by-domain-face
  }
  igbdy_  .resize(nind); // indices of bdy nodes in volume

  nind = 0;
  for ( auto j=0; j<igbdy_bdyface_.size(); j++ ) { // over can. bdy faces
    for ( auto i=0; i<igbdy_bdyface_[j].size(); i++ ) {
      igbdy_  [nind] = igbdy_bdyface_ [j][i];
      nind++;
    }
  }

  // Create bdy type bins for each domain bdy:
  //   [Dom bdy][bdytype][volume index]:
  GSIZET          *ind=NULLPTR;
  GSIZET           n, nbdy, nind=0;
  GTVector<GSIZET> iunique;
    
  // Compute 'binned' structures for global indices
  // defining bdys, and to which index in bdy normal vectors
  // a given bdy node corresponds:
  n = 0; // cycle over all bdy nodes
  igbdy_binned_.resize(igbdy_bdyface_.size();
//ilbdy_binned_.resize(igbdy_bdyface_.size();
  for ( auto k=0; k<igbdy_bdyface_.size(); k++ ) { // cycle over canonical bdy face
    nbdy = igbdyt_bdyface[k].size();
    igbdy_binned_ [k].resize(GBDY_MAX);
//  ilbdy_binned_ [k].resize(GBDY_MAX); // index into bdy arrays for each bdy type
    for ( auto j=0; j<GBDY_MAX; j++ ) { // cycle over each bc type
      itype = static_cast<GBdyType>(j);
      val  = igbdyt_bdyface[k][itype];
      nbdy = igbdyt_bdyface[k].multiplicity(val, ind, nind);
      igbdy_binned_[k][j].resize(nbdy);
//    ilbdy_binned_[k][j].resize(nbdy);
      for ( auto i=0; i<nbdy; i++ ) { // assign comp. volume index
        igbdy_binned_[k][j][i] = igbdy_bdyface_[k][j][ind[i]];
//      ilbdy_binned_[k][j]    = n;     
        n++;
      }
    } // end, bdy cond type loop
  } // end, can. bdy loop
  if ( ind != NULLPTR ) delete [] ind;

  // Compute mask matrix from bdy vector:
  mask_.resize(ndof()); 
  mask_ = 1.0;
  for ( auto k=0; k<igbdy_binned_.size(); k++ ) {
    for ( auto j=0; j<GBDY_MAX; j++ ) {
      btype = std::static_cast<BdyType>(mbdy);
      if ( btype != GBDY_NONE && btype != GBDY_PERIODIC ) {
        for ( auto i=0; i<igbdy_binned_[k][j].size(); i++ ) {
          mask_[igbdy_binned_[k][j][i]] = 0.0;
        }
      }
    }
  }


} // end of method init_bc_info


//**********************************************************************************
//**********************************************************************************
// METHOD : add_terrain
// DESC   : Add bdy terrain by deforming coordinates. 
//          Works with only non-regular grids, and must be called 
//          after 'base' grid is completely constructed.
// ARGS   : xb  : boundary coordinates
//          utmp: tmp array pool
// RETURNS: none.
//**********************************************************************************
void GGrid::add_terrain(const State &xb, State &utmp)
{
   assert(gtype_ == GE_2DEMBEDDED
       || gtype_ == GE_DEFORMED && "Invalid element type");

   assert(utmp.size() > 5);

  // Set up temp arrays from pool:
  StateComp *x0=utmp[0];
  StateComp *b =utmp[1];
  State      tmp(utmp.size()-2);

  for ( auto j=0; j<tmp.size(); j++ ) tmp[j] = utmp[j+2];

  assert(ggfx_ != NULLPTR && "GGFX not set");

  // Construct solver and weak Laplacian operator:
  GCG<CGTypes>   cg(cgtraits_, *this, *ggfx_, tmp);
  GHelmholtz     H(*this);

  H.use_metric(TRUE); // do Laplacian in real space


  // Solve Nabla^2 (Xnew + Xb - XNodes) = 0 
  // for new (homgogeneous) grid solution, Xnew, 
  // given terrain, Xb, and // 'base' grid, XNodes:
  GTimerStart("GGrid::add_terrain: Solve");
  for ( auto j=0; j<xNodes_.size(); j++ ) {
   *b  = 0.0;
   *x0 = 0.0; // first guess
    cg.solve(H, *b, *xb[j], *x0);
    cout << "GGrid::add_terrain: xb_new[" << j << "]=" << *x0 << endl;
    xNodes_[j] = *x0;             // Reset XNodes = x0
//GPP(comm_,"GGrid::add_terrain: new_xNodes[" << j << "]=" << xNodes_[j]);
  }
  GTimerStop("GGrid::add_terrain: Solve");
 
  // Before computing new metric, Jacobian, etc, must set
  // new coordinates in elements that have already been initialized:
   GSIZET ibeg, iend; // beg, end indices for global arrays
   for ( auto e=0; e<gelems_.size(); e++ ) {
     ibeg  = gelems_[e]->igbeg(); iend  = gelems_[e]->igend();
     for ( auto j=0; j<xNodes_.size(); j++ ) xNodes_[j].range(ibeg, iend);
     gelems_[e]->set_nodes(xNodes_);
   }
   for ( auto j=0; j<xNodes_.size(); j++ ) xNodes_[j].range_reset();


  // Now, with new coordinates, recompute metric terms, 
  // Jacobian, normals:
  GTimerStart("GGrid::add_terrain: def_geom_init");
  def_geom_init();
  GTimerStop("GGrid::add_terrain: def_geom_init");

  // Compute new minimum node distance:
  GTimerStart("GGrid::grid_init: find_min_dist");
  minnodedist_ = find_min_dist();
  GTimerStop("GGrid::grid_init: find_min_dist");

  // Compute new (global) grid volume:
  *utmp[0] = 1.0;
  volume_  = integrate(*utmp[0], *utmp[1]);
  ivolume_ = 1.0 / volume_;

} // end of method add_terrain


