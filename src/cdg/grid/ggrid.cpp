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
#include "gelem_base.hpp"
#include "ggrid.hpp"
#include "gcomm.hpp"


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
nprocs_        (GComm::WorldSize(comm)),
irank_         (GComm::WorldRank(comm)),
bdycallback_                  (NULLPTR),
comm_                            (comm)
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
  for ( GSIZET j=0; j<gelems_.size(); j++ ) itmp[j] = gelems_[j]->elemtype();

  itype_.resize(GE_MAX);
  ntype_.resize(GE_MAX);
  ntype_ = 0;
  for ( GSIZET j=0; j<GE_MAX; j++ ) {
    itmp.contains(static_cast<GElemType>(j),ind,nd);
    for ( GSIZET i=0; i<nd; i++ ) itype_[j].push_back(ind[i]);
    ntype_[j] = nd;
  }
  if ( ind != NULLPTR ) delete [] ind;

} // end of method do_typing


#if 0
//**********************************************************************************
//**********************************************************************************
// METHOD : build
// DESC   : Do build and return of GGrid object
// ARGS   : ptree : property tree
//          gbasis: basis object
//          comm  : GC_Comm object
// RETURNS: GGrid object ptr
//**********************************************************************************
GGrid *GGrid::build(const geoflow::tbox::PropertyTree& ptree, GTVector<GNBasis<GCTYPE,GFTYPE>*> gbasis, GC_COMM comm)
{

  GString gname = ptree.getValue<GString>("grid_name");

  if      ( gname.compare  ("grid_icos") == 0   // 2d or 3d Icos grid
      ||    gname.compare("grid_sphere") == 0 ) {
    return new GGridIcos(ptree, gbasis, comm);
  }
  else if ( gname.compare("grid_box") == 0 ) { // 2d or 3D Cart grid
    return new GGridBox(ptree, gbasis, comm);
  }
  else {
    assert(FALSE && "Invalid PropertyTree grid specification");
  }

} // end, factory method build
#endif


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
    for ( GSIZET i=0; i<gelems_.size(); i++ ) { // for each element
       
       for ( GSIZET j=0; j<gelems_[i]->nvertices(); j++ ) { // for each vertex of element
        vert = &gelems_[i]->xVertices(j);
        for ( GSIZET k=0; k<vert->size()-1; k++ ) {
            ios << (*vert)[k] << " ";
        }
        ios << (*vert)[vert->size()-1] << std::endl;
         
      }
    }
    return;
  }

  // Print internal dofs too:
  for ( GSIZET i=0; i<gelems_.size(); i++ ) { // for each element
    xnodes  = &gelems_[i]->xNodes();  
    N       = gelems_[i]->size();
    n       = (*xnodes)[0].size();
    if ( GDIM == 1 ) {
      for ( GSIZET k=0; k<n; k++ ) {
        ios << (*xnodes)[0][k] << " " << std::endl;
      }
    }
    else if ( GDIM == 2 && gelems_[i]->elemtype() != GE_2DEMBEDDED ) {
      xnodes  = &gelems_[i]->xNodes();  
      for ( GSIZET k=0; k<n; k++ ) {
        ios << (*xnodes)[0][k] << " "
            << (*xnodes)[1][k] << std::endl;
      }
    }
    else if (GDIM==2 && gelems_[i]->elemtype() == GE_2DEMBEDDED ) {
      // Lay down in separate 'sub-quads':
      for ( GSIZET k=0; k<N[1]-1; k++ ) {
        for ( GSIZET j=0; j<N[0]-1; j++ ) {
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
      for ( GSIZET k=0; k<n-1; k++ ) {
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
   for ( GSIZET i=0; i<gelems_.size(); i++ ) Ntot += gelems_[i]->nnodes();

   return Ntot;
} // end of method ndof


//**********************************************************************************
//**********************************************************************************
// METHOD : nsurfdof
// DESC   : Find number of surface dof in grid
// ARGS   : none
// RETURNS: GSIZET number surface dof
//**********************************************************************************
GSIZET GGrid::nsurfdof()
{
   assert(gelems_.size() > 0 && "Elements not set");

   GSIZET Ntot=0;
   for ( GSIZET i=0; i<gelems_.size(); i++ ) {
      for ( GSIZET j=0; j<gelems_[i]->nfaces(); j++ ) 
        Ntot += gelems_[i]->face_indices(j).size();
   }

   return Ntot;
} // end of method nsurfdof


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
   for ( GSIZET i=0; i<gelems_.size(); i++ ) {
     xverts = &gelems_[i]->xVertices();
     #if defined(_G_IS2D)
     for ( GSIZET j=0; j<xverts->size(); j++ ) {
       dr = (*xverts)[(j+1)%xverts->size()] - (*xverts)[j];
       lmin = MIN(lmin,dr.norm());
     }
     #elif defined(_G_IS3D)
     for ( GSIZET j=0; j<4; j++ ) { // bottom
       dr = (*xverts)[(j+1)%xverts->size()] - (*xverts)[j];
       lmin = MIN(lmin,dr.norm());
     }
     for ( GSIZET j=4; j<8; j++ ) { // top
       dr = (*xverts)[(j+1)%xverts->size()] - (*xverts)[j];
       lmin = MIN(lmin,dr.norm());
     }
     for ( GSIZET j=0; j<4; j++ ) { // vertical edges
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
// RETURNS: GFTYPE separation
//**********************************************************************************
GFTYPE GGrid::maxlength()
{
   assert(gelems_.size() > 0 && "Elements not set");

   GFTYPE lmax, gmax;
   GTPoint<GFTYPE> dr;
   GTVector<GTPoint<GFTYPE>> *xverts;

   lmax = 0.0;
   for ( GSIZET i=0; i<gelems_.size(); i++ ) {
     xverts = &gelems_[i]->xVertices();
     #if defined(_G_IS2D)
     for ( GSIZET j=0; j<xverts->size(); j++ ) {
       dr = (*xverts)[(j+1)%xverts->size()] - (*xverts)[j];
       lmax = MAX(lmax,dr.norm());
     }
     #elif defined(_G_IS3D)
     for ( GSIZET j=0; j<4; j++ ) { // bottom
       dr = (*xverts)[(j+1)%xverts->size()] - (*xverts)[j];
       lmax = MAX(lmax,dr.norm());
     }
     for ( GSIZET j=4; j<8; j++ ) { // top
       dr = (*xverts)[(j+1)%xverts->size()] - (*xverts)[j];
       lmax = MAX(lmax,dr.norm());
     }
     for ( GSIZET j=0; j<4; j++ ) { // vertical edges
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
// METHOD : grid_init
// DESC   : Initialize global (metric) variables. All elements are assumed to be
//          of the same type.
// ARGS   : none
// RETURNS: none
//**********************************************************************************
void GGrid::grid_init()
{

  // Restrict grid to a single element type:
  assert(ntype_.multiplicity(0) == GE_MAX-1
        && "Only a single element type allowed on grid");


  // Have elements been set yet?
  assert(gelems_.size() > 0 && "Elements not set");

  if      ( itype_[GE_2DEMBEDDED].size() > 0 ) gtype_ = GE_2DEMBEDDED;
  else if ( itype_[GE_DEFORMED]  .size() > 0 ) gtype_ = GE_DEFORMED;
  else if ( itype_[GE_REGULAR]   .size() > 0 ) gtype_ = GE_REGULAR;

  if ( itype_[GE_2DEMBEDDED].size() > 0
    || itype_  [GE_DEFORMED].size() > 0 ) {
    def_init();
  }

  if ( itype_[GE_REGULAR].size() > 0 ) {
    reg_init();
  }

  find_min_dist();
  init_bc_info();

  bInitialized_ = TRUE;

} // end of method grid_init


//**********************************************************************************
//**********************************************************************************
// METHOD : def_init
// DESC   : Initialize global (metric) variables for deformed elems. 
//          All elements are assumed to be of the same type.
// ARGS   : none
// RETURNS: none
//**********************************************************************************
void GGrid::def_init()
{
   assert(gelems_.size() > 0 && "Elements not set");

   GSIZET nxy = itype_[GE_2DEMBEDDED].size() > 0 ? GDIM+1 : GDIM;
   GTMatrix<GTVector<GFTYPE>> rijtmp;

   // Resize geometric quantities to global size:
   dXidX_.resize(nxy,nxy);
   rijtmp.resize(nxy,nxy);
   for ( GSIZET j=0; j<nxy; j++ ) {
     for ( GSIZET i=0; i<nxy; i++ )  {
       dXidX_(i,j).resize(ndof());
     }
   }
   Jac_.resize(ndof());
   faceJac_.resize(nsurfdof());

   // Resize surface-point-wise normals:
   faceNormal_.resize(nxy); // no. coords for each normal at each face point
   for ( GSIZET i=0; i<faceNormal_.size(); i++ ) faceNormal_[i].resize(nsurfdof());

   // Now, set the geometry/metric quanties from the elements:
   GSIZET nfnode; // number of face nodes
   GSIZET ibeg, iend; // beg, end indices for global arrays
   GSIZET ifbeg, ifend; // beg, end indices for global arrays for face quantities
   for ( GSIZET e=0; e<gelems_.size(); e++ ) {
     ibeg  = gelems_[e]->igbeg(); iend  = gelems_[e]->igend();
     ifbeg = gelems_[e]->ifbeg(); ifend = gelems_[e]->ifend();

     // Restrict global data to local scope:
     for ( GSIZET j=0; j<nxy; j++ ) {
       faceNormal_[j].range(ifbeg, ifend); // set range for each coord, j
       for ( GSIZET i=0; i<nxy; i++ )  {
         dXidX_(i,j).range(ibeg, iend);
       }
     }
     Jac_.range(ibeg, iend);
     faceJac_.range(ifbeg, ifend);

     // Set the geom/metric quantities:
     if ( GDIM == 2 ) {
       gelems_[e]->dogeom2d(rijtmp, dXidX_, Jac_, faceJac_, faceNormal_);
     } else if ( GDIM == 3 ) {
       gelems_[e]->dogeom3d(rijtmp, dXidX_, Jac_, faceJac_, faceNormal_);
     }
     
   } // end, element loop

   // Reset global scope:
   for ( GSIZET j=0; j<faceNormal_.size(); j++ ) {
     faceNormal_[j].range_reset();
   }
   for ( GSIZET j=0; j<dXidX_.size(2); j++ )  {
     for ( GSIZET i=0; i<dXidX_.size(1); i++ )  {
       dXidX_(i,j).range_reset();
     }
   }
   Jac_.range_reset();
   
} // end of method def_init


//**********************************************************************************
//**********************************************************************************
// METHOD : reg_init
// DESC   : Initialize global (metric) variables for regular elemes. 
//          All elements are assumed to be
//          of the same type.
// ARGS   : none
// RETURNS: none
//**********************************************************************************
void GGrid::reg_init()
{
   assert(gelems_.size() > 0 && "Elements not set");

   GSIZET nxy = GDIM;
   GTMatrix<GTVector<GFTYPE>>  rijtmp;
   GTVector<GTVector<GFTYPE>> *xe;

   // Resize geometric quantities to global size:
   dXidX_.resize(nxy,1);
   rijtmp.resize(nxy,1);
   for ( GSIZET i=0; i<nxy; i++ ) {
     dXidX_(i,0).resize(ndof());
   }
   Jac_.resize(ndof());
   faceJac_.resize(nsurfdof());

   xNodes_.resize(nxy);
   for ( GSIZET j=0; j<nxy; j++ ) xNodes_[j].resize(ndof());

   // Resize surface-point-wise normals:
   faceNormal_.resize(nxy); // no. coords for each normal at each face point
   for ( GSIZET i=0; i<faceNormal_.size(); i++ ) faceNormal_[i].resize(nsurfdof());

   // Now, set the geometry/metric quanties from the elements:
   GSIZET nfnode; // number of face nodes
   GSIZET ibeg, iend; // beg, end indices for global arrays
   GSIZET ifbeg, ifend; // beg, end indices for global arrays for face quantities
   for ( GSIZET e=0; e<gelems_.size(); e++ ) {
     ibeg  = gelems_[e]->igbeg(); iend  = gelems_[e]->igend();
     ifbeg = gelems_[e]->ifbeg(); ifend = gelems_[e]->ifend();
     xe    = &gelems_[e]->xNodes();
   
     // Restrict global data to local scope:
     for ( GSIZET j=0; j<nxy; j++ ) {
       faceNormal_[j].range(ifbeg, ifend); // set range for each coord, j
       for ( GSIZET i=0; i<nxy; i++ )  {
         dXidX_(i,j).range(ibeg, iend);
       }
     }
     Jac_.range(ibeg, iend);
     faceJac_.range(ifbeg, ifend);

     // Set nodal Cart coords
     for ( GSIZET j=0; j<nxy; j++ ) {
       xNodes_[j].range(ibeg, iend);
       xNodes_[j] = (*xe)[j];
        // 0-out local xNodes; only global allowed now:
       (*xe)[j].clear(); 
     }

     // Set the geom/metric quantities:
     if ( GDIM == 2 ) {
       gelems_[e]->dogeom2d(rijtmp, dXidX_, Jac_, faceJac_, faceNormal_);
     } else if ( GDIM == 3 ) {
       gelems_[e]->dogeom3d(rijtmp, dXidX_, Jac_, faceJac_, faceNormal_);
     }
    
      
   } // end, element loop

   // Reset global scope:
   for ( GSIZET j=0; j<faceNormal_.size(); j++ ) {
     faceNormal_[j].range_reset();
   }
   for ( GSIZET j=0; j<dXidX_.size(2); j++ )  {
     for ( GSIZET i=0; i<dXidX_.size(1); i++ )  {
       dXidX_(i,j).range_reset();
     }
   }
   Jac_.range_reset();
   for ( GSIZET j=0; j<nxy; j++ ) xNodes_[j].range_reset();
   
} // end of method reg_init


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
// METHOD : faceNormal
// DESC   : Return global vector of normals at face nodes
// ARGS   : none
// RETURNS: GTVector<GTVector<GFTYPE>> &
//**********************************************************************************
GTVector<GTVector<GFTYPE>> &GGrid::faceNormal()
{
   assert(bInitialized_ && "Object not inititaized");
   return faceNormal_;

} // end of method faceNormal

//**********************************************************************************
//**********************************************************************************
// METHOD : find_min_dist
// DESC   : Compute min node distance for each element
// ARGS   : none
// RETURNS: none
//**********************************************************************************
void GGrid::find_min_dist()
{
  assert(gelems_.size() > 0 && "Elements not set");

  
 
  GSIZET            Nx, Ny, Nz, nxy;
  GFTYPE            del, dr;
  GTVector<GINT>   *N, I(7);
  GTPoint<GFTYPE>  dx(3), p0(3), p1(3);
  GTVector<GTVector<GFTYPE>> *xn;

  // Find min node distance for each elem. 
  dx[2] = 0.0;
  p0[2] = 0.0;
  p1[2] = 0.0;

  minnodedist_.resize(gelems_.size()); // one dr for each elem

  // For GE_2DEMBEDDED grids:
  if ( itype_[GE_2DEMBEDDED].size() > 0 ) {
    for ( auto e=0; e<gelems_.size(); e++ ) {
      dr = std::numeric_limits<GFTYPE>::max();
      xn = &gelems_[e]->xNodes(); 
      N  = &gelems_[e]->dim(); Nx = (*N)[0]; Ny = (*N)[1];
      for ( auto j=0; j<Ny; j++ ) { // loop over all subcells
        for ( auto i=0; i<Nx; i++ ) {
          I[0] = i+1+j*Nx; I[1] = i+j*Nx; I[2] = i+1+(j+1)*Nx;
          for ( auto n=0; n<GDIM; n++ ) p0[n] = (*xn)[n][i+j*Nx];
          for ( auto l=0; l<GDIM*(GDIM-1)+1; l++ ) {
            for ( auto n=0; n<GDIM+1; n++ ) p1[n] = (*xn)[n][I[l]];
            dx = p1 - p0;
            del   = dx.norm();
            if ( del > 0.0 ) dr = MIN(dr,del);
          } // end, cell vertices
        }
      }
      minnodedist_[e] = dr;
    } // end, elem loop
  } 

  // For GE_DEFROMED grids:
  if ( itype_[GE_DEFORMED].size() > 0  ) {
    if ( GDIM == 2 ) {
      for ( auto e=0; e<gelems_.size(); e++ ) {
        dr = std::numeric_limits<GFTYPE>::max();
        xn = &gelems_[e]->xNodes(); 
        N  = &gelems_[e]->dim(); Nx = (*N)[0]; Ny = (*N)[1];
        for ( auto j=0; j<Ny; j++ ) { // loop over all subcells
          for ( auto i=0; i<Nx; i++ ) {
            I[0] = i+1+j*Nx; I[1] = i+j*Nx; I[2] = i+1+(j+1)*Nx;
            for ( auto n=0; n<GDIM; n++ ) p0[n] = (*xn)[n][i+j*Nx];
            for ( auto l=0; l<GDIM*(GDIM-1)+1; l++ ) {
              for ( auto n=0; n<GDIM; n++ ) p1[n] = (*xn)[n][I[l]];
              dx = p1 - p0;
              del= dx.norm();
              if ( del > 0.0 ) dr = MIN(dr,del);
            } // end, cell vertices
          }
        }
        minnodedist_[e] = dr;
      } // end, elem loop
    }
    else if ( GDIM == 3 ) {
      for ( auto e=0; e<gelems_.size(); e++ ) {
        dr = std::numeric_limits<GFTYPE>::max();
        xn = &gelems_[e]->xNodes(); 
        N  = &gelems_[e]->dim();
        Nx = (*N)[0]; Ny = (*N)[1]; Nz = (*N)[2]; nxy = Nx*Ny;
        for ( auto k=0; k<Nz-1; k++ ) { // loop over all subcells
          for ( auto j=0; j<Ny-1; j++ ) {
            for ( auto i=0; i<Nx-1; i++ ) {
              I[0] = i+1+j*Nx+k*nxy; I[1] = i+j*Nx+k*nxy; I[2] = i+1+(j+1)*Nx+k*nxy;
              I[3] = i+1+j*Nx+k*nxy; I[4] = i+j*Nx+k*nxy; I[5] = i+1+(j+1)*Nx+k*nxy; I[6] = i+j*Nx+k*nxy;
              for ( auto n=0; n<GDIM; n++ ) p0[n] = (*xn)[n][i+j*Nx];
              for ( auto l=0; l<GDIM*(GDIM-1)+1; l++ ) {
                for ( auto n=0; n<GDIM; n++ ) p1[n] = (*xn)[n][I[l]];
                dx = p1 - p0;
                del= dx.norm();
                if ( del > 0.0 ) dr = MIN(dr,del);
              } // end, cell vertices, l loop
            } // end, i loop
          } // end, j loop
        } // end, k loop
        minnodedist_[e] = dr;
      } // end, element loop
    } // end of GDIM == 3 test
  } // end/ of GE_DEFORMED test

  // For GE_REGULAR grids:
  if ( itype_[GE_REGULAR].size() > 0 ) {
    for ( auto e=0; e<gelems_.size(); e++ ) {
      dr = std::numeric_limits<GFTYPE>::max();
      xn = &gelems_[e]->xNodes(); 
      N  = &gelems_[e]->dim();
      for ( auto j=0; j<(*N)[0]-1; j++ ) { // x-direction
        del = fabs((*xn)[0][j+1] - (*xn)[0][j]);
        if ( del > 0.0 ) dr = MIN(dr,del);
      }
      for ( auto j=0; j<(*N)[1]-1; j++ ) { // y-direction
        del = fabs((*xn)[1][(j+1)*(*N)[0]] - (*xn)[1][j*(*N)[0]]);
        if ( del > 0.0 ) dr = MIN(dr,del);
      }
      if ( GDIM > 2 ) { // z-direction
        nxy = (*N)[0] * (*N)[1];
        for ( auto j=0; j<(*N)[1]-1; j++ ) {
          del = fabs((*xn)[2][(j+1)*nxy] - (*xn)[2][j*nxy]);
          if ( del > 0.0 ) dr = MIN(dr,del);
        }
      }
      minnodedist_[e] = dr;
    }
  } // end/ of GE_REGULAR test

} // end of method find_min_dist


//**********************************************************************************
//**********************************************************************************
// METHOD : integrate
// DESC   : Compute spatial integral of global input vector. Result 
//          is a sum over all MPI tasks. NOTE: This method extracts weights
//          from the elements, so don't use this method very often.
// ARGS   : u  : 'global' integral argument
//          tmp: tmp vector, same size as u
// RETURNS: GFTYPE integral
//**********************************************************************************
GFTYPE GGrid::integrate(GTVector<GFTYPE> &u, GTVector<GFTYPE> &tmp)
{
  assert(bInitialized_ && "Object not inititaized");

  GSIZET                       ibeg, iend; // beg, end indices for global array
  GSIZET                       n;
  GFTYPE                       xint, xgint;
  GTVector<GINT>               N(GDIM);    // coord node sizes
  GTVector<GTVector<GFTYPE>*>  W(GDIM);    // element weights

  tmp = u;
#if defined(_G_IS2D)
  for ( GSIZET e=0; e<gelems_.size(); e++ ) {
    ibeg  = gelems_[e]->igbeg(); iend  = gelems_[e]->igend();

    for ( GSIZET k=0; k<GDIM; k++ ) {
      W[k] = gelems_[e]->gbasis(k)->getWeights();
      N[k] = gelems_[e]->size(k);
    }
    n = 0;
    for ( GSIZET k=0; k<N[1]; k++ ) {
      for ( GSIZET j=0; j<N[0]; j++,n++ ) {
        tmp[n] = (*W[1])[k]*(*W[0])[j] * u[n];
      }
    }
    // Restrict global data to local scope:
    tmp.range(ibeg, iend);
  } // end, element loop
#elif defined(_G_IS3D)
  for ( GSIZET e=0; e<gelems_.size(); e++ ) {
    ibeg  = gelems_[e]->igbeg(); iend  = gelems_[e]->igend();

    for ( GSIZET k=0; k<GDIM; k++ ) {
      W[k] = gelems_[e]->gbasis(k)->getWeights();
      N[k] = gelems_[e]->size(k);
    }
    n = 0;
    for ( GSIZET k=0; k<N[2]; k++ ) {
      for ( GSIZET j=0; j<N[1]; j++ ) {
        for ( GSIZET i=0; i<N[0]; i++,n++ ) {
          tmp[n] = (*W[2])[k]*(*W[1])[k]*(*W[0])[j] * u[n];
        }
      }
    }
    // Restrict global data to local scope:
    tmp.range(ibeg, iend);
  } // end, element loop
#endif
  tmp.range_reset(); 

  // Multiply by Jacobian:
  tmp.pointProd(Jac_);
  xint = tmp.sum();  
  GComm::Allreduce(&xint, &xgint, 1, T2GCDatatype<GFTYPE>() , GC_OP_SUM, comm_);

  return xgint;

} // end of method integrate


//**********************************************************************************
//**********************************************************************************
// METHOD : init_bc_info
// DESC   : Set global bdy condition data from the element bdy data.
// ARGS   : none
// RETURNS: none.
//**********************************************************************************
void GGrid::init_bc_info()
{
  GSIZET                       ibeg, iend; // beg, end indices for global array
  GTVector<GINT>              *iebdy;  // domain bdy indices
  GTVector<GBdyType>          *iebdyt; // domain bdy types

  // Collect all element bdy types and indicection indices
  // into global vectors (so we can use the vector 
  // to do sorting):
  GSIZET              nn=0; 
  GSIZET                ig; 
  GTVector <GBdyType> btmp; 
  GTVector   <GSIZET> itmp; 
  for ( GSIZET e=0; e<gelems_.size(); e++ ) {
    ibeg   = gelems_[e]->igbeg(); iend  = gelems_[e]->igend();
    iebdy  = &gelems_[e]->bdy_indices(); 
    iebdyt = &gelems_[e]->bdy_types(); 
    for ( GSIZET j=0; j<iebdyt->size(); j++ ) {
      ig = nn + (*iebdy)[j];
      itmp.push_back(ig); // index in full global arrays
      btmp.push_back((*iebdyt)[j]);
      nn += gelems_[e]->nnodes();
    }
  } // end, element loop
 
  // Now, create type-bins (one bin for each GBdyType), and
  // for each type, set the indirection indices into global
  // vectors that have that type:
  GBdyType         itype;
  GSIZET    *ind=NULLPTR;
  GSIZET      nind, nw=0;
  igbdy_.resize(GBDY_NONE); // set of bdy indices for each type
  for ( GSIZET k=0; k<GBDY_NONE; k++ ) { // cycle over each bc type
    itype = static_cast<GBdyType>(k);
    nind = btmp.contains(itype, ind, nw);
    igbdy_[k].resize(nind);
    for ( GSIZET j=0; j<nind; j++ ) igbdy_[k][j] = itmp[ind[j]];
    nind = 0;
  } // end, element loop

  if ( ind != NULLPTR ) delete [] ind;

} // end of method init_bc_info





