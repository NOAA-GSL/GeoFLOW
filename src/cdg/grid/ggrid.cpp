//==================================================================================
// Module       : ggrid_
// Date         : 8/31/18 (DLR)
// Description  : GeoFLOW grid object. Is primarily a container for
//                a list of elements, as an indirection buffer
//                that breaks down elements by type.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include "gelem_base.hpp"
#include "ggrid.hpp"


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Default constructor
// ARGS   : 
// RETURNS: none
//**********************************************************************************
GGrid::GGrid()
:
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
  for ( GSIZET j=0; j<GE_MAX; j++ ) {
    itmp.contains(static_cast<GElemType>(j),ind,nd);
    for ( GSIZET i=0; i<nd; i++ ) itype_[j].push_back(ind[i]);
    ntype_[j] = nd;
  }
  if ( ind != NULLPTR ) delete [] ind;

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
// METHOD : minsep
// DESC   : Find min nodal separation
// ARGS   : none
// RETURNS: GFTYPE separation
//**********************************************************************************
GFTYPE GGrid::minsep()
{
   assert(gelems_.size() > 0 && "Elements not set");

   GFTYPE msep = std::numeric_limits<GFTYPE>::max();
   GFTYPE dr;
   GTVector<GTVector<GFTYPE>> *xnodes;

   for ( GSIZET i=0; i<gelems_.size(); i++ ) {
     xnodes = &gelems_[i]->xNodes();
     if ( xnodes->size() == 1 ) { 

       for ( GSIZET k=0; k<(*xnodes)[0].size()-1; k++ ) {
         dr   = fabs( (*xnodes)[0][k+1]-(*xnodes)[0][k] ); 
         msep = MIN(msep, dr);
       }

     } else if ( xnodes->size() == 2 ) {

       for ( GSIZET l=0; l<(*xnodes)[1].size()-1; l++ ) {
         for ( GSIZET k=0; k<(*xnodes)[0].size()-1; k++ ) {
           dr   = sqrt( pow((*xnodes)[0][k+1]-(*xnodes)[0][k],2) 
                      + pow((*xnodes)[1][l+1]-(*xnodes)[1][l],2) );
           msep = MIN(msep, dr);
         }
       }

     } else if ( xnodes->size() == 3 ) {

       for ( GSIZET m=0; m<(*xnodes)[2].size()-1; m++ ) {
         for ( GSIZET l=0; l<(*xnodes)[1].size()-1; l++ ) {
           for ( GSIZET k=0; k<(*xnodes)[0].size()-1; k++ ) {
             dr   = sqrt( pow((*xnodes)[0][k+1]-(*xnodes)[0][k],2) 
                        + pow((*xnodes)[1][l+1]-(*xnodes)[1][l],2)
                        + pow((*xnodes)[2][m+1]-(*xnodes)[2][m],2) );
             msep = MIN(msep, dr);
           }
         }
       }

     }
   }

   return msep;
} // end of method minsep


//**********************************************************************************
//**********************************************************************************
// METHOD : maxsep
// DESC   : Find max nodal separation
// ARGS   : none
// RETURNS: GFTYPE separation
//**********************************************************************************
GFTYPE GGrid::maxsep()
{
   assert(gelems_.size() > 0 && "Elements not set");

   GFTYPE msep = 0.0;
   GFTYPE dr;
   GTVector<GTVector<GFTYPE>> *xnodes;

   for ( GSIZET i=0; i<gelems_.size(); i++ ) {
     xnodes = &gelems_[i]->xNodes();
     if ( xnodes->size() == 1 ) { 

       for ( GSIZET k=0; k<(*xnodes)[0].size()-1; k++ ) {
         dr   = fabs( (*xnodes)[0][k+1]-(*xnodes)[0][k] ); 
         msep = MAX(msep, dr);
       }

     } else if ( xnodes->size() == 2 ) {

       for ( GSIZET l=0; l<(*xnodes)[1].size()-1; l++ ) {
         for ( GSIZET k=0; k<(*xnodes)[0].size()-1; k++ ) {
           dr   = sqrt( pow((*xnodes)[0][k+1]-(*xnodes)[0][k],2) 
                      + pow((*xnodes)[1][l+1]-(*xnodes)[1][l],2) );
           msep = MAX(msep, dr);
         }
       }

     } else if ( xnodes->size() == 3 ) {

       for ( GSIZET m=0; m<(*xnodes)[2].size()-1; m++ ) {
         for ( GSIZET l=0; l<(*xnodes)[1].size()-1; l++ ) {
           for ( GSIZET k=0; k<(*xnodes)[0].size()-1; k++ ) {
             dr   = sqrt( pow((*xnodes)[0][k+1]-(*xnodes)[0][k],2) 
                        + pow((*xnodes)[1][l+1]-(*xnodes)[1][l],2)
                        + pow((*xnodes)[2][m+1]-(*xnodes)[2][m],2) );
             msep = MAX(msep, dr);
           }
         }
       }

     }
   }

   return msep;
} // end of method maxsep


//**********************************************************************************
//**********************************************************************************
// METHOD : init
// DESC   : Initialize global (metric) variables. All elements are assumed to be
//          of the same type.
// ARGS   : none
// RETURNS: none
//**********************************************************************************
void GGrid::init()
{

   assert(gelems_.size() > 0 && "Elements not set");

  // Restrict grid to a single element type:
  assert(grid_->ntype().multiplicity(0) == GE_MAX-1
        && "Only a single element type allowed on grid");

  if ( grid_->itype(GE_2DEMBEDDED).size() > 0
    || grid_->itype  (GE_DEFORMED).size() > 0 ) {
    def_init();
  }

  if ( grid_->itype(GE_REGULAR).size() > 0 ) {
    reg_init();
  }


} // end of method init


//**********************************************************************************
//**********************************************************************************
// METHOD : def_init
// DESC   : Initialize global (metric) variables. All elements are assumed to be
//          of the same type.
// ARGS   : none
// RETURNS: none
//**********************************************************************************
void GGrid::def_init()
{
   assert(gelems_.size() > 0 && "Elements not set");

   GSIZET nxy = grid_->itype(GE_2DEMBEDDED) > 0 > GDIM+1 : GDIM;
   GMatrix<GTVector<GFTYPE>> rijtmp;

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
   faceNormal.resize(nxy); // no. coords for each normal at each face point
   for ( GSIZET i=0; i<faceNormal.size(); i++ ) faceNormal[i].resize(nsurfdof());

   // Now, set the geometry/metric quanties from the elements:
   GSIZET nfnode; // number of face nodes
   GSIZET ibeg, iend; // beg, end indices for global arrays
   GSIZET ifbeg, ifend; // beg, end indices for global arrays for face quantities
   for ( GSIZET i=0; i<gelems_.size(); i++ ) {
     ibeg  = (*gelems)[e]->igbeg(); iend  = (*gelems)[e]->igend();
     ifbeg = (*gelems)[e]->ifbeg(); ifend = (*gelems)[e]->ifend();

     // Restrict global data to local scope:
     for ( GSIZET j=0; j<nxy; j++ ) {
       faceNormal[j].range(ifbeg, ifend); // set range for each coord, j
       for ( GSIZET i=0; i<nxy; i++ )  {
         dXidX_(i,j).range(ibeg, iend);
       }
     }
     Jac_.range(ibeg, iend);
     faceJac_.range(ifbeg, ifend);

     // Set the geom/metric quantities:
     if ( GDIM == 2 ) {
       gelems_[i]->dogeom2d(rijtmp, dXidX_, Jac_, faceJac_, faceNormal_);
     } else if ( GDIM == 3 ) {
       gelems_[i]->dogeom3d(rijtmp, dXidX_, Jac_, faceJac_, faceNormal_);
     }
     
   } // end, element loop

   // Reset global scope:
   ibeg  = 0; iend  = ndof()-1;;
   ifbeg = 0; ifend = nsurfdof()-1;
   for ( GSIZET j=0; j<nxy; j++ ) {
     faceNormal[j].range(ifbeg, ifend);
     for ( GSIZET i=0; i<nxy; i++ )  {
       dXidX_(i,j).range(ibeg, iend);
     }
   }
   Jac_.range(ibeg, iend);
   
} // end of method def_init


