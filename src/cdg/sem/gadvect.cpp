//==================================================================================
// Module       : gadvect.hpp
// Date         : 11/11/18 (DLR)
// Description  : Represents the SEM discretization of the advection operator:
//                u.Grad p  This is a nonlinear operator, so should not derive 
//                from GLinOp. This operator requires that grid consist of
//                elements of only one type.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none
//==================================================================================

#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include "gadvect.hpp"
#include "gtmatrix.hpp"
#include "gmtk.hpp"
#include "gmass.hpp"

using namespace std;

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Constructor accepting grid operator only
// ARGS   : grid  : GGrid object
// RETURNS: none
//**********************************************************************************
GAdvect::GAdvect(GGrid &grid)
:
bInitialized_   (FALSE),
grid_           (&grid)
{
  init();
} // end of constructor method (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
GAdvect::~GAdvect()
{
} // end, destructor


//**********************************************************************************
//**********************************************************************************
// METHOD : apply
// DESC   : Compute application of this operator to input vector:
//            po = u.Grad p
//          NOTE: Require GDIM velocity components for GE_DEFORMED or
//                GE_REGULAR elements, and require GDIM+1 components
//                for GE_2DEMBEDDED elements.
// ARGS   : p   : input p field
//          u   : input vector field
//          utmp: tmp arrays
//          po  : output (result) vector
//             
// RETURNS:  none
//**********************************************************************************
void GAdvect::apply(GTVector<GFTYPE> &p, const GTVector<GTVector<GFTYPE>*> &u, GTVector<GTVector<GFTYPE>*> &utmp, GTVector<GFTYPE> &po) 
{
  assert(bInitialized_ && "Operator not initialized");
    
  GSIZET nxy = grid_->gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;

  assert(u.size() >= nxy && "Insufficient number of velocity components");

  if ( grid_->gtype() == GE_REGULAR ) {
    reg_prod(p, u, utmp, po);
  }
  else if ( grid_->gtype() == GE_DEFORMED 
         || grid_->gtype() == GE_2DEMBEDDED ) {
    def_prod(p, u, utmp, po);
  }

} // end of method apply


//**********************************************************************************
//**********************************************************************************
// METHOD : def_prod
// DESC   : Compute this term for 2d & 3d GE_DEFORMED elements.
//          NOTE: must have 5 utmp vectors set via call if
//                using GDIM=3 or if GDIM=2 & are using GE_2DEMBEDDED elements;
//                If using 2D deformed elements, need just 4 utmp vectors.
// ARGS   : p : input vector (global)
//          u : input vector field (global)
//          po: output (result) field (global)
//             
// RETURNS:  none
//**********************************************************************************
void GAdvect::def_prod(GTVector<GFTYPE> &p, const GTVector<GTVector<GFTYPE>*> &u, GTVector<GTVector<GFTYPE>*> &utmp, GTVector<GFTYPE> &po) 
{
  assert( utmp.size() >= GDIM+1
       && "Insufficient temp space specified");


// Must compute:
//    po = u . Grad p = Dj Rji (pu^i) W J  
// Can be cast into:
//                      --              -- --  --
//                      | R11   R12  R13 | | D1 |
//       = [u1 u2 u3]   | R21   R22  R23 | | D2 | p W J
//                      | R31   R32  R33 | | D3 |
//                      --              -- --  --
//
// where
//    D1 f = (I_X_I_X_Dx)f; D2 f = (I_X_Dy_X_I)f; D3 f = (Dz_X_I_X_I)f
// and Rij derivative matrix computed in the element, dX^j/dX^i, which
// doesn't include the weights or the Jacobian.

  GTMatrix<GTVector<GFTYPE>> *dXidX = &grid_->dXidX(); // get Rij
  GTVector<GFTYPE>           *Jac   = &grid_->Jac  (); // get J

  // Get derivatives with weights:
  GMTK::compute_grefderivsW(*grid_, p, etmp1_, FALSE, utmp); // utmp stores tensor-prod derivatives, Dj p

  GSIZET nxy = grid_->gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;

  // Compute po += ui Rij (D^j p): 
  po = 0.0;
  for ( GSIZET j=0; j<nxy; j++ ) { 
    *utmp[nxy+1]=0.0;
    for ( GSIZET i=0; i<nxy; i++ ) {
      (*dXidX)(i,j).pointProd(*utmp[i],*utmp[nxy]); // Rij Dj p
      *utmp[nxy+1] += *utmp[nxy];
    }
    utmp[nxy+1]->pointProd(*Jac); // J Rij  Dj p
    utmp[nxy+1]->pointProd(*u[j]); // do uj J Rij Dj p
    po += *utmp[nxy+1];
  }

} // end of method def_prod


//**********************************************************************************
//**********************************************************************************
// METHOD : reg_prod
// DESC   : Compute application of this operator to input vector for 
//          GE_REGULAR elements
//          NOTE: must have GDIM utmp vectors allocated
// ARGS   : p   : input vector (global)
//          u   : input vector field (global)
//          utmp: temp arrays
//          po  : output (result) field (global)
//             
// RETURNS:  none
//**********************************************************************************
void GAdvect::reg_prod(GTVector<GFTYPE> &p, const GTVector<GTVector<GFTYPE>*> &u, GTVector<GTVector<GFTYPE>*> &utmp, GTVector<GFTYPE> &po) 
{

  assert( utmp.size() >= GDIM
       && "Insufficient temp space specified");

// Must compute:
//    po = u . Grad p => W J (G1 u1 D1 + G2 u2 D2 + G3 u3 D3 ) p
//
// where
//    D1u = (I_X_I_X_Dx)u; D2u = (I_X_Dy_X_I)u; D3u = (Dz_X_I_X_I)u
// and Gj are the 'metric' terms computed for the grid, dXi^j/dX^j.
// For regular elements, these are diagnonal, and we include in G_j
// Jacobians, and weights are inclued in the computation of 
// the derivative:

  // Get reference derivatives:
  GMTK::compute_grefderivsW(*grid_, p, etmp1_, FALSE, utmp); // utmp stores tensor-prod derivatives, Dj p

  // Compute po += Gj uj D^j p): 
  po = *utmp[0];
  po.pointProd(*G_[0]);// remember, mass & Jac included in G
  po.pointProd(*u[0]); // do uj * (Gj * Dj p)
  for ( GSIZET j=1; j<GDIM; j++ ) { 
    utmp [j]->pointProd(*G_[j]);// remember, mass & Jac included in G
    utmp [j]->pointProd(*u[j]); // do uj * (Gj * Dj p)
    po += *utmp[j];
  }

} // end of method reg_prod


//**********************************************************************************
//**********************************************************************************
// METHOD : init
// DESC   : Compute metric components for all elements. Split up this way
//          in case additional data must be initialized.
// ARGS   : none
// RETURNS: none
//**********************************************************************************
void GAdvect::init()
{
  assert(grid_->ntype().multiplicity(0) == GE_MAX-1 
        && "Only a single element type allowed on grid");

  if ( grid_->gtype() == GE_2DEMBEDDED 
    || grid_->gtype() == GE_DEFORMED ) {
    def_init();
  }
  if ( grid_->gtype() == GE_REGULAR ) {
    reg_init();
  }

  bInitialized_ = TRUE;

} // end, method init


//**********************************************************************************
//**********************************************************************************
// METHOD : def_init
// DESC   : Compute metric components if there are deformed elements.
// ARGS   : none
// RETURNS: none
//**********************************************************************************
void GAdvect::def_init()
{

  if ( grid_->gtype() == GE_2DEMBEDDED 
    || grid_->gtype() == GE_DEFORMED ) return;


  GTVector<GSIZET>             N(GDIM);
  GTMatrix<GTVector<GFTYPE>>  *dXidX;    // element-based dXi/dX matrix
  GTVector<GTVector<GFTYPE>*>  W(GDIM);  // element-based weights
  GTVector<GFTYPE>            *Jac;      // element-based Jacobian
  GElemList                   *gelems = &grid_->elems();

  // Compute 'metric' components:
  // Gi = [dxi/dx, deta/dy, dzeta/dz]; 
  GSIZET nxy = grid_->gtype() == GE_2DEMBEDDED ? GDIM+1: GDIM;
  GSIZET ibeg, iend; // beg, end indices for global array
  G_ .resize(nxy);
  G_ = NULLPTR;
  for ( GSIZET j=0; j<nxy; j++ ) {
    G_ [j] = new GTVector<GFTYPE>(grid_->ndof());
  }

  Jac   = &grid_->Jac();
  dXidX = &grid_->dXidX();

  // Cycle through all elements; fill metric elements
  for ( GSIZET e=0; e<grid_->elems().size(); e++ ) {
    if ( (*gelems)[e]->elemtype() != GE_DEFORMED 
      && (*gelems)[e]->elemtype() != GE_2DEMBEDDED ) continue;

    ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
    for ( GSIZET j=0; j<GDIM; j++ ) {
      W[j]= (*gelems)[e]->gbasis(j)->getWeights();
      N[j]= (*gelems)[e]->size(j);
    }
    Jac->range(ibeg, iend); // Restrict range of global vectors:
    for ( GSIZET j=0; j<nxy; j++ )
      for ( GSIZET i=0; i<nxy; i++ ) (*dXidX)(i,j).range(ibeg, iend);

#if defined(_G_IS2D)

    for ( GSIZET j=0; j<nxy; j++ ) { // G vector element 
      (*G_[j]).range(ibeg, iend); // restrict global vec to local range
      for ( GSIZET m=0, n=0; m<N[1]; m++ ) {
        for ( GSIZET l=0; l<N[0]; l++,n++ ) {
          (*G_[j])[n] = (*dXidX)(j,j)[n] 
                      * (*W[0])[l] * (*W[1])[m] * (*Jac)[n];
        }
      }
      (*G_[j]).range_reset(); // reset to global range
    }

#else

    for ( GSIZET j=0; j<nxy; j++ ) { // G vector element 
      (*G_[j]).range(ibeg, iend); // restrict global vec to local range
        for ( GSIZET p=0, n=0; p<N[2]; p++ ) {
          for ( GSIZET m=0; m<N[1]; m++ ) {
            for ( GSIZET l=0; l<N[0]; l++,n++ ) {
              (*G_[j])[n] = (*dXidX)(j,j)[n] 
                          * (*W[0])[l] * (*W[1])[m] * (*W[2])[p] * (*Jac)[n];
            }
          }
        }
      (*G_[j]).range_reset(); // reset global vec to global range
    }

#endif
  } // end, element list

  // Reset ranges to global scope:
  Jac->range_reset();
  for ( GSIZET j=0; j<nxy; j++ )
    for ( GSIZET i=0; i<nxy; i++ ) (*dXidX)(i,j).range_reset();

} // end of method def_init


//**********************************************************************************
//**********************************************************************************
// METHOD : reg_init
// DESC   : Compute metric components for regular elements.
// ARGS   : none
// RETURNS: none
//**********************************************************************************
void GAdvect::reg_init()
{
  if ( grid_->gtype() != GE_REGULAR ) return; 


  GTVector<GSIZET>             N(GDIM);
  GTMatrix<GTVector<GFTYPE>>  *dXidX;    // dXi/dX matrix
  GTVector<GFTYPE>            *Jac;      // Jacobian
//GMass                        mass(*grid_);

  // Compute 'metric' components:
  // Gi = [dxi/dx, deta/dy, dzeta/dz]; 
  //
  // We have
  // 
  //    Jac = L1 L2 L3/8 , 3d
  //    Jac = L1 L2 /4 , 2d
  //
  //
  //  G1 = dXi^1/dX^1 = 2 / L1
  //  G2 = dXi^2/dX^2 = 2 / L2
  //  G3 = dXi^3/dX^3 = 2 / L3
  // These values should already be set in grid object
  Jac   = &grid_->Jac();
  dXidX = &grid_->dXidX();

  GSIZET nxy = GDIM;
  GSIZET ibeg, iend; // beg, end indices for global array
  G_ .resize(nxy);
  G_ = NULLPTR;
  for ( GSIZET j=0; j<nxy; j++ ) {
    G_ [j] = new GTVector<GFTYPE>(grid_->ndof());
  }

  // Set G_, and include Jacobian & mass:
  for ( GSIZET j=0; j<GDIM; j++ ) {
   *G_[j] = (*dXidX)(j,0);
    G_[j]->pointProd(*Jac);
//  G_[j]->pointProd(*(mass.data()));
  }


} // end of method reg_init


