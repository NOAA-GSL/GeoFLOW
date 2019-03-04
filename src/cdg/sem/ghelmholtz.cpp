//==================================================================================
// Module       : glaplacian.hpp
// Date         : 11/1/18 (DLR)
// Description  : Represents the SEM generalized Helmholtz operator:
//                  H = qM + pL,
//                where M = mass operator, L is Laplacian operator, and
//                q, p are scalars that may or may not be constant. 
//                The mass term is added only if the mass operator is set
//                via call to set_mass, and the scalars applied only
//                if calls are made to set_Lap_scalar and set_mass_scalar,
//                respectively. In fact, only if the mass operator is set
//                is this operator a real Helmholtz operator; otherwise, it's
//                really just a weak Laplacian operator. 
//                Note: this operator will fail if the grid contains more 
//                      than a single element type.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : GLinOp
//==================================================================================

#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include "ghelmholtz.hpp"
#include "gtmatrix.hpp"
#include "gmtk.hpp"

using namespace std;

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Default constructor
// ARGS   : none
// RETURNS: none
//**********************************************************************************
GHelmholtz::GHelmholtz(GGrid &grid)
: GLinOp(grid),
bown_p_        (TRUE),
bown_q_       (FALSE),
bown_mass_    (FALSE),
bcompute_helm_(FALSE),
p_          (NULLPTR),
q_          (NULLPTR),
massop_     (NULLPTR)
{
  grid_ = &grid;
  p_ = new GTVector<GFTYPE>(1);
 *p_ =  1.0; // diffusion defaults to a scalar, with value 1
  init();
} // end of constructor method (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
GHelmholtz::~GHelmholtz()
{
  for ( GSIZET j=0; j<G_.size(2); j++ ) {
    for ( GSIZET i=j; i<G_.size(1); i++ ) {
      if ( G_ (i,j) != NULLPTR ) delete G_(i,j);
    }
  }
  if ( bown_q_ && q_ != NULLPTR ) delete q_;
  if ( bown_p_ && p_ != NULLPTR ) delete p_;
  if ( bown_mass_ && massop_ != NULLPTR ) delete massop_;

} // end, destructor


//**********************************************************************************
//**********************************************************************************
// METHOD : opVec_prod
// DESC   : Compute application of this operator to input vector.
// ARGS   : u   : input field (component)
//          uo  : output (result) field
//          utmp: tmp space
//             
// RETURNS:  none
//**********************************************************************************
void GHelmholtz::opVec_prod(GTVector<GFTYPE> &u, 
                            GTVector<GTVector<GFTYPE>*> &utmp,
                            GTVector<GFTYPE> &uo)
{
  assert(bInitialized_ && "Operator not initialized");
    
  if ( grid_->igtype() == GE_REGULAR) ) {
    reg_prod(u, utmp, uo);
  }
  else if ( grid_->gtype() == GE_DEFORMED ) {
    def_prod(u, utmp, uo);
  }
  if ( grid_->gtype() == GE_2DEMBEDDED ) {
     embed_prod(u, utmp, uo);
  }

} // end of method opVec_prod


//**********************************************************************************
//**********************************************************************************
// METHOD : def_prod
// DESC   : Compute application of this operator to input vector for 
//          for 2d & 3d GE_DEFORMED elements.
//          NOTE: must have 2*GDIM utmp vectors set in utmp
// ARGS   : u : input vector
//          uo: output (result) vector
//          utmp: tmp space
//             
// RETURNS:  none
//**********************************************************************************
void GHelmholtz::def_prod(GTVector<GFTYPE> &u, 
                          GTVector<GTVector<GFTYPE>*> &utmp,
                          GTVector<GFTYPE> &uo)
{
  assert( utmp.size() >= GDIM+3
       && "Insufficient temp space specified");

  GTVector<GTVector<GFTYPE>*> gdu(GDIM);
  GElemList *gelems=&grid_->elems();

  // Must compute:

// uo = ( p L + q M ) u
//
// where 
//                        --              -- --  --
//                        | G11   G12  G13 | | D1 |
// p L u = p [D1 D2 D3]^T | G21   G22  G23 | | D2 | u + q M
//                        | G31   G32  G33 | | D3 |
//                        --              -- --  --
// where
//    D1u = (I_X_I_X_Dx)u; D2u = (I_X_Dy_X_I)u; D3u = (Dz_X_I_X_I)u
// and Gij are the metric terms computed in the element,that include
// the weights and Jacobian; p is 'viscosity', M is mass matrix, 
// and q a factor for M. If GDIM=2, then we remove last column and last row of
// metric matrix.

  // Re-arrange temp space for divergence:
  for ( GSIZET i=0; i<GDIM; i++ ) gdu[i] = utmp[i+GDIM];

  // Compute derivatives of u:
  GMTK::compute_grefderivs(*grid_, u, etmp1_, FALSE, utmp); // utmp stores tensor-prod derivatives, Dj u

  // Compute Gij (D^j u): 
  for ( GSIZET i=0; i<GDIM; i++ ) { 
    *utmp[GDIM+i] = 0.0;
    for ( GSIZET j=0; j<GDIM; j++ ) {
      utmp[j]->pointProd(*G_(i,j),uo); // Gij * du^j
      *utmp[GDIM+i] += uo;
    }
  }

  // Apply p parameter ('viscosity') if necessary to Laplacian:
  if ( p_ != NULLPTR ) {
    if ( p_->size() >= grid_->ndof() ) uo.pointProd(*p_);
    else if ( (*p_)[0] != 1.0 ) uo *= (*p_)[0];
  }

  // utmp[GDIM+1], utmp[GDIM+2', uo now hold Ti = Gij D^j u, i = 0, GDIM-1
  // Now compute DT^j ( t^j ):
  GMTK::compute_grefdiv(*grid_, gdu, etmp1_, TRUE, uo); // Compute 'divergence' with DT_j

  // At this point, we have uo = L u
 

  // Apply Mass operator and q parameter if necessary:
  if ( bcompute_helm_ ) {
    massop_->opVec_prod(u, utmp, *utmp[0]);
    if ( q_ != NULLPTR ) {
      if ( q_->size() >= grid_->ndof() ) utmp[0]->pointProd(*q_);
      else if ( (*q_)[0] != 1.0)        *utmp[0] *= (*q_)[0];
    }
    uo += *utmp[0];
  }

} // end of method def_prod


//**********************************************************************************
//**********************************************************************************
// METHOD : embed_prod
// DESC   : Compute application of this operator to input vector for
//          GE_2DEMBEDDED elements
//          NOTE: must have GDIM+3 utmp vectors set via set_tmp method
// ARGS   : u  : input vector
//          uo : output (result) vector
//          utmp: tmp space
//             
// RETURNS:  none
//**********************************************************************************
void GHelmholtz::embed_prod(GTVector<GFTYPE> &u, 
                            GTVector<GTVector<GFTYPE>*> &utmp,
                            GTVector<GFTYPE> &uo)
{
  assert( GDIM == 2 && utmp.size() >= GDIM+3
       && "Insufficient temp space specified");

  GTVector<GTVector<GFTYPE>*> gdu(GDIM);
  GElemList *gelems=&grid_->elems();

  // Must compute:
// uo = ( p L + q M ) u
//                     --              -- --  --
//                     | G11   G12  G13 | | D1 |
//p L u = [D1 D2 D3]^T | G21   G22  G23 | | D2 | u
//                     | G31   G32  G33 | | D3 |
//                     --              -- --  --
// where
//    D1u = (I_X_I_X_Dx)u; D2u = (I_X_Dy_X_I)u; D3u = (Dz_X_I_X_I)u
// and Gij are the metric terms computed in init method, that include
// the weights and Jacobian; p is 'viscosity', M is mass matrix, 
// and q a factor for M.

// Note: for 2d embedded elements, the basis functions are just
//       Phi_J = h_i(xi) h_j(eta) zeta
//       where zeta is the embedding coordinate, so the Dz is
//       just the identity matrix, as is D3; also the # tensor
//       products are reduced to 2; e.g.:
//    D1u = (I_X_Dx)u; D2u = (Dy_X_I)u; D3u = (I_X_I)u
// p is 'viscosity', M is mass matrix, and q a factor for M

  // Re-arrange local temp space for divergence:
  for ( GSIZET i=0; i<GDIM; i++ ) gdu[i] = utmp[i+GDIM];

  // Compute derivatives of u:
  GMTK::compute_grefderivs(*grid_, u, etmp1_, FALSE, utmp); // utmp stores tensor-prod derivatives, Dj u
  
  // Compute Gij (D^j u): 
  for ( GSIZET i=0; i<GDIM; i++ ) { 
    *utmp[GDIM+i] = 0.0;
    for ( GSIZET j=0; j<GDIM; j++ ) {
      utmp[j]->pointProd(*G_(i,j), uo); // Gij * du^j
      *utmp[GDIM+i] += uo;
    }
  }

  // Apply p parameter ('viscosity') if necessary to Laplacian:
  if ( p_ != NULLPTR ) {
    if ( p_->size() >= grid_->ndof() ) uo.pointProd(*p_);
    else if ( (*p_)[0] != 1.0 )  uo *= (*p_)[0];
  }

  // {utmp[GDIM+1], utmp[GDIM+2], uo} now hold Ti = p Gij D^j u, i = 0, GDIM-1
  // Now compute DT^j ( t^j ):
  GMTK::compute_grefdiv(*grid_, gdu, etmp1_, TRUE, uo); // Compute 'divergence' with DT_j


  // At this point, we have uo = L u
 

  // Apply Mass operator and q parameter if necessary:
  if ( bcompute_helm_ ) {
    massop_->opVec_prod(u, utmp, *utmp[0]); // final arg unused
    if ( q_ != NULLPTR ) {
      if ( q_->size() >= grid_->ndof() ) utmp[0]->pointProd(*q_);
      else if ( (*q_)[0] != 1.0)        *utmp[0] *= (*q_)[0];
    }
    uo += *utmp[0];
  }

} // end of method embed_prod


//**********************************************************************************
//**********************************************************************************
// METHOD : reg_prod
// DESC   : Compute application of this operator to input vector for 
//          GE_REGULAR elements
//          NOTE: must have at least GDIM utmp vectors set via set_tmp method
// ARGS   : u  : input vector
//          utmp: tmp space
//          uo : output (result) vector
//             
// RETURNS:  none
//**********************************************************************************
void GHelmholtz::reg_prod(GTVector<GFTYPE> &u, 
                          GTVector<GTVector<GFTYPE>*> &utmp,
                          GTVector<GFTYPE> &uo)
{

  assert( utmp.size() >= GDIM
       && "Insufficient temp space specified");

  GSIZET            ibeg, iend;
  GTVector<GTVector<GFTYPE>*> gdu(GDIM);
  GElemList        *gelems=&grid_->elems();
  GElem_base       *elem;

  // Compute:
  //   uo = ( p L + q M ) u
  // where
  // p L u = W ( W^-1 D1T G11 W  D1 + W^-1 D2T G22 W D2 + W^-1 D3T G33 W D3 )u,
  // and
  //   D1 = I_X_I_X_Dx u ; D2 = I_X_Dy_X_I u; D2 = Dz_X_I_X_I u
  // and M is just the mass matrix and p, and q are Laplacian and Mass factors that
  // may be position dependent. Here, W are the 1d mass weights, and M is the GDIM-d
  // tensor product form.

  // Re-arrange local temp space for divergence:
  for ( GSIZET i=0; i<GDIM; i++ ) gdu[i] = utmp[i];

  // Compute weighted deriviatives of u:
#if 1
  GMTK::compute_grefderivsW(*grid_, u, etmp1_, FALSE, gdu); // utmp stores tensor-prod derivatives
#else
  GMTK::compute_grefderivs(*grid_, u, etmp1_, FALSE, gdu); // utmp stores tensor-prod derivatives
#endif

  // Multiply by (element-size const) metric factors, possibly x-dependent 
  // 'viscosity':
  for ( GSIZET k=0; k<GDIM; k++ ) {
    gdu[k]->pointProd(*G_(k,0));
    // Apply p parameter ('viscosity') if necessary to Laplacian:
    if ( p_ != NULLPTR ) {
      if ( p_->size() >= grid_->ndof() ) gdu[k]->pointProd(*p_);
      else if ( (*p_)[0] != 1.0 )  *gdu[k] *= (*p_)[0];
    }
  }

  // Take 'divergence' with D:
  GMTK::compute_grefdivW(*grid_, gdu, etmp1_, FALSE, uo); // Compute 'divergence' with W^-1 D_j

  // Apply mass operator (includes Jacobian already):
  *gdu[0] = uo;
  massop_->opVec_prod(*gdu[0], gdu, uo); // tmp array does nothing

  // Add q X mass operator (this really defines the Helmholtz op)  if necessary:
  if ( bcompute_helm_ ) {
    if ( q_ != NULLPTR ) {
      massop_->opVec_prod(u, gdu, *utmp[0]);
      if ( q_->size() >= grid_->ndof() ) *utmp[0].pointProd(*q_);
      else if ( (*q_)[0] != 1.0)         *utmp[0] *= (*q_)[0];
    }
    uo += *utmp[0];
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
void GHelmholtz::init()
{
  assert(grid_->ntype().multiplicity(0) == GE_MAX-1 
        && "Only a single element type allowed on grid");

  if ( grid_->itype(GE_2DEMBEDDED).size() > 0 
    || grid_->itype  (GE_DEFORMED).size() > 0 ) {
    def_init();
  }
  if ( grid_->itype(GE_REGULAR).size() > 0 ) {
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
void GHelmholtz::def_init()
{

  if ( grid_->itype(GE_2DEMBEDDED).size() == 0 
    && grid_->itype  (GE_DEFORMED).size() == 0 ) return;

  if ( grid_->gtype() != GE_2DEMBEDDED
    && grid_->gtype() != GE_DEFORMED ) return; 

  GTVector<GSIZET>             N(GDIM);
  GTMatrix<GTVector<GFTYPE>>  *dXidX;    // element-based dXi/dX matrix
  GTVector<GTVector<GFTYPE>*>  W(GDIM);  // element-based weights
  GTVector<GFTYPE>            *Jac;      // element-based Jacobian
  GElemList                   *gelems = &grid_->elems();

  // Compute metric components:
  // G = Sum_k dxi^i/dx^k dxi^j/dx^k * Jac * W.
  //
  // Metric is symmetric, so we use pointers to 
  // ensure repeated elements aren't duplicated
  // in memory. Each element points to a vector:
  GSIZET nxy = grid_->gtype() == GE_2DEMBEDDED ? GDIM+1: GDIM;
  GSIZET ibeg, iend; // beg, end indices for global array
  G_ .resize(nxy,nxy);
  G_ = NULLPTR;
  for ( GSIZET j=0; j<nxy; j++ ) {
    for ( GSIZET i=j; i<nxy; i++ ) {
      G_ (i,j) = new GTVector<GFTYPE>(grid_->ndof());
    }
    for ( GSIZET i=0; i<j; i++ ) {
      G_ (i,j) = G_(j,i);
    }
  }


  // Cycle through all elements; fill metric elements


  Jac = &grid_->Jac();
  dXidX = &grid_->dXidX();


  for ( GSIZET e=0; e<grid_->elems().size(); e++ ) {
    if ( (*gelems)[e]->elemtype() != GE_DEFORMED 
      && (*gelems)[e]->elemtype() != GE_2DEMBEDDED ) continue;

    // Restrict global data to element range:
    ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
    for ( GSIZET j=0; j<GDIM; j++ ) {
      W[j]= (*gelems)[e]->gbasis(j)->getWeights();
      N[j]= (*gelems)[e]->size(j);
    }
    Jac->range(ibeg, iend);
    for ( GSIZET j=0; j<nxy; j++ ) 
      for ( GSIZET i=0; i<nxy; i++ ) (*dXidX)(i,j).range(ibeg, iend);
    

#if defined(_G_IS2D)

    for ( GSIZET j=0; j<nxy; j++ ) { // G matrix element col
      for ( GSIZET i=j; i<nxy; i++ ) { // G matrix element row
        (*G_(i,j)).range(ibeg, iend); // restrict global vec to local range
        for ( GSIZET k=0; k<nxy; k++ ) {
          for ( GSIZET m=0, n=0; m<N[1]; m++ ) {
            for ( GSIZET l=0; l<N[0]; l++,n++ ) {
              (*G_(i,j))[n] = (*dXidX)(i,k)[n] * (*dXidX)(j,k)[n]
                            * (*W[0])[l] * (*W[1])[m] * (*Jac)[n];
            }
          }
        }
        (*G_(i,j)).range_reset(); // reset to global range
      }
    }

#else

    for ( GSIZET j=0; j<nxy; j++ ) { // G matrix element col
      for ( GSIZET i=0; i<nxy; i++ ) { // G matrix element row
        (*G_(i,j)).range(ibeg, iend); // restrict global vec to local range
        for ( GSIZET k=0; k<nxy; k++ ) {
          for ( GSIZET p=0, n=0; p<N[2]; p++ ) {
            for ( GSIZET m=0; m<N[1]; m++ ) {
              for ( GSIZET l=0; l<N[0]; l++,n++ ) {
                (*G_(i,j))[n] = (*dXidX)(i,k)[n] * (*dXidX)(j,k)[n]
                              * (*W[0])[l] * (*W[1])[m] * (*W[2])[p] * (*Jac)[n];
              }
            }
          }
        }
        (*G_(i,j)).range_reset(); // reset global vec to global range
      }
    }

#endif
  }

  // Reset ranges to global scope:
  Jac->range_reset();
  for ( GSIZET j=0; j<dXidX->size(2); j++ ) 
    for ( GSIZET i=0; i<dXidX->size(1); i++ ) (*dXidX)(i,j).range_reset();

} // end of method def_init


//**********************************************************************************
//**********************************************************************************
// METHOD : reg_init
// DESC   : Compute metric components for regular elements. Mass matrix may be 
//          instantiated here if it hasn't already been set.
// ARGS   : none
// RETURNS: none
//**********************************************************************************
void GHelmholtz::reg_init()
{
  if ( grid_->gtype() != GE_REGULAR ) return; 

  GTVector<GSIZET>             N(GDIM);
  GTVector<GTVector<GFTYPE>*>  W(GDIM);  // element-based weights
  GTMatrix<GTVector<GFTYPE>>  *dXidX;    // element-based dXi/dX matrix
  GElemList                   *gelems = &grid_->elems();

  // Compute metric components:
  // G = Sum_k dxi^i/dx^k dxi^j/dx^k 
  // Note: weights & Jac not included here, so mass matrix
  //       must be set prior to entry, or created here. Mass 
  //       is multiplied in the reg_prod method
  // For regular elements, Gij = 0 for i!= j, and
  // Are constant for i=j. So we make the valid
  // metric components constant vectors, and
  // set those.
  //
  // For reg elems, we have
  // 
  //    Jac = L1 L2 L3/8 (3d) or L1 L2 / 4 (2d)
  //  G11 = (dXi_1/dx)^2 = 4 / L1^2 (2d, 3d)
  //  G22 = (dXi_2/dy)^2 = 4 / L2^2 (2d, 3d)
  //  G33 = (dXi_3/dz)^2 = 4 / L3^2 (2d, 3d)
  // 
  // Metric elements dXi_j/dX*i should already have
  // been computed with the grid.
  
  dXidX = &grid_->dXidX();

  GSIZET nxy = grid_->gtype() == GE_2DEMBEDDED ? GDIM+1: GDIM;
  GSIZET ibeg, iend; // beg, end indices for global array
  G_ .resize(nxy,1);
  G_ = NULLPTR;
  for ( GSIZET j=0; j<nxy; j++ ) {
    G_ (j,0) = new GTVector<GFTYPE>(grid_->ndof());
  }

  // Don't include Jacobian in this quantity for reg elements:
  for ( GSIZET j=0; j<GDIM; j++ ) {
    for ( GSIZET i=0; i<(*G_(j,0)).size(); i++ ) {  
      (*G_(j,0))[i] = (*dXidX)(j,0)[i] * (*dXidX)(j,0)[i];
    }
  }

  massop_ = new GMass(*grid_);
  bown_mass_ = TRUE;

} // end of method reg_init


//**********************************************************************************
//**********************************************************************************
// METHOD : set_Lap_scalar
// DESC   : Set scalar multiplying Laplacian operator
// ARGS   : 
//          p     : 'viscosity' parameter, global
//             
// RETURNS:  none
//**********************************************************************************
void GHelmholtz::set_Lap_scalar(GTVector<GFTYPE> &p)
{
  assert(p.size() == 1 || p.size() >= grid_->ndof() 
       && "Viscosity parameter of insufficient size");
  
  if ( p_ != NULLPTR && bown_p_ ) delete p_;

  p_ = &p;
  bown_p_ = FALSE;

} // end of method set_Lap_scalar


//**********************************************************************************
//**********************************************************************************
// METHOD : set_mass_scalar
// DESC   : Set scalar multiplying mass operator
// ARGS   : 
//          q     : Mass parameter, global
//             
// RETURNS:  none
//**********************************************************************************
void GHelmholtz::set_mass_scalar(GTVector<GFTYPE> &q)
{
  assert(q.size() >= grid_->ndof() 
       && "Mass parameter of insufficient size");
  
  if ( q_ != NULLPTR && bown_q_ ) delete q_;

  q_ = &q;
  bown_q_ = FALSE;
  bcompute_helm_ = TRUE;

} // end of method set_mass_scalar


//**********************************************************************************
//**********************************************************************************
// METHOD : set_mass
// DESC   : Set mass operator. This method tell operator to compute full Helmholtz
//          operator, and not just Laplacian
// ARGS   : 
//          mass: global mass operator
//             
// RETURNS:  none
//**********************************************************************************
void GHelmholtz::set_mass(GMass &m)
{

  if ( massop_ != NULL  && bown_mass_ ) delete massop_;
  
  if ( bown_mass_ && massop_ != NULLPTR ) delete massop_;
  massop_ = &m;
  bown_mass_ = FALSE;

} // end of method set_mass


