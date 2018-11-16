//==================================================================================
// Module       : gmass.hpp
// Date         : 10/19/18 (DLR)
// Description  : Represents the SEM mass operator. Mass-lumping is assumed.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : GLinOp
//==================================================================================

#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include "gmass.hpp"
#include "gtmatrix.hpp"


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Default constructor
// ARGS   : none
// RETURNS: none
//**********************************************************************************
GMass::GMass(GGrid &grid)
: GLinOp(grid),
bmasslumped_      (TRUE)
{
  grid_ = &grid;
  init();
} // end of constructor method (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
GMass::~GMass()
{
} // end, destructor

#if 0
//**********************************************************************************
//**********************************************************************************
// METHOD : do_mass_lumping
// DESC   : Set flag to do mass lumping
// ARGS   : bml : boolean flag (TRUE or FALSE)
// RETURNS:  none
//**********************************************************************************
void GMass::do_mass_lumping(GBOOL bml) 
{
  bmasslumped_ = bml;

} // end of method do_mass_lumping
#endif


//**********************************************************************************
//**********************************************************************************
// METHOD : init
// DESC   : initialize operator
// ARGS   : none.
// RETURNS: none.
//**********************************************************************************
void GMass::init()
{

  #if defined(_G_IS1D)
    init1d();
  #elif defined(_G_IS2D)
    init2d();
  #elif defined(_G_IS3D)
    init3d();
  #endif
  bInitialized_ = TRUE;

} // end of method init


//**********************************************************************************
//**********************************************************************************
// METHOD : init1d
// DESC   : initialize 1d  operator
// ARGS   : none.
// RETURNS: none.
//**********************************************************************************
void GMass::init1d()
{

  // Fill mass vector with weights. This would have to be
  // a matrix if we were not doing mass lumping:
  GSIZET n;
  GTVector<GINT> N(GDIM);
  GTVector<GTVector<GFTYPE>*> W(GDIM);
  GTVector<GFTYPE> *Jac;
  GElemList *gelems;

 
  gelems = &grid_->elems(); 
  mass_.resize(grid_->ndof());
  mass_ = 0.0;
  for ( GSIZET i=0, n=0; i<gelems->size(); i++ ) {
    for ( GSIZET j=0; j<GDIM; j++ ) {
      Jac     = &(*gelems)[i]->Jac();
      W[i]    = (*gelems)[i]->gbasis(j)->getWeights();
      N[j]    = (*gelems)[i]->size(j);
    }
    for ( GSIZET j=0; j<N[0]; j++,n++ ) {
      mass_[n] = (*W[0])[j]*(*Jac)[j];
    }
  }

} // end of method init1d


//**********************************************************************************
//**********************************************************************************
// METHOD : init2d
// DESC   : initialize 2d  operator
// ARGS   : none.
// RETURNS: none.
//**********************************************************************************
void GMass::init2d()
{


  GSIZET n;
  GTVector<GINT> N(GDIM);
  GTVector<GTVector<GFTYPE>*> W(GDIM);
  GTVector<GFTYPE> *Jac;
  GElemList *gelems;


  // Fill Fill mass vector with tensor product of weights and
  // Jacobian: 
  gelems = &grid_->elems();
  mass_.resize(grid_->ndof());
  for ( GSIZET i=0, n=0; i<gelems->size(); i++ ) {
    Jac     = &(*gelems)[i]->Jac();
    for ( GSIZET j=0; j<GDIM; j++ ) {
      W[j]    = (*gelems)[i]->gbasis(j)->getWeights();
      N[j]    = (*gelems)[i]->size(j);
    }
    for ( GSIZET k=0; k<N[1]; k++ ) {
      for ( GSIZET j=0; j<N[0]; j++,n++ ) {
        mass_[n] = (*W[1])[k]*(*W[0])[j] 
                 * (*Jac)[j+k*N[0]];
      }
    }
  }

} // end of method init2d


//**********************************************************************************
//**********************************************************************************
// METHOD : init3d
// DESC   : initialize 3d  operator
// ARGS   : none.
// RETURNS: none.
//**********************************************************************************
void GMass::init3d()
{


  GSIZET n;
  GTVector<GINT> N(GDIM);
  GTVector<GTVector<GFTYPE>*> W(GDIM);
  GTVector<GFTYPE> *Jac;
  GElemList *gelems;


  // Fill mass vector with tensor product of weights and
  // Jacobian: 
  gelems = &grid_->elems();
  mass_.resize(grid_->ndof());
  mass_ = 0.0;
  for ( GSIZET i=0, n=0; i<gelems->size(); i++ ) {
    Jac     = &(*gelems)[i]->Jac();
    for ( GSIZET j=0; j<GDIM; j++ ) {
      W[i]    = (*gelems)[i]->gbasis(j)->getWeights();
      N[j]    = (*gelems)[i]->size(j);
    }
    for ( GSIZET l=0; l<N[2]; l++ ) {
      for ( GSIZET k=0; k<N[1]; k++) {
        for ( GSIZET j=0; j<N[0]; j++,n++ ) {
          mass_[n] = (*W[2])[l]*(*W[1])[k]*(*W[0])[j]
                   * (*Jac)[j+k*N[0]+l*N[0]*N[1]];
        }
      }
    }
  }

} // end of method init3d


//**********************************************************************************
//**********************************************************************************
// METHOD : opVec_prod
// DESC   : Compute application of this operator to input vector.
// ARGS   : input : input vector
//          output: output (result) vector
//             
// RETURNS:  none
//**********************************************************************************
void GMass::opVec_prod(GTVector<GFTYPE> &input, GTVector<GFTYPE> &output) 
{

  mass_.pointProd(input, output);

} // end of method opVec_prod
