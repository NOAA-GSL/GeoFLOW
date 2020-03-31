//==================================================================================
// Module       : gpoisson.hpp
// Date         : 3/27/20 (DLR)
// Description  : Encapsulates the methods and data associated with
//                a Poisson solver.The (continuous) Poisson equation 
//                is solved:
//                        Nabla^2 (u + ub) = f,
//                where ub is the continuous boundary solution. If
//                disc_rhs_==TRUE, then solver will 'discretize'
//                the RHS, assuming the one provided is smooth/continuous.
//                'Discretizing' RHS means multiplying by -M_L, where
//                M_L is the local mass matrix. Otherwise, it is assumed
//                this is done by caller.
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#include <cmath>
#include <type_traits>
#include <cassert>
#include <limits>
#include "ggfx.hpp"
#include "gcomm.hpp"
#include "tbox/error_handler.hpp"

using namespace std;

//************************************************************************************
//************************************************************************************
// METHOD : Constructor
// DESC   : 
// ARGS   : traits : GCG traits
//          grid   : Grid object
//          Lap    : Laplacian operator
//          precond: Preconditioner
//          ggfx   : Connectivity object
//          tmppack: temp vector list
// RETURNS: GPoisson
//************************************************************************************
template<typename Types>
GPoisson<Types>::GPoisson(Traits& traits, Grid& grid, LapOperator& Lap, Preconditioner* precond,  ConnectivityOp& ggfx, State& tmppack)
comm_            (ggfx.getComm()),
disc_rhs_                  (TRUE),
tmppack_               (&tmppack),
grid_                    (&grid_),
Lap_                       (&Lap),
precond_                (precond),
cg_                     (NULLPTR)
{
  // Save one element of tmppack for computation:
  tmpcg_.resize(tmppack_->size());
  for ( auto j=0; j<tmpcg_.size(); j++ ) tmpcg_[j] = tmppack_[j+1];

  cg_ = new CG<Types>(traits, *grid_, *ggfx_, *tmpcg_);
  if ( precond_ != NULLPTR ) cg_->set_precond(*precond_);

  gmass_ = new GMass(grid);

} // end of constructor (2) method


//************************************************************************************
//************************************************************************************
// METHOD : Destructor
// DESC   : 
// ARGS   : none.
// RETURNS: none.
//************************************************************************************
template<typename Types>
GPoisson<Types>::~GPoisson()
{
  if ( cg_    != NULLPTR ) delete cg_;
  if ( gmass_ != NULLPTR ) delete gmass_;
}


//************************************************************************************
//************************************************************************************
// METHOD : Copy constructor
// DESC   : 
// ARGS   : GPoisson
// RETURNS: none
//************************************************************************************
// Copy constructor method
template<typename Types>
GPoisson<Types>::GPoisson(const GPoisson<Types> &a)
{

} // end of copy constructor method


//************************************************************************************
//************************************************************************************
// METHOD : Assignment operatior
// DESC   : 
// ARGS   : GPoisson
// RETURNS: GPoisson
//************************************************************************************
template<typename Types>
GPoisson<Types>& GPoisson<Types>::operator=(const GPoisson<Types> &a)
{

  return *this;
 
} // end of = operator


//************************************************************************************
//************************************************************************************
// METHOD : solve (1)
// DESC   : Solve implementation to find homogeneous solution. 
//           
// ARGS   :
//          b    : right-hand side vector
//          x    : solution, returned
// RETURNS: integer error code; 0 on success
//************************************************************************************
template<typename Types>
GINT GPoisson<Types>::solve(const StateComp& b, StateComp& x)
{
  GINT iret; 

  // 'Discretize' rhs by multiplying by
  //   b <- -M_L b, wherer M_L is local (not assembled)
  // mass matrix:
  (*tmppack_)[0] = b;
  if ( disc_rhs_ ) { 
    (*tmppack_)[0].pointProd(*(gmass_->data()));
    (*tmppack_)[0] *= -1.0;
  }

  iret = cg.solve(*L_, (*tmppack_)[0], x);

  return iret;

} // end of method solve (1)


//************************************************************************************
//************************************************************************************
// METHOD : solve (2)
// DESC   : Solve implementation, with boundary solution, 
//          assuming an initial guess, x
// ARGS   : 
//          b    : right-hand side vector
//          xb   : boundary solution. Must be smooth, but not assembled
//          x    : solution, x0+xb, returned, where x0=homogeneous solution
// ARGS   : 
// RETURNS: integer error code; 0 on success
//************************************************************************************
template<typename Types>
GINT GPoisson<Types>::solve(const StateComp& b, const StateComp& xb, StateComp& x)
{
  GINT iret;

  // 'Discretize' rhs by multiplying by
  //   b <- -M_L b, wherer M_L is local (not assembled)
  // mass matrix:
  (*tmppack_)[0] = b;
  if ( disc_rhs_ ) { 
    (*tmppack_)[0].pointProd(*(gmass_->data()));
    (*tmppack_)[0] *= -1.0;
  }

  iret = cg.solve(*L_, (*tmppack_)[0], xb, x);
  
  return iret;

} // end of method solve (2)

