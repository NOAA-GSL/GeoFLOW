//==================================================================================
// Module       : gpoisson.hpp
// Date         : 3/27/20 (DLR)
// Description  : Encapsulates the methods and data associated with
//                a Poisson solver.The (continuous) Poisson equation 
//                is solved:
//                        Nabla^2 (u + ub) = f,
//                where ub is the continuous boundary solution.
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
bInit_                    (FALSE),
tmppack_               (&tmppack),
grid_                    (&grid_),
Lap_                       (&Lap),
precond_                (precond),
cg_                     (NULLPTR)
{
  cg_ = new CG<Types>(traits, *grid_, *ggfx_, *tmppack_);
  if ( precond_ != NULLPTR ) cg_->set_precond(*precond_);

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
  if ( cg_ != NULLPTR ) delete cg_;
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

  iret = cg.solve(*L_, b, x);

  return iret;

} // end of method solve (1)


//************************************************************************************
//************************************************************************************
// METHOD : solve (2)
// DESC   : Solve implementation, with boundary solution, 
//          assuming an initial guess, x
// ARGS   : 
//          b    : right-hand side vector
//          xb   : boundary solution
//          x    : solution, x0+xb, returned, where x0=homogeneous solution
// ARGS   : 
// RETURNS: integer error code; 0 on success
//************************************************************************************
template<typename Types>
GINT GPoisson<Types>::solve(const StateComp& b, const StateComp& xb, StateComp& x)
{
  GINT iret;

  iret = cg.solve(*L_, b, xb, x);
  
  return iret;

} // end of method solve (2)

