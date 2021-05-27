//==================================================================================
// Module       : gcg.hpp
// Date         : 3/7/20 (DLR)
// Description  : Encapsulates the methods and data associated with
//                a conjugate Gradient (CG) solver
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#include <cmath>
#include <type_traits>
#include <cassert>
#include <limits>
#include "ggfx.hpp"
#include "gcomm.hpp"
#include "gmtk.hpp"
#include "tbox/error_handler.hpp"

using namespace std;

//************************************************************************************
//************************************************************************************
// METHOD : Constructor
// DESC   : 
// ARGS   : GD_COMM object, defaults to GD_COMM_WORLD
// RETURNS: GCG
//************************************************************************************
template<typename Types>
GCG<Types>::GCG(Traits& traits, Grid& grid, ConnectivityOp& ggfx, State& tmppack)
: SolverBase(traits, grid, ggfx, tmppack),
comm_            (ggfx.getComm()),
bInit_                    (FALSE),
bbv_                      (FALSE),
nprocs_                       (1),
irank_                        (0),
residmax_                   (0.0),
residmin_
  (std::numeric_limits<GFTYPE>::max()),
precond_            (NULLPTR)
{
  irank_   = GComm::WorldRank(comm_);
  nprocs_ = GComm::WorldSize(comm_);

} // end of constructor (2) method


//************************************************************************************
//************************************************************************************
// METHOD : Destructor
// DESC   : 
// ARGS   : none.
// RETURNS: none.
//************************************************************************************
template<typename Types>
GCG<Types>::~GCG()
{
}


//************************************************************************************
//************************************************************************************
// METHOD : Copy constructor
// DESC   : 
// ARGS   : GCG
// RETURNS: none
//************************************************************************************
// Copy constructor method
template<typename Types>
GCG<Types>::GCG(const GCG<Types> &a)
{

} // end of copy constructor method


//************************************************************************************
//************************************************************************************
// METHOD : Assignment operatior
// DESC   : 
// ARGS   : GCG
// RETURNS: GCG
//************************************************************************************
template<typename Types>
GCG<Types>& GCG<Types>::operator=(const GCG<Types> &a)
{

  return *this;
 
} // end of = operator


//************************************************************************************
//************************************************************************************
// METHOD : init
// DESC   : Performs initialization of global GCG operator
// ARGS   : none.
// RETURNS: none.
//************************************************************************************
template<typename Types>
void GCG<Types>::init()
{
  if ( bInit_ ) return;

  residuals_.resize(this->traits_.maxit);
  bInit_ = TRUE;

} // end of method init


#if 0
//************************************************************************************
//************************************************************************************
// METHOD : solve_impl (1)
// DESC   : Solve implementation to find homogeneous solution. 
//          Taken from van der Vorst text.
//           
// ARGS   : A    : linear operator to invert
//          b    : right-hand side vector
//          x    : solution, returned
// RETURNS: integer error code
//************************************************************************************
template<typename Types>
GINT GCG<Types>::solve_impl(Operator& A, const StateComp& b, StateComp& x)
{
  GINT       iret=GCGERR_NONE;
  GFTYPE     alpha, beta, rnorm, rho, rhom;
  GFTYPE     rtol = this->traits_.tol;
  GFTYPE     eps = 100.0*std::numeric_limits<GFTYPE>::epsilon();
  StateComp *q, *r, *w, *z;
  State      tmp(this->tmp_->size()-4);
  StateComp *mask  = &this->grid_->get_mask();

  StateComp imult_actual(this->grid_->ndof());
  this->ggfx_->get_imult(imult_actual);
  StateComp* imult = &imult_actual;


  assert(this->tmp_->size() > 5);
  init();

  // Set some pointers:
  tmp.resize(this->tmp_->size()-4);
  q  = (*this->tmp_)[0];
  r  = (*this->tmp_)[1];
  w  = (*this->tmp_)[2];
  z  = (*this->tmp_)[3];
  for ( auto j=0; j<this->tmp_->size()-4; j++ ) {
    tmp[j] = (*this->tmp_)[j+4];
  }
  residuals_ = 0.0;

 // Initialize CG loop, and enter loop:
 *r = b;

  A.opVec_prod(x, tmp, *w);             // Ax

 *r -= (*w);                            // r = b - Ax, initial residual

  this->ggfx_->doOp(*r, typename GGFX<decltype(r[0])>::Sum());// DSS r


  // Create effective initial residual
  // and tolerance:
  rnorm = compute_norm(b, tmp);
  rtol = rnorm < 1.0 ? this->traits_.tol : this->traits_.tol * rnorm;

  iter_ = 0; rnorm = 10.0*rtol;

//cout << "solve_impl: rnorm_0=" << rnorm << " traits.tol=" << this->traits_.tol <<  " rtol=" << rtol << endl;

  while ( iret == GCGERR_NONE 
       && iter_ < this->traits_.maxit 
       && rnorm > rtol ) {

    if ( precond_ != NULLPTR ) {        // z = P^-1 r for z,
      iret = precond_->solve(*r, *z);   // where P^-1 is precond
      if ( iret >  0 ) {
        iret = GCGERR_PRECOND; break;
      }
    }
    else {
      *z = *r;                          // use identity preconditioner
    }

    rho  = r->gdot(*z, *imult, comm_);  // rho = r^T imult z
    if ( iter_ == 0 ) {
      *w = *z;
    }
    else {
      beta = rho / rhom;
      GMTK::saxpby<Ftype>(*w, beta, *z, 1.0);  // w = beta w +  z
    }

    A.opVec_prod(*w, tmp, *q);          // q = A w

    this->ggfx_->doOp(*q, typename GGFX<decltype(q[0])>::Sum()); // q <- DSS q

    if ( bbv_ ) q->pointProd(*mask);    // Mask(q)

    alpha = rho /
     (w->gdot(*q, *imult, comm_));      // alpha=rho/w^T imult q

    GMTK::saxpby<Ftype>( x, 1.0, *w, alpha);   // x = x + alpha w
    GMTK::saxpby<Ftype>(*r, 1.0, *q,-alpha);   // r = r - alpha q

    rhom = rho;

    rnorm = compute_norm(*r, tmp);      // residual norm
    residuals_[iter_] = rnorm;
    residmax_ = MAX(rnorm,residmax_);
    residmin_ = MIN(rnorm,residmin_);
    iter_++;

  } // end, CG loop

  if ( bbv_ ) x.pointProd(*mask);

  if ( iret == GCGERR_NONE 
    && iter_ >= this->traits_.maxit 
    && rnorm > this->traits_.tol ) iret = GCGERR_NOCONVERGE;

  return iret;

} // end of method solve_impl (1)

#else

//************************************************************************************
//************************************************************************************
// METHOD : solve_impl (1)
// DESC   : Solve implementation to find homogeneous solution. 
//          Taken from High Order Methods for Incompressible Fluid 
//          Flow, Deville, Fischer, Mund, Cambridge 2002, p 197-198
//           
// ARGS   : A    : linear operator to invert
//          b    : right-hand side vector
//          x    : solution, returned
// RETURNS: integer error code
//************************************************************************************
template<typename Types>
GINT GCG<Types>::solve_impl(Operator& A, const StateComp& b, StateComp& x)
{
  GINT       iret=GCGERR_NONE;
  GFTYPE     alpha, beta, rnorm, rho, rhom;
  GFTYPE     rtol = this->traits_.tol;
  GFTYPE     eps = 100.0*std::numeric_limits<GFTYPE>::epsilon();
  StateComp *q, *r, *w, *z;
  State      tmp(this->tmp_->size()-4);
  StateComp *mask  = &this->grid_->get_mask();

  StateComp imult_actual(this->grid_->ndof());
  this->ggfx_->get_imult(imult_actual);
  StateComp* imult = &imult_actual;

  assert(this->tmp_->size() > 5);
  init();


  // Set some pointers:
  tmp.resize(this->tmp_->size()-4);
  q  = (*this->tmp_)[0];
  r  = (*this->tmp_)[1];
  w  = (*this->tmp_)[2];
  z  = (*this->tmp_)[3];
  for ( auto j=0; j<this->tmp_->size()-4; j++ ) {
    tmp[j] = (*this->tmp_)[j+4];
  }
  residuals_ = 0.0;

 // Initialize CG loop, and enter loop:
 *r = b;

  A.opVec_prod(x, tmp, *w);             // Ax

 *r -= (*w);                            // r = b - Ax, initial residual

  this->ggfx_->doOp(*r, typename GGFX<decltype(r[0])>::Sum());   // DSS r

  if ( bbv_ ) r->pointProd(*mask);      // Mask DSS r
  if ( precond_ != NULLPTR ) {          // solve P z = r for z
    iret = precond_->solve(*r, *z);  
    if ( iret >  0 ) iret = GCGERR_PRECOND; 
  }
  else {
    *z = *r;                            // use identity preconditioner
  }
  *w = *z;                              // w = z, for initial w

  rho = r->gdot(*z, *imult, comm_);     // rho = r^T imult z

  // Create effective initial residual
  // and tolerance:
  rnorm = compute_norm(b, tmp);
  rtol = rnorm < 1.0 ? this->traits_.tol : this->traits_.tol * rnorm;

  iter_ = 0; rnorm = 10.0*rtol;

//cout << "solve_impl: rnorm_0=" << rnorm << " traits.tol=" << this->traits_.tol <<  " rtol=" << rtol << endl;

  while ( iret == GCGERR_NONE 
       && iter_ < this->traits_.maxit 
       && rnorm > rtol ) {

    A.opVec_prod(*w, tmp, *q);          // q = A w
    this->ggfx_->doOp(*q, typename GGFX<decltype(q[0])>::Sum()); // q <- DSS q

    if ( bbv_ ) q->pointProd(*mask);    // Mask(q)

    alpha = rho /
     (w->gdot(*q, *imult, comm_));      // alpha=rho/w^T imult q

    GMTK::saxpby<Ftype>( x, 1.0, *w, alpha);   // x = x + alpha w
    GMTK::saxpby<Ftype>(*r, 1.0, *q,-alpha);   // r = r - alpha q

    if ( precond_ != NULLPTR ) {        // z = P^-1 r for z,
      iret = precond_->solve(*r, *z);   // where P^-1 is precond
      if ( iret >  0 ) {
        iret = GCGERR_PRECOND; break;
      }
    }
    else {
      *z = *r;                          // use identity preconditioner
    }

    rhom = rho;
    rho  = r->gdot(*z, *imult, comm_);  // rho = r^T imult z
    beta = rho / rhom;

    GMTK::saxpby<Ftype>(*w, beta, *z, 1.0);    // w = z + beta w

    rnorm = compute_norm(*r, tmp);      // residual norm
    residuals_[iter_] = rnorm;
    residmax_ = MAX(rnorm,residmax_);
    residmin_ = MIN(rnorm,residmin_);
    iter_++;
//cout << "GCG::solve_impl: iter: ..............." << iter_ << " rnorm=" << rnorm << endl;

  } // end, CG loop

  if ( bbv_ ) x.pointProd(*mask);

cout << "GCG::solve_impl: .....................bv=" << bbv_ << " iterations: " << iter_ << " rnorm=" << rnorm << " traits.tol=" << this->traits_.tol <<  " rtol=" << rtol << " iret=" << iret << endl;

  if ( iret == GCGERR_NONE 
    && iter_ >= this->traits_.maxit 
    && rnorm > this->traits_.tol ) iret = GCGERR_NOCONVERGE;

  return iret;

} // end of method solve_impl (1)

#endif

//************************************************************************************
//************************************************************************************
// METHOD : solve_impl (2)
// DESC   : Solve implementation, with boundary solution, 
//          assuming an initial guess, x
// ARGS   : A    : linear operator to invert
//          b    : right-hand side vector
//          xb   : boundary solution
//          x    : solution, x0+xb, returned, where x0=homogeneous solution
// ARGS   : 
// RETURNS: integer error code
//************************************************************************************
template<typename Types>
GINT GCG<Types>::solve_impl(Operator& A, const StateComp& b, const StateComp& xb, StateComp& x)
{
  GINT iret;

  bbv_ = TRUE; // there is a boundary vector/solution

  // Add in boundary solution, so 
  // that it will form RHS in homogeneous solve:

  if ( bbv_ ) x.pointProd(this->grid_->get_mask()); // apply mask
  x += xb; 

  // Compute homogeneous solution:
  iret = solve_impl(A, b, x);

  // Add back boundary solution:
  x += xb;
 
  return iret;

} // end of method solve_impl (2)


//************************************************************************************
//************************************************************************************
// METHOD : compute_norm
// DESC   : Compute required norm of input vector, based
//          on traits
// ARGS   : x    : vector
//          tmp  : tmp vector
// RETURNS: required norm 
//************************************************************************************
template<typename Types>
GFTYPE GCG<Types>::compute_norm(const StateComp& x, State& tmp)
{
  GFTYPE lret[2], gret[2];

  switch (this->traits_.normtype) {
    case LinSolverBase<Types>::GCG_NORM_INF:
     *tmp[0]  = x; 
//    lret[0] = tmp[0]->infnorm();
      lret[0] = tmp[0]->amax();
      GComm::Allreduce(lret, gret, 1, T2GCDatatype<GFTYPE>() , GC_OP_MAX, comm_);
      break;
    case LinSolverBase<Types>::GCG_NORM_EUC:
     *tmp[0]  = x; tmp[0]->rpow(2);
      lret[0]  = tmp[0]->sum();
      lret[1]  = x.size();
      GComm::Allreduce(lret, gret, 2, T2GCDatatype<GFTYPE>() , GC_OP_SUM, comm_);
      gret[0]  = sqrt(gret[0]/gret[1]);
      break;
    case LinSolverBase<Types>::GCG_NORM_L2:
     *tmp[0]  = x; tmp[0]->rpow(2);
      gret[0]  = this->grid_->integrate(*tmp[0],*tmp[1]);
      gret[0] *= this->grid_->ivolume();
      gret[0]  = sqrt(gret[0]);
      break;
    case LinSolverBase<Types>::GCG_NORM_L1:
     *tmp[0]  = x; tmp[0]->abs();
      gret[0]  = this->grid_->integrate(*tmp[0],*tmp[1]);
      gret[0] *= this->grid_->ivolume();
      break;
    default:
      assert(FALSE);
  }

  return gret[0];

} // end, compute_norm


