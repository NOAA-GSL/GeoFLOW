//==================================================================================
// Module       : gcg.hpp
// Date         : 3/7/20 (DLR)
// Description  : Encapsulates the methods and data associated with
//                a Poisson solver. The (continuous) Poisson equation 
//                is solved:
//                        Nabla^2 (u + ub) = f,
//                where ub is the continuous boundary solution.
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================

#if !defined(_GCG_HPP)
#define _GCG_HPP

#include <sstream>
#include "gtvector.hpp"
#include "gtmatrix.hpp"
#include "glinop.hpp"
#include "gcomm.hpp"
#include "ggfx.hpp"
#include "pdeint/lin_solver_base.hpp"

using namespace geoflow::pdeint;
using namespace std;

template<typename TypePack> 
class GCG : public LinSolverBase<TypePack>
{
public:
                      enum  GCGERR         {GCGERR_NONE=0, GCGERR_NOCONVERGE, GCGERR_SOLVE, GCGERR_PRECOND,  GCGERR_BADITERNUM};
                      using Types          = TypePack;
                      using SolverBase     = LinSolverBase<Types>;
                      using Operator       = typename Types::Operator;
                      using Preconditioner = typename Types::Preconditioner;
                      using State          = typename Types::State;
                      using StateComp      = typename Types::StateComp;
                      using Grid           = typename Types::Grid;
                      using Value          = typename Types::Value;
                      using ConnectivityOp = typename Types::ConnectivityOp;
                      using Traits         = typename SolverBase::Traits;

/*
                      static_assert(std::is_same<Operator,GLinOp>::value,
                                    "Operator is of incorrect type");
*/

                      static_assert(std::is_same<State,GTVector<GTVector<GFTYPE>*>>::value,
                                    "State is of incorrect type");
                      static_assert(std::is_same<StateComp,GTVector<GFTYPE>>::value,
                                    "StateComp is of incorrect type");
                      static_assert(std::is_same<Value,GFTYPE>::value,
                                    "Value is of incorrect type");
                      static_assert(std::is_same<ConnectivityOp,GGFX<Value>>::value,
                                    "ConnectivityOp is of incorrect type");


                      GCG() = delete;
                      GCG(Traits& traits, Grid& grid, ConnectivityOp& ggfx, State& tmppack);
                     ~GCG();
                      GCG(const GCG &a);
                      GCG  &operator=(const GCG &);
                   

                      GINT         solve_impl(Operator& A, const StateComp& b, StateComp& x);
                      GINT         solve_impl(Operator& A, const StateComp& b, 
                                              const StateComp& xb, StateComp& x);
                      GINT         solve_impl(const StateComp& b, StateComp& x) {assert(FALSE);}
                      void         set_precond(Preconditioner& precond)
                                   {precond_ = &precond;}
                      StateComp&   get_residuals() { return residuals_; }  
                      GFTYPE       get_resid_max() { return residmax_; }
                      GFTYPE       get_resid_min() { return residmin_; }
                      GINT         get_iteration_count() { return iter_+1; }  


private:
// Private methods:
                       void        init();
                       GFTYPE      compute_norm(const StateComp& x, State& tmp);
// Private data:
     GC_COMM           comm_;        // communicator
     GBOOL             bInit_;       // initialization flag
     GBOOL             bbv_;         // is there a bdy vector?
     GINT              irank_;       // rank
     GINT              iter_;        // iteration number
     GINT              nprocs_;      // no. tasks
     GFTYPE            residmax_;    // max residual
     GFTYPE            residmin_;    // min residual
     StateComp         residuals_;   // list of resituals for each iteration
     LinSolverBase<TypePack>
                      *precond_;     // preconditioner


};

#include "gcg.ipp"

#endif

