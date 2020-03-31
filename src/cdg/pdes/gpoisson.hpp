//==================================================================================
// Module       : gpoisson.hpp
// Date         : 3/27/20 (DLR)
// Description  : Encapsulates the methods and data associated with
//                a Poisson solver. he (continuous) Poisson equation
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

#if !defined(_GPOISSON_HPP)
#define _GPOISSON_HPP

#include <sstream>
#include "gtvector.hpp"
#include "gtmatrix.hpp"
#include "gmass.hpp"
#include "gcg.hpp"
#include "gcomm.hpp"
#include "ggfx.hpp"
#include "pdeint/lin_solver_base.hpp"

using namespace geoflow::pdeint;
using namespace std;

template<typename TypePack> 
class GPoisson 
{
public:
                      using Types          = TypePack;
                      using SolverBase     = LinSolverBase<Types>;
                      using LapOperator    = typename types::LapOperator
                      using Preconditioner = typename Types::Preconditioner;
                      using State          = typename Types::State;
                      using StateComp      = typename Types::StateComp;
                      using Grid           = typename Types::Grid;
                      using Value          = typename Types::Value;
                      using ConnectivityOp = typename Types::ConnectivityOp;
                      using Traits         = typename SolverBase::Traits;

                      static_assert(std::is_same<State,GTVector<GTVector<GFTYPE>*>>::value,
                                    "State is of incorrect type");
                      static_assert(std::is_same<StateComp,GTVector<GFTYPE>>::value,
                                    "StateComp is of incorrect type");
                      static_assert(std::is_same<Value,GFTYPE>::value,
                                    "Value is of incorrect type");
                      static_assert(std::is_same<ConnectivityOp,GGFX<Value>>::value,
                                    "ConnectivityOp is of incorrect type");


                      GPoisson() = delete;
                      GPoisson(Traits& traits, Grid& grid, LapOperator& L, Preconditioner* P, ConnectivityOp& ggfx, State& tmppack);
                     ~GPoisson();
                      GPoisson(const GPoisson &a);
                      GPoisson  &operator=(const GPoisson &);
                   

                      GINT         solve(const StateComp& b, StateComp& u);
                      GINT         solve(const StateComp& b, 
                                     const StateComp& ub, StateComp& u);
                      void         discretize_rhs(GBOOL flag) { discrete_rhs_=flag; };

                      StateComp&   get_residuals() { return cg_->get_residuals(); }  
                      GFTYPE       get_resid_max() { return cg_->get_resid_max(); }
                      GFTYPE       get_resid_min() { return cg_->get_resod_min(); }
                      GINT         get_iteration_count() { return cg_->get_iteration_count(); }  


private:
// Private methods:
                      
// Private data:
     GC_COMM          comm_;        // communicator
     GBOOL            disc_rhs_;    // discretize RHS (multiply by -M_L)
     State           *tmppack_;     // tmp vector pool
     State            tmpcg_;       // tmp vector pool given to linear solver
     Grid            *grid_;        // grid object
     GMass           *gmass_;       // mass matrix operator
     Preconditioner  *precond_;     // preconditioner
     LapOperator     *L_;           // Laplacian operator
     GCG<Types>      *cg_;          // Conjugate gradient operator


};

#include "gpoisson.ipp"

#endif

