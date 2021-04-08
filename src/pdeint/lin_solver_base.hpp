/*
 * linsolver_base.hpp
 *
 *  Created on: Mar 7, 2020 
 *      Author: d.rosenberg
 */

#ifndef SRC_PDEINT_LINSOLVERBASE_HPP_
#define SRC_PDEINT_LINSOLVERBASE_HPP_


#include <memory>
#include <vector>

namespace geoflow {
namespace pdeint {

/**
 * IO base to carry out state IO 
 *
 */
template<typename TypePack>
class LinSolverBase {

public:
        enum GNormType        {GCG_NORM_INF=0, GCG_NORM_EUC , GCG_NORM_L2 , GCG_NORM_L1, GCG_NORM_NONE};
        // Note: if GNormType is modified, must also modify
        //       static str2normtype method!
        using Types          = TypePack;
	using Operator       = typename Types::Operator;
	using Preconditioner = typename Types::Preconditioner;
	using State          = typename Types::State;
	using StateComp      = typename Types::StateComp;
	using Grid           = typename Types::Grid;
	using Ftype          = typename Types::Ftype;
	using ConnectivityOp = typename Types::ConnectivityOp;

        struct Traits {
          bool         bbdycond = false;    // is there a bdy condition? 
          int          maxit    = 512;      // max no. iterations
          GNormType    normtype = GCG_NORM_INF;
                                            // norm type
          double       tol      = 1e-6;     // tolerance
        };

      
	LinSolverBase() = delete;

	/**
	 * Constructor to initialize everything needed to do linear solve
	 *
	 */
	LinSolverBase(Traits& traits, Grid& grid, ConnectivityOp& ggfx, State& tmppack):
          traits_ (traits), grid_(&grid), ggfx_(&ggfx), tmp_(&tmppack) {}
	LinSolverBase(const LinSolverBase& LS) = default;
	~LinSolverBase() = default;
	LinSolverBase& operator=(const LinSolverBase& LS) = default;

	/**
	 * Do solve LinSolveClass x = b, for x; 
         * akin to matrix inversion:
	 *
	 * @param[in]     b  RHS
	 * @param[in,out] x  (homogeneous) solution
	 */
	int  solve( const StateComp& b,
                    StateComp&       x){
               return this->solve_impl(b, x);
             }
	/**
	 * Do solve, with arbitrary operator,
         * without boundary solution:
	 *
	 * @param[in]     A  Linear operator
	 * @param[in]     b  RHS
	 * @param[in,out] x  solution
	 */
	int  solve( Operator&        A,
                    const StateComp& b,
                    StateComp&       x){
               return this->solve_impl(A, b, x);
             }
	/**
	 * Do solve, with arbitrary operator,
         * with boundary solution
	 *
	 * @param[in]     A  Linear operator
	 * @param[in]     b  RHS
	 * @param[in,out] xb boundary solution
	 * @param[in,out] x  (homogeneous) solution
	 */
	int  solve( Operator&        A,
                    const StateComp& b,
                    const StateComp& xb,
                    StateComp&       x){
               return this->solve_impl(A, b, xb, x);
             }
        /**
	 * Get traits
	 *
	 */
	Traits &get_traits( ) { 
               return traits_;
             }

        /**
	 * str2normtype: get GNormType from string
	 * @param[in] : snorm
	 */
        static GNormType str2normtype(std::string snorm);

protected:
        virtual int solve_impl(const StateComp& b,    
                               StateComp&       x) = 0;

        virtual int solve_impl(Operator&        A,
                               const StateComp& b,    
                               StateComp&       x) = 0;

        virtual int solve_impl(Operator&        A,
                               const StateComp& b,    
                               const StateComp& xb,    
                               StateComp&       x) = 0;

        Traits          traits_;
        Grid           *grid_;
        ConnectivityOp *ggfx_;
        State          *tmp_;
};


} // namespace pdeint
} // namespace geoflow

#include "pdeint/lin_solver_base.ipp"

#endif /* SRC_PDEINT_LINSOLVERBASE_HPP_ */
