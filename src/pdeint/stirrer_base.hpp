/*
 * stirrer_base.hpp
 *
 *  Created on: Mar 27, 2019
 *      Author: bflynt, d.rosemnerg
 */

#ifndef SRC_PDEINT_STIRRER_BASE_HPP_
#define SRC_PDEINT_STIRRER_BASE_HPP_

#include <memory>
#include <vector>
#include <string>

namespace geoflow {
namespace pdeint {

/**
 * Base class to provided an interface for all 'stirrers'
 *
 * The base class provides the interface using the strategy pattern
 * for all stirrer types.  Stirrers take the forcing after
 * every time step and determine if it should 'randomized', and
 * how.
 */
template<typename EquationType>
class StirrerBase {

public:
        enum StirType     {STIR_CYCLE=0, STIR_TIME};
	using Equation    = EquationType;
	using State       = typename Equation::State;
	using Grid        = typename Equation::Grid;
	using Value       = typename Equation::Value;
	using Derivative  = typename Equation::Derivative;
	using Time        = typename Equation::Time;
	using Jacobian    = typename Equation::Jacobian;
	using Size        = typename Equation::Size;
	using EquationPtr = std::shared_ptr<Equation>;

        /**
         * Data structure to hold user selected parameters
         */
        struct Traits {
                StirType  itype          = STIR_TIME; // stir cadence type
                std::vector<int>     
                          state_index;                // which state members to 'stir'
                size_t    corr_cycle     = 10;        // correlation cycle interval 
                double    corr_time      = 1.0;       // correlation time interval 
        };

        StirrerBase() = default;
	StirrerBase(Traits& traits, Grid& grid){traits_=traits; grid_= &grid;}
	StirrerBase(const StirrerBase& obs) = default;
	virtual ~StirrerBase() = default;
	StirrerBase& operator=(const StirrerBase& obs) = default;

	/**
	 * Stir the forcing state at time
	 *
	 * @param[in] t Time of current state
	 * @param[in] u State at the current time
	 * @param[in] uf forcing tendency at the current time
	 */
	void stir(const Time t, const State& u, State& uf){
		this->stir_impl(t,u,uf);
	}

protected:
        Traits traits_;
        Grid  *grid_;
	/**
	 * Must be provided by implementation
	 */
	virtual void stir_impl(const Time& t, const State& u, State &uf) = 0;

};

} // namespace pdeint
} // namespace geoflow

#endif /* SRC_PDEINT_STIRRER_BASE_HPP_ */
