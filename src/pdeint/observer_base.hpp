/*
 * observer_base.hpp
 *
 *  Created on: Nov 27, 2018
 *      Author: bflynt
 */

#ifndef SRC_PDEINT_OBSERVER_BASE_HPP_
#define SRC_PDEINT_OBSERVER_BASE_HPP_

#include <memory>
#include <vector>
#include <string>

namespace geoflow {
namespace pdeint {

/**
 * Base class to provided an interface for all observations.
 *
 * The base class provides the interface using the strategy pattern
 * for all observation types.  Observers take the solution after
 * every time step and can extract any user requested information,
 * display, save states, etc.
 */
template<typename EquationType>
class ObserverBase {

public:
        enum ObsType     {OBS_CYCLE=0, OBS_TIME};
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
                ObsType   itype      = OBS_CYCLE; // obs output type
                std::vector<int>     
                          state_index;            // which state members to observe
                std::vector<std::string>
                          state_names;            // file/ref names for each state member
                size_t    cycle_interval = 10;    // cycle interval for observation
                Time      time_interval = 1.0;    // time interval for observation
        };

	ObserverBase() = delete;
	ObserverBase(Traits& traits, Grid& grid){traits_=traits; grid_= &grid;};
	ObserverBase(const ObserverBase& obs) = default;
	virtual ~ObserverBase() = default;
	ObserverBase& operator=(const ObserverBase& obs) = default;

	/**
	 * Observe the current state at time
	 *
	 * @param[in] Traits structure
	 * @param[in] t Time of current state
	 * @param[in] u State at the current time
	 */
	void observe(const Time t, const State& u){
		this->observe_impl(t,u);
	}

protected:
        Traits traits_;
        Grid  *grid_;
	/**
	 * Must be provided by implementation
	 */
	virtual void observe_impl(const Traits& traits, const Time& t, const State& u) = 0;

};

} // namespace pdeint
} // namespace geoflow

#endif /* SRC_PDEINT_OBSERVER_BASE_HPP_ */
