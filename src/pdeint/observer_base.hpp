/*
 * observer_base.hpp
 *
 *  Created on: Nov 27, 2018
 *      Author: bflynt
 */

#ifndef SRC_PDEINT_OBSERVER_BASE_HPP_
#define SRC_PDEINT_OBSERVER_BASE_HPP_

#include <memory>

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
	using Equation    = EquationType;
	using State       = typename Equation::State;
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
                int       istate;     // which state member to observe (-1 for all)
                char**    state_names;// file/ref names for each state member
                size_t    icycle;     // time cycle, e.g.
        };

	ObserverBase();
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
	void observe(const Traits& traits, const Time& t, const State& u){
		this->observe_impl(t,u);
	}

protected:

	/**
	 * Must be provided by implementation
	 */
	virtual void observe_impl(const Traits& traits, const Time& t, const State& u) = 0;

};

} // namespace pdeint
} // namespace geoflow

#endif /* SRC_PDEINT_OBSERVER_BASE_HPP_ */
