/*
 * null_observer.hpp
 *
 *  Created on: Nov 27, 2018
 *      Author: bflynt
 */

#ifndef SRC_PDEINT_NULL_OBSERVER_HPP_
#define SRC_PDEINT_NULL_OBSERVER_HPP_


#include "pdeint/observer_base.hpp"

namespace geoflow {
namespace pdeint {

/**
 * Do nothing Observer implementation
 *
 * The NullObserver is the default observer when nothing else is
 * provided and does nothing when called.
 *
 * @see ObserverBase
 */
template<typename EquationType>
class NullObserver : public ObserverBase<EquationType> {

public:
	using Interface  = ObserverBase<EquationType>;
	using State      = typename Interface::State;
        using StateInfo  = typename Interface::StateInfo;
	using Grid       = typename Interface::Grid;
	using Time       = typename Interface::Time;
        using Equation   = EquationType;
        using EqnBase    = EquationBase<Equation>;
        using EqnBasePtr = std::shared_ptr<EqnBase>;

        
        NullObserver(EqnBasePtr& equation, Grid& grid){}

protected:

	void observe_impl(const Time& /* t */, const Time& dt /* dt */, const State& /* u */, const State& /* uf */){
		// Do nothing ...
	}
	void init_impl(StateInfo &info) {}

};


} // namespace pdeint
} // namespace geoflow



#endif /* SRC_PDEINT_NULL_OBSERVER_HPP_ */
