/*
 * null_observer.hpp
 *
 *  Created on: Nov 27, 2018
 *      Author: bflynt
 */

#ifndef SRC_ODEINT_NULL_OBSERVER_HPP_
#define SRC_ODEINT_NULL_OBSERVER_HPP_

#include <memory>

#include "odeint/observer_base.hpp"


namespace geoflow {
namespace odeint {

template<typename EquationType>
class NullObserver : public ObserverBase<EquationType> {

public:
	using Interface  = ObserverBase<EquationType>;
	using State      = typename Interface::State;
	using Time       = typename Interface::Time;

protected:

	void observe_impl(const Time& /* t */, const State& /* u */){
	}

};


} // namespace odeint
} // namespace geoflow



#endif /* SRC_ODEINT_NULL_OBSERVER_HPP_ */
