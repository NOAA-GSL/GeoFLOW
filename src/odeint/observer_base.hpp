/*
 * observer_base.hpp
 *
 *  Created on: Nov 27, 2018
 *      Author: bflynt
 */

#ifndef SRC_ODEINT_OBSERVER_BASE_HPP_
#define SRC_ODEINT_OBSERVER_BASE_HPP_

#include <memory>

namespace geoflow {
namespace odeint {

template<typename EquationType>
class ObserverBase {

public:
	using Equation    = EquationType;
	using State       = typename Equation::State;
	//using Value       = typename Equation::Value;
	//using Derivative  = typename Equation::Derivative;
	using Time        = typename Equation::Time;
	//using Jacobian    = typename Equation::Jacobian;
	//using Size        = typename Equation::Size;
	//using EquationPtr = std::shared_ptr<Equation>;


	void observe( const Time& t, const State& u){
		this->observe_impl(t,u);
	}

protected:

	virtual void observe_impl(const Time& t, const State& u) = 0;

};


} // namespace odeint
} // namespace geoflow


#endif /* SRC_ODEINT_OBSERVER_BASE_HPP_ */
