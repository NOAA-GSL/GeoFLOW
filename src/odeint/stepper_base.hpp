/*
 * stepper_interface.hpp
 *
 *  Created on: Nov 16, 2018
 *      Author: bflynt
 */

#ifndef SRC_ODEINT_STEPPER_BASE_HPP_
#define SRC_ODEINT_STEPPER_BASE_HPP_

#include <cstddef>
#include <memory>

namespace geoflow {
namespace odeint {

/**
 * \brief
 * Base class to provided an interface strategy for all
 * time stepping schemes.
 *
 * \details
 * The base class provides the strategy for all time step
 * implementations.
 */
template<typename EquationType>
class StepperBase {

public:
	using Equation    = EquationType;
	using State       = typename Equation::State;
	using Value       = typename Equation::Value;
	using Derivative  = typename Equation::Derivative;
	using Time        = typename Equation::Time;
	using Jacobian    = typename Equation::Jacobian;
	using EquationPtr = std::shared_ptr<Equation>;
	using Size        = std::size_t;


	StepperBase() = default;
	StepperBase(const StepperBase& si) = default;
	virtual ~StepperBase() = default;
	StepperBase& operator=(const StepperBase& si) = default;

	/**
	 * Return the order of temporal accuracy for the method
	 */
	Size order() const{
		return this->order_impl();
	}

	/** Take one step from time t to time t+dt.
	 *
	 * The simplest time step method which takes exactly one
	 * time step from t to t+dt.
	 *
	 * \param[in,out] eqn Is the eqntem of equations were are solving
	 * \param[in,out] u Is the state of the eqntem
	 * \param[in] t Current time of state u before taking step
	 * \param[in] dt Size of time step to take
	 */
	void step(EquationPtr& eqn, State& u, const Time& t, const Time& dt){
		this->step_impl(eqn,u,t,dt);
	}

	void step(EquationPtr& eqn, const State& uin, const Time& t, State& uout, const Time& dt){
		this->step_impl(eqn,uin,t,uout,dt);
	}

	void step(EquationPtr& eqn, State& u, const Derivative& dudt, const Time& t, const Time& dt){
		this->step_impl(eqn,u,dudt,t,t,dt);
	}

	void step(EquationPtr& eqn, const State& uin, const Derivative& dudt, const Time& t, State& uout, const Time& dt){
		this->step_impl(eqn,uin,dudt,t,uout,dt);
	}



protected:

	virtual Size order_impl() const = 0;


	virtual void step_impl(EquationPtr& eqn, State& u, const Time& t, const Time& dt) = 0;


	virtual void step_impl(EquationPtr& eqn, const State& uin, const Time& t, State& uout, const Time& dt) = 0;


	virtual void step_impl(EquationPtr& eqn, State& u, const Derivative& dudt, const Time& t, const Time& dt) = 0;


	virtual void step_impl(EquationPtr& eqn, const State& uin, const Derivative& dudt, const Time& t, State& uout, const Time& dt) = 0;

};

} // namespace odeint
} // namespace geoflow



#endif /* SRC_ODEINT_STEPPER_BASE_HPP_ */
