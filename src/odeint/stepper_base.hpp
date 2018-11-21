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
	using Size        = typename Equation::Size;
	using EquationPtr = std::shared_ptr<Equation>;


	StepperBase() = default;
	StepperBase(const StepperBase& si) = default;
	virtual ~StepperBase() = default;
	StepperBase& operator=(const StepperBase& si) = default;


	StepperBase(const EquationPtr& eqn) :
		eqn_ptr_(eqn){
	}

	/**
	 * Set the system of equations we are stepping
	 */
	void setEquation(EquationPtr& eqn){
		this->eqn_ptr_ = eqn;
	}

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
	 * \param[in,out] u Is the state of the system of equations
	 * \param[in] t Current time of state u before taking step
	 * \param[in] dt Size of time step to take
	 */
	void step(State& u, const Time& t, const Time& dt){
		this->step_impl(u,t,dt);
	}

	void step(State& uin, const Time& t, State& uout, const Time& dt){
		this->step_impl(uin,t,uout,dt);
	}

	void step(State& u, const Derivative& dudt, const Time& t, const Time& dt){
		this->step_impl(u,dudt,t,t,dt);
	}

	void step(State& uin, const Derivative& dudt, const Time& t, State& uout, const Time& dt){
		this->step_impl(uin,dudt,t,uout,dt);
	}



protected:

	EquationPtr eqn_ptr_;

	virtual Size order_impl() const = 0;


	virtual void step_impl(State& u, const Time& t, const Time& dt) = 0;


	virtual void step_impl(State& uin, const Time& t, State& uout, const Time& dt) = 0;


	virtual void step_impl(State& u, const Derivative& dudt, const Time& t, const Time& dt) = 0;


	virtual void step_impl(State& uin, const Derivative& dudt, const Time& t, State& uout, const Time& dt) = 0;

};

} // namespace odeint
} // namespace geoflow



#endif /* SRC_ODEINT_STEPPER_BASE_HPP_ */
