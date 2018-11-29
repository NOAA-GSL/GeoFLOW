/*
 * stepper_base.hpp
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
	 * Get the system of equations we are stepping
	 */
	const EquationPtr& getEquationPtr() const{
		return this->eqn_ptr_;
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
	 * \param[out] uerr Error of current step
	 */
	void step(const Time& t, const Time& dt, State& u){
		this->step_impl(t,dt,u);
	}

	void step(const Time& t, const State& uin, const Time& dt, State& uout){
		this->step_impl(t,uin,dt,uout);
	}

	void step(const Time& t,  const Time& dt, const Derivative& dudt, State& u){
		this->step_impl(t,dt,dudt,u);
	}

	void step(const Time& t, const State& uin, const Derivative& dudt, const Time& dt, State& uout){
		this->step_impl(t,uin,dudt,dt,uout);
	}



protected:

	EquationPtr eqn_ptr_;


	virtual Size order_impl() const = 0;


	virtual void step_impl(const Time& t, const Time& dt, State& u) = 0;


	virtual void step_impl(const Time& t, const State& uin, const Time& dt, State& uout) = 0;


	virtual void step_impl(const Time& t, const Time& dt, const Derivative& dudt, State& u) = 0;


	virtual void step_impl(const Time& t, const State& uin, const Derivative& dudt, const Time& dt, State& uout) = 0;

};

} // namespace odeint
} // namespace geoflow

#endif /* SRC_ODEINT_STEPPER_BASE_HPP_ */
