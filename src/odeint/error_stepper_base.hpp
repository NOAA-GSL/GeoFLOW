/*
 * error_stepper_base.hpp
 *
 *  Created on: Nov 16, 2018
 *      Author: bflynt
 */

#ifndef SRC_ODEINT_ERROR_STEPPER_BASE_HPP_
#define SRC_ODEINT_ERROR_STEPPER_BASE_HPP_

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
template<typename EquationType, typename ErrorType = double>
class ErrorStepperBase {

public:
	using Equation    = EquationType;
	using State       = typename Equation::State;
	using Value       = typename Equation::Value;
	using Derivative  = typename Equation::Derivative;
	using Time        = typename Equation::Time;
	using Jacobian    = typename Equation::Jacobian;
	using EquationPtr = std::shared_ptr<Equation>;
	using Size        = std::size_t;
	using Error       = ErrorType;

	ErrorStepperBase() = default;
	ErrorStepperBase(const ErrorStepperBase& si) = default;
	virtual ~ErrorStepperBase() = default;
	ErrorStepperBase& operator=(const ErrorStepperBase& si) = default;

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
	 * \param[in,out] sys Is the system of equations were are solving
	 * \param[in,out] u Is the state of the system
	 * \param[in] t Current time of state u before taking step
	 * \param[in] dt Size of time step to take
	 * \param[out] uerr Error of current step
	 */
	void step(EquationPtr& sys, State& u, const Time& t, const Time& dt, Error& uerr){
		this->step_impl(sys,u,t,dt,uerr);
	}

	void step(EquationPtr& sys, const State& uin, const Time& t, State& uout, const Time& dt, Error& uerr){
		this->step_impl(sys,uin,t,uout,dt,uerr);
	}

	void step(EquationPtr& sys, State& u, const Derivative& dudt, const Time& t, const Time& dt, Error& uerr){
		this->step_impl(sys,u,dudt,t,t,dt,uerr);
	}

	void step(EquationPtr& sys, const State& uin, const Derivative& dudt, const Time& t, State& uout, const Time& dt, Error& uerr){
		this->step_impl(sys,uin,dudt,t,uout,dt,uerr);
	}



protected:

	virtual Size order_impl() const = 0;


	virtual void step_impl(EquationPtr& sys, State& u, const Time& t, const Time& dt, Error& uerr) = 0;


	virtual void step_impl(EquationPtr& sys, const State& uin, const Time& t, State& uout, const Time& dt, Error& uerr) = 0;


	virtual void step_impl(EquationPtr& sys, State& u, const Derivative& dudt, const Time& t, const Time& dt, Error& uerr) = 0;


	virtual void step_impl(EquationPtr& sys, const State& uin, const Derivative& dudt, const Time& t, State& uout, const Time& dt, Error& uerr) = 0;

};

} // namespace odeint
} // namespace geoflow

#endif /* SRC_ODEINT_ERROR_STEPPER_BASE_HPP_ */
