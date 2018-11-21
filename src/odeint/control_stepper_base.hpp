/*
 * control_stepper_base.hpp
 *
 *  Created on: Nov 18, 2018
 *      Author: bflynt
 */

#ifndef SRC_ODEINT_CONTROL_STEPPER_BASE_HPP_
#define SRC_ODEINT_CONTROL_STEPPER_BASE_HPP_



#include <cstddef>
#include <memory>

namespace geoflow {
namespace odeint {

/**
 *
 */
template<typename StepperType>
class ControlStepperBase {

public:
	using Equation    = typename StepperType::Equation;
	using State       = typename Equation::State;
	using Value       = typename Equation::Value;
	using Derivative  = typename Equation::Derivative;
	using Time        = typename Equation::Time;
	using Jacobian    = typename Equation::Jacobian;
	using EquationPtr = std::shared_ptr<Equation>;
	using Size        = std::size_t;

	ControlStepperBase() = default;
	ControlStepperBase(const ControlStepperBase& si) = default;
	virtual ~ControlStepperBase() = default;
	ControlStepperBase& operator=(const ControlStepperBase& si) = default;

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
	bool step(State& u, const Time& t, Time& dt){
		return this->step_impl(u,t,dt);
	}

	bool step(State& uin, const Time& t, State& uout, Time& dt){
		return this->step_impl(uin,t,uout,dt);
	}

	bool step(State& u, const Derivative& dudt, const Time& t, Time& dt){
		return this->step_impl(u,dudt,t,t,dt);
	}

	bool step(State& uin, const Derivative& dudt, const Time& t, State& uout, Time& dt){
		return this->step_impl(uin,dudt,t,uout,dt);
	}



protected:

	virtual Size order_impl() const = 0;


	virtual bool step_impl(State& u, const Time& t, Time& dt) = 0;


	virtual bool step_impl(State& uin, const Time& t, State& uout, Time& dt) = 0;


	virtual bool step_impl(State& u, const Derivative& dudt, const Time& t, Time& dt) = 0;


	virtual bool step_impl(State& uin, const Derivative& dudt, const Time& t, State& uout, Time& dt) = 0;

};

} // namespace odeint
} // namespace geoflow


#endif /* SRC_ODEINT_CONTROL_STEPPER_BASE_HPP_ */
