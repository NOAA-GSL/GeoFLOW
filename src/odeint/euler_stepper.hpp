/*
 * euler_stepper.hpp
 *
 *  Created on: Nov 16, 2018
 *      Author: bflynt
 */

#ifndef SRC_ODEINT_EULER_STEPPER_HPP_
#define SRC_ODEINT_EULER_STEPPER_HPP_


#include "odeint/stepper_base.hpp"


namespace geoflow {
namespace odeint {


template<typename EquationType>
class EulerStepper : public StepperBase<EquationType> {

public:
	// Inheriting these from the Base we insure our
	// interface is consistent and we can perform runtime
	// polymorphism or use them directly.
	using Interface   = StepperBase<EquationType>;
	using Equation    = typename Interface::Equation;
	using State       = typename Interface::State;
	using Value       = typename Interface::Value;
	using Derivative  = typename Interface::Derivative;
	using Time        = typename Interface::Time;
	using Jacobian    = typename Interface::Jacobian;
	using EquationPtr = typename Interface::EquationPtr;
	using Size        = typename Interface::Size;

	EulerStepper() = default;
	EulerStepper(const EulerStepper& eb) = default;
	virtual ~EulerStepper() = default;
	EulerStepper& operator=(const EulerStepper& eb) = default;


	EulerStepper(const EquationPtr& eqn) : Interface(eqn){
	}



protected:


	Size order_impl() const{
		return 1;
	}


	void step_impl(State& u, const Time& t, const Time& dt){
		Derivative dudt(u.size());
		this->eqn_ptr_->dudt(u, dudt, t);
		this->step_impl(u,dudt,t,dt);
	}


	void step_impl(const State& uin, const Time& t, State& uout, const Time& dt){
		uout = uin;
		this->step_impl(uout,t,dt);
	}


	void step_impl(State& u, const Derivative& dudt, const Time& t, const Time& dt){
		u += dt * dudt;
	}


	void step_impl(const State& uin, const Derivative& dudt, const Time& t, State& uout, const Time& dt){
		uout = uin;
		this->step_impl(uout,dudt,t,dt);
	}

};

} // namespace odeint
} // namespace geoflow

#endif /* SRC_ODEINT_EULER_STEPPER_HPP_ */
