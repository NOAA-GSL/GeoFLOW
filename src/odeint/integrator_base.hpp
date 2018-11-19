/*
 * integrator_base.hpp
 *
 *  Created on: Nov 18, 2018
 *      Author: bflynt
 */

#ifndef SRC_ODEINT_INTEGRATOR_BASE_HPP_
#define SRC_ODEINT_INTEGRATOR_BASE_HPP_



template<typename EquationType, typename StepperType, typename ObserverType>
class IntegratorBase {

public:
	using Stepper     = StepperType;
	using Equation    = EquationType;
	using State       = typename Stepper::State;
	using Value       = typename Stepper::Value;
	using Derivative  = typename Stepper::Derivative;
	using Time        = typename Stepper::Time;
	using Jacobian    = typename Stepper::Jacobian;
	using EquationPtr = typename Stepper::EquationPtr;
	using StepperPtr  = std::shared_ptr<Stepper>;
	using Size        = std::size_t;

	IntegratorBase() = default;
	IntegratorBase(const IntegratorBase& ti) = default;
	virtual ~IntegratorBase() = default;
	IntegratorBase& operator=(const IntegratorBase& ti) = default;

	void setEquation(EquationPtr& eqn){
		eqn_ptr_ = eqn;
	}

	void setStepper(StepperPtr& stp){
		stp_ptr_ = stp;
	}

	void setObserver(ObserverPtr& obs){
		obs_ptr_ = obs;
	}

	void integrate(State& u0, const Time& t0, const Time& t1, const Time& dt){
		this->integrate_impl(stp_ptr_,eqn_ptr_,u0,t0,t1,dt);
	}

protected:
	EquationPtr eqn_ptr_;
	StepperPtr  stp_ptr_;
	ObserverPtr obs_ptr_;


	virtual void integrate_impl(std::shared_ptr<StepperBase>, std::shared_ptr<EquationBase>, State& u0, const Time& t0, const Time& t1, const Time& dt) = 0;

	virtual void integrate_impl(std::shared_ptr<ErrorStepperBase>, std::shared_ptr<EquationBase>, State& u0, const Time& t0, const Time& t1, const Time& dt) = 0;

	virtual void integrate_impl(std::shared_ptr<ControlStepperBase>, std::shared_ptr<EquationBase>, State& u0, const Time& t0, const Time& t1, const Time& dt) = 0;

	virtual void integrate_impl(std::shared_ptr<StepperBase>, std::shared_ptr<EquationBase>, State& u0, const Time& t0, const Time& t1, const Time& dt) = 0;

};


#endif /* SRC_ODEINT_INTEGRATOR_BASE_HPP_ */
