/*
 * integrate_time.hpp
 *
 *  Created on: Nov 18, 2018
 *      Author: bflynt
 */

#ifndef SRC_ODEINT_INTEGRATE_TIME_HPP_
#define SRC_ODEINT_INTEGRATE_TIME_HPP_




#include <memory>

namespace geoflow {
namespace odeint {

template<typename EquationType, typename StepperType>
class TimeIntegrator {

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

	TimeIntegrator() = default;
	TimeIntegrator(const TimeIntegrator& ti) = default;
	virtual ~TimeIntegrator() = default;
	TimeIntegrator& operator=(const TimeIntegrator& ti) = default;

	//void setEquation(EquationPtr& eqn);

	//void setStepper(StepperPtr& stp);

	//void setObserver(ObserverPtr& obs);

protected:

	void integrate_impl(std::shared_ptr<StepperBase>, std::shared_ptr<EquationBase>, State& u0, const Time& t0, const Time& t1, const Time& dt);

	void integrate_impl(std::shared_ptr<ErrorStepperBase>, std::shared_ptr<EquationBase>, State& u0, const Time& t0, const Time& t1, const Time& dt);

	void integrate_impl(std::shared_ptr<ControlStepperBase>, std::shared_ptr<EquationBase>, State& u0, const Time& t0, const Time& t1, const Time& dt);

	void integrate_impl(std::shared_ptr<StepperBase>, std::shared_ptr<EquationBase>, State& u0, const Time& t0, const Time& t1, const Time& dt);

};

} // namespace odeint
} // namespace geoflow



#endif /* SRC_ODEINT_INTEGRATE_TIME_HPP_ */
