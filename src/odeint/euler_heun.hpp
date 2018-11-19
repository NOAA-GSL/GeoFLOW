/*
 * euler_heun_error_stepper.hpp
 *
 *  Created on: Nov 18, 2018
 *      Author: bflynt
 */

#ifndef SRC_ODEINT_ERROR_STEPPER_EULER_HEUN_HPP_
#define SRC_ODEINT_ERROR_STEPPER_EULER_HEUN_HPP_



#include <cstddef>
#include <memory>

#include "odeint/error_stepper/error_stepper_base.hpp"

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
class EulerHeun : public ErrorStepperBase<EquationType,ErrorType> {

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

	EulerHeun() = default;
	EulerHeun(const EulerHeun& si) = default;
	virtual ~EulerHeun() = default;
	EulerHeun& operator=(const EulerHeun& si) = default;


protected:

	Size order_impl() const = 0{
		return 2;
	}

	void step_impl(EquationPtr& sys, State& u, const Time& t, const Time& dt, Error& uerr){
		Derivative dudt(u.size());
		eqn->dudt(u, dudt, t);
		this->step_impl(eqn,u,dudt,t,dt,uerr);
	}

	void step_impl(EquationPtr& sys, const State& uin, const Time& t, State& uout, const Time& dt, Error& uerr){
		uout = uin;
		this->step_impl(eqn,uout,t,dt,uerr);
	}

	void step_impl(EquationPtr& sys, State& u, const Derivative& dudt, const Time& t, const Time& dt, Error& uerr){
		u += dt * dudt;
		Derivative dudt1(dudt.size());
		eqn->dudt(u, dudt1, t+dt);
		uerr = 0.5 * dt * (dudt1 - dudt);
		u += uerr;
	}

	void step_impl(EquationPtr& sys, const State& uin, const Derivative& dudt, const Time& t, State& uout, const Time& dt, Error& uerr){
		uout = uin;
		this->step_impl(sys,uout,dudt,t,dt,uerr);
	}

};

} // namespace odeint
} // namespace geoflow


#endif /* SRC_ODEINT_ERROR_STEPPER_EULER_HEUN_HPP_ */
