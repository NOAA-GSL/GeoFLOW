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

#include "odeint/error_stepper_base.hpp"

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
class EulerHeun : public ErrorStepperBase<EquationType> {

public:
	// Inheriting these from the Base we insure our
	// interface is consistent and we can perform runtime
	// polymorphism or use them directly.
	using Equation    = EquationType;
	using Interface   = ErrorStepperBase<EquationType>;
	using State       = typename Interface::State;
	using Value       = typename Interface::Value;
	using Derivative  = typename Interface::Derivative;
	using Time        = typename Interface::Time;
	using Jacobian    = typename Interface::Jacobian;
	using Size        = typename Interface::Size;
	using EquationPtr = std::shared_ptr<Equation>;

	EulerHeun() = default;
	EulerHeun(const EulerHeun& si) = default;
	virtual ~EulerHeun() = default;
	EulerHeun& operator=(const EulerHeun& si) = default;

	EulerHeun(const EquationPtr& eqn) : Interface(eqn){
	}

protected:

	Size order_impl() const = 0{
		return 2;
	}

	void step_impl(State& u, const Time& t, const Time& dt, Error& uerr){
		Derivative dudt(u.size());
		this->eqn_ptr_->dudt(u, dudt, t);
		this->step_impl(u,dudt,t,dt,uerr);
	}

	void step_impl(State& uin, const Time& t, State& uout, const Time& dt, Error& uerr){
		uout = uin;
		this->step_impl(uout,t,dt,uerr);
	}

	void step_impl(State& u, const Derivative& dudt, const Time& t, const Time& dt, Error& uerr){
		u += dt * dudt;
		Derivative dudt1(dudt.size());
		this->eqn_ptr_->dudt(u, dudt1, t+dt);
		uerr = 0.5 * dt * (dudt1 - dudt);
		u += uerr;
	}

	void step_impl(State& uin, const Derivative& dudt, const Time& t, State& uout, const Time& dt, Error& uerr){
		uout = uin;
		this->step_impl(uout,dudt,t,dt,uerr);
	}

};

} // namespace odeint
} // namespace geoflow


#endif /* SRC_ODEINT_ERROR_STEPPER_EULER_HEUN_HPP_ */
