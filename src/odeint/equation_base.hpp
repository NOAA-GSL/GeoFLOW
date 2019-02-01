/*
 * equation_base.hpp
 *
 *  Created on: Nov 16, 2018
 *      Author: bflynt
 */

#ifndef SRC_ODEINT_EQUATION_BASE_HPP_
#define SRC_ODEINT_EQUATION_BASE_HPP_


namespace geoflow {
namespace odeint {



/**
 * Base class to provided an interface for all equation types.
 *
 * The base class provides the interface using the strategy pattern
 * for all equation implementations.
 */
template<typename TypePack>
class EquationBase {

public:
	using State      = typename TypePack::State;
	using Value      = typename TypePack::Value;
	using Derivative = typename TypePack::Derivative;
	using Time       = typename TypePack::Time;
	using Jacobian   = typename TypePack::Jacobian;
	using Size       = typename TypePack::Size;

	EquationBase() = default;
	EquationBase(const EquationBase& eb) = default;
	virtual ~EquationBase() = default;
	EquationBase& operator=(const EquationBase& eb) = default;

	/**
	 * Return if the equation system can determine a time step.
	 *
	 * Query if the equation system can return a time step value
	 * which is equivalent to a CFL of one.
	 *
	 * @return Boolean value
	 */
	bool has_dt() const{
		return this->has_dt_impl();
	}

	/**
	 * Time difference corresponding to a CFL of one.
	 *
	 * Calculate the time step size which equates to a CFL of one
	 * using the equations you're modeling.
	 *
	 * @param[in]     t  Current time of state u before taking step
	 * @param[in,out] u  Is the state of the system of equations
	 * @param[out]    dt Size of time step corresponding to CFL of one
	 */
	void dt(const Time& t, State& u, Time& dt){
		this->dt_impl(t,u,dt);
	}

	/**
	 * Calculate dudt for the system of equations.
	 *
	 * Using the time and current state calculate the
	 * derivative of the state with respect to time.
	 *
	 * @param[in]     t    Current time of state u before taking step
	 * @param[in,out] u    State of the system of equations at time t
	 * @param[out]    dudt Derivative of the state u with respect to time t
	 */
	void dudt(const Time& t, State& u, Derivative& dudt){
		this->dudt_impl(t,u,dudt);
	}

	/**
	 * Calculate dfdu for the system of equations
	 *
	 * @param[in]     t    Current time of state u before taking step
	 * @param[in,out] u    State of the system of equations at time t
	 * @param[out]    dfdu Derivative of the function with respect to time u
	 */
	void dfdu(const Time& t, State& u, Jacobian& dfdu){
		this->dfdu_impl(t,u,dfdu);
	}

	/**
<<<<<<< HEAD
	 * Apply boundary conditions to u at time t
	 *
         * @param[in]     t Current time of state u before taking step
         * @param[in,out] u State of the system of equations at time t
	 */
	bool apply_bc(const Time& t, State& u){
		return this->apply_bc_impl(t,u);
	}
=======
	 * Calculate non-linear component of dudt for the system of equations.
	 *
	 * Using the time and current state calculate the
	 * non-linear portion of the derivative for the state with
	 * respect to time.
	 *
	 * @param[in]     t    Current time of state u before taking step
	 * @param[in,out] u    State of the system of equations at time t
	 * @param[out]    dudt Non-Linear derivative of the state u with respect to time t
	 */
	void nonlinear_dudt(const Time& t, State &u, Derivative& dudt){
		this->nonlinear_dudt_impl(t,u,nlin);
	}

	/**
	 * Calculate linear component of dudt for the system of equations.
	 *
	 * Using the time and current state calculate the
	 * linear portion of the derivative for the state with
	 * respect to time.
	 *
	 * @param[in]     t    Current time of state u before taking step
	 * @param[in,out] u    State of the system of equations at time t
	 * @param[out]    dudt Linear derivative of the state u with respect to time t
	 */
	void linear_dudt(const Time& t, State &u, Derivative& dudt){
		this->linear_dudt_impl(t,u,lin);
	}

>>>>>>> develop

protected:

	/**
	 * Must be provided by implementation
	 */
	virtual bool has_dt_impl() const = 0;

	/**
	 * Must be provided by implementation
	 */
	virtual void dt_impl(const Time& t, State& u, Time& dt) = 0;

	/**
	 * Must be provided by implementation
	 */
	virtual void dudt_impl(const Time& t, State& u, Derivative& dudt) = 0;

	/**
	 * Must be provided by implementation
	 */
	virtual void dfdu_impl(const Time& t, State& u, Jacobian& dfdu) = 0;

	/**
	 * Must be provided by implementation
	 */
<<<<<<< HEAD
	virtual void apply_bc_impl(const Time& t, State& u) = 0;
=======
	virtual void nonlinear_dudt_impl(const Time& t, State &u, Derivative& dudt) = 0;

	/**
	 * Must be provided by implementation
	 */
	virtual void linear_dudt_impl(const Time& t, State &u, Derivative& dudt) = 0;

>>>>>>> develop
};


} // namespace odeint
} // namespace geoflow


#endif /* SRC_ODEINT_EQUATION_BASE_HPP_ */
