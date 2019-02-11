/*
 * equation_base.hpp
 *
 *  Created on: Nov 16, 2018
 *      Author: bflynt
 */

#ifndef SRC_PDEINT_EQUATION_BASE_HPP_
#define SRC_PDEINT_EQUATION_BASE_HPP_


namespace geoflow {
namespace pdeint {



/**
 * Base class to provided an interface for all equation types.
 *
 * The base class provides the interface using the strategy pattern
 * for all equation implementations.
 */
template<typename TypePack>
class EquationBase {

public:
	using Types      = TypePack;
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
	 * using the modeled equation.
	 *
	 * @param[in]     t  Current time of state u before taking step
	 * @param[in,out] u  Is the state of the system of equations
	 * @param[out]    dt Size of time step corresponding to CFL of one
	 */
	void dt(const Time& t, State& u, Time& dt){
		this->dt_impl(t,u,dt);
	}

	/**
	 * Apply Boundary Condition.
	 *
	 * Apply the provided ub boundary condition to the state u.
	 *
	 * @param[in]     t  Current time of state u before taking step
	 * @param[in,out] u  Is the state of the system of equations
	 * @param[in]     ub Boundary condition values
	 */
	void apply_bc(const Time &t, State &u, State &ub){
		this->apply_bc_impl(t,u,ub);
	}


	/** Take one step from time t to time t+dt.
	 *
	 * Take exactly one time step from t with current solution u to
	 * time t + dt and return new solution within u.
	 *
	 * \param[in]     t  Current time of state u before taking step
	 * \param[in]     dt Size of time step to take
	 * \param[in]     ub Boundary conditions for u
	 * \param[in,out] u  State of the system of equations
	 */
	void step(const Time& t, const Time& dt, State& u, State& ub){
		this->step_impl(t,dt,u,ub);
	}

	/** Take one step from time t to time t + dt.
	 *
	 * Take exactly one time step from t with current solution uin to
	 * time t + dt and return new solution within uout.
	 *
	 * \param[in]  t     Current time of state u before taking step
	 * \param[in]  uin   State of the system of equations at time t
	 * \param[in]  ub    Boundary conditions for u
	 * \param[in]  dt    Size of time step to take
	 * \param[out] uout  New state of the system of equations at t + dt
	 */
	void step(const Time& t, const State& uin, const State& ub, const Time& dt, State& uout){
		this->step_impl(t,uin,ub,dt,uout);
	}

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
	virtual void apply_bc_impl(const Time &t, State &u, State &ub) = 0;

	/**
	 * Must be provided by implementation
	 */
	virtual void step_impl(const Time& t, const Time& dt, State& u, State& ub) = 0;

	/**
	 * Must be provided by implementation
	 */
	virtual void step_impl(const Time& t, const State& uin, const State& ub, const Time& dt, State& uout) = 0;
};


} // namespace pdeint
} // namespace geoflow


#endif /* SRC_PDEINT_EQUATION_BASE_HPP_ */
