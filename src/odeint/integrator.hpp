/*
 * integrator.hpp
 *
 *  Created on: Nov 26, 2018
 *      Author: bryan.flynt
 */

#ifndef SRC_ODEINT_INTEGRATOR_HPP_
#define SRC_ODEINT_INTEGRATOR_HPP_


#include <memory>
#include <vector>

#include "odeint/null_observer.hpp"
#include "odeint/stepper_base.hpp"
#include "stepper_base.hpp"

namespace geoflow {
namespace odeint {

/**
 * Integrator to perform multiple time steps on system of equations.
 *
 * The Integrator class performs multiple time steps and observations
 * to a provided system of equations.
 *
 * @see EquationBase
 * @see StepperBase
 */
template<typename EquationType>
class Integrator {

public:
	using Equation     = EquationType;
	using State        = typename Equation::State;
	using Value        = typename Equation::Value;
	using Derivative   = typename Equation::Derivative;
	using Time         = typename Equation::Time;
	using Jacobian     = typename Equation::Jacobian;
	using Size         = typename Equation::Size;
	using StepBasePtr  = std::shared_ptr<StepperBase<Equation>>;
	using ObsBasePtr   = std::shared_ptr<ObserverBase<Equation>>;

	/**
	 * Data structure to hold user selected parameters
	 */
	struct Traits {
		Value cfl_min;
		Value cfl_max;
		Time  dt_min;
		Time  dt_max;
	};

	Integrator() = delete;

	/**
	 * Constructor to initialize everything needed to run
	 *
	 * @param[in] stepper  Pointer to the time stepping algorithm
	 * @param[in] observer Pointer to the observer to use
	 * @param[in] traits   Use selected traits for time step options
	 */
	Integrator(const StepBasePtr& stepper,
			   const ObsBasePtr&  observer,
			   const Traits& traits);

	Integrator(const Integrator& I) = default;
	~Integrator() = default;
	Integrator& operator=(const Integrator& I) = default;

	/**
	 * Take as many steps required to progress from time t0 to t1.
	 *
	 * The time method takes as many steps required to progress
	 * from time t0 to time t1.  It will attempt to use the provided
	 * time step size dt unless one of the user provided traits limits
	 * the step size or a smaller time step is required to terminate
	 * at time t1.
	 *
	 * @param[in]     t0 Initial time at start
	 * @param[in]     t1 Final time an completion of time stepping
	 * @param[in]     dt Recommend time step size
	 * @param[in,out] u  Current and final equation state values
	 */
	void time( const Time& t0,
			   const Time& t1,
			   const Time& dt,
			   State&      u );

	/**
	 * Take an exact number of time steps.
	 *
	 * The steps method takes a specified number of time steps
	 * regardless of the final time it arrives at.  It will
	 * attempt to use the provided time step size unless one
	 * of the user provided traits limits the step size.
	 *
	 * @param[in]     t0 Initial time at start
	 * @param[in]     dt Recommend time step size
	 * @param[in]     n  Number of steps to take
	 * @param[in,out] u  Current and final equation state values
	 * @param[out]    t  Final time resulting from taking n steps
	 */
	void steps( const Time& t0,
			    const Time& dt,
				const Size& n,
			    State&      u,
				Time&       t );

	/**
	 * Take steps that correspond to the provided list of times
	 *
	 * The list method takes steps required to reach the provided
	 * list of time steps.  The method attempt to reach the next time
	 * using one time step equal to the difference between the times.
	 * If the user provided traits limits the step size then multiple
	 * steps will be taken to reach the listed times exactly.
	 *
	 * @param[in]     tvec Vector of time points to march through
	 * @param[in,out] u  Current and final equation state values
	 */
	void list( const std::vector<Time>& tvec,
		       State&                   u );


protected:
	Traits      traits_;
	StepBasePtr stp_ptr_;
	ObsBasePtr  obs_ptr_;

	/**
	 * Used to calculate the limited time step size if required.
	 */
	void init_dt(const Time& t, State& u, Time& dt) const;

};


} // namespace odeint
} // namespace geoflow


#include "odeint/integrator.ipp"

#endif /* SRC_ODEINT_INTEGRATOR_HPP_ */
