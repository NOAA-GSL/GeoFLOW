/*
 * integrator.hpp
 *
 *  Created on: Nov 26, 2018
 *      Author: bryan.flynt
 */

#ifndef SRC_PDEINT_INTEGRATOR_HPP_
#define SRC_PDEINT_INTEGRATOR_HPP_


#include <memory>
#include <vector>

#include "pdeint/null_observer.hpp"
#include "pdeint/equation_base.hpp"

namespace geoflow {
namespace pdeint {

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
        enum IntegType     {INTEG_CYCLE=0, INTEG_TIME};
	using Equation     = EquationType;
        using EqnBase      = EquationBase<EquationType>;
	using State        = typename Equation::State;
	using Grid         = typename Equation::Grid;
	using Value        = typename Equation::Value;
	using Derivative   = typename Equation::Derivative;
	using Time         = typename Equation::Time;
	using Jacobian     = typename Equation::Jacobian;
	using Size         = typename Equation::Size;
	using EqnBasePtr   = std::shared_ptr<EqnBase>;
	using ObsBasePtr   = std::shared_ptr<std::vector<std::shared_ptr<ObserverBase<Equation>>>>;
      
	/**
	 * Data structure to hold user selected parameters
	 */
	struct Traits {
		IntegType integ_type;
		size_t    cycle_end   = 1;
		Value     cfl_min     = std::numeric_limits<Value>::min();
		Value     cfl_max     = std::numeric_limits<Value>::max();
		Time      dt_min      = std::numeric_limits<Time>::min();
		Time      dt_max      = std::numeric_limits<Time>::max();   
		Time      dt          = static_cast<Time>(1.0e-2);
		Time      time_end    = static_cast<Time>(1.0);
	};

	Integrator() = delete;

	/**
	 * Constructor to initialize everything needed to run
	 *
	 * @param[in] stepper  Pointer to the time stepping algorithm
	 * @param[in] observer Pointer to the observer to use
	 * @param[in] traits   Use selected traits for time step options
	 */
	Integrator(const EqnBasePtr&  equation,
		   const ObsBasePtr&  observer,
                         Grid&        grid,
		   const Traits&      traits); 

	Integrator(const Integrator& I) = default;
	~Integrator() = default;
	Integrator& operator=(const Integrator& I) = default;

	/**
	 * Integrate state at t, as required by prop tree input
	 *
	 * The time_integrate method integrates from current state 
         * using either 'time' or 'steps' methods according to the 
         * property tree input. The state time, t, will be updated
         * with the final integration time.
	 *
	 * @param[in,out] t  Initial time at start, and final time
	 * @param[in,out] uf Forcing tendency
	 * @param[in,out] ub Boundary condition vectors
	 * @param[in,out] u  Current and final equation state values
	 */
	void time_integrate( Time&        t,
                             State&       uf,
                             State&       ub,
			     State&       u );
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
	 * @param[in,out] uf Forcing tendency
	 * @param[in,out] ub Boundary condition vectors
	 * @param[in,out] u  Current and final equation state values
	 */
	void time( const Time&  t0,
		   const Time&  t1,
		   const Time&  dt,
                   State&       uf,
                   State&       ub,
		   State&       u );

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
	 * @param[in,out] uf Forcing tendency
	 * @param[in,out] ub Boundary condition vectors
	 * @param[in,out] u  Current and final equation state values
	 * @param[out]    t  Final time resulting from taking n steps
	 */
	void steps( const Time&  t0,
		    const Time&  dt,
		    const Size&  n,
                    State&       uf,
                    State&       ub,
		    State&       u,
	    	    Time&        t );

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
	 * @param[in,out] uf Forcing tendency
	 * @param[in,out] ub Boundary condition vectors
	 * @param[in,out] u  Current and final equation state values
	 */
	void list( const std::vector<Time>& tvec,
                   State&                   uf,
                   State&                   ub,
		   State&                   u );


protected:
        size_t      cycle_; // no. time cycles taken so far
	Traits      traits_;
	EqnBasePtr  eqn_ptr_;
	ObsBasePtr  obs_ptr_;
        Grid*       grid_;

	/**
	 * Used to calculate the limited time step size if required.
	 */
	void init_dt(const Time& t, State& u, Time& dt) const;

};


} // namespace pdeint
} // namespace geoflow


#include "pdeint/integrator.ipp"

#endif /* SRC_PDEINT_INTEGRATOR_HPP_ */
