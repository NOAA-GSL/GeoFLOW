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


	struct Traits {
		Value cfl_min;
		Value cfl_max;
		Time  dt_min;
		Time  dt_max;
	};

	Integrator() = delete;

	Integrator(const StepBasePtr& stepper,
			   const ObsBasePtr&  observer,
			   const Traits& traits);

	Integrator(const Integrator& I) = default;
	~Integrator() = default;
	Integrator& operator=(const Integrator& I) = default;


	void time( const Time& t0,
			   const Time& t1,
			   const Time& dt,
			   State&      u );


	void steps( const Time& t0,
			    const Time& dt,
				const Size& n,
			    State&      u,
				Time&       t );

	void list( const std::vector<Time>& tlist,
		       State&                   u );



protected:
	Traits      traits_;
	StepBasePtr stp_ptr_;
	ObsBasePtr  obs_ptr_;

	void init_dt(const Time& t, State& u, Time& dt) const;

};


} // namespace odeint
} // namespace geoflow


#include "odeint/integrator.ipp"

#endif /* SRC_ODEINT_INTEGRATOR_HPP_ */
