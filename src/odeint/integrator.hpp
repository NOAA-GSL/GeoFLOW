/*
 * integrator.hpp
 *
 *  Created on: Nov 26, 2018
 *      Author: bryan.flynt
 */

#ifndef SRC_ODEINT_INTEGRATOR_HPP_
#define SRC_ODEINT_INTEGRATOR_HPP_


#include <memory>

namespace geoflow {
namespace odeint {

template<typename StepperType>
struct Integrator {
	using Stepper     = StepperType;
	using Equation    = typename Stepper::Equation;
	using State       = typename Stepper::State;
	using Value       = typename Stepper::Value;
	using Derivative  = typename Stepper::Derivative;
	using Time        = typename Stepper::Time;
	using Jacobian    = typename Stepper::Jacobian;
	using EquationPtr = typename Stepper::EquationPtr;
	using StepperPtr  = std::shared_ptr<Stepper>;
	using Size        = typename Stepper::Size;


	static void time(StepperPtr&  stp_ptr,
			         EquationPtr& eqn_ptr,
					 State& u,
					 Time& t0,
					 Time& t1,
					 Time& dt
					 ) //, ObserverPtr& obs_ptr = std::nullptr);
	{

		// Build list of times
		const Size sz = static_cast<Size>((t1-t0)/dt);
		std::vector<Time> tlist(sz);
		for(Size i = 0; i < sz; ++i){
			tlist[i] = t0 + (i * dt);
		}
		tlist.push_back(t1);

		// Call the list integrator
		list(stp_ptr,eqn_ptr,u,tlist);
	}

	static void step(StepperPtr&  stp_ptr,
			         EquationPtr& eqn_ptr,
					 State& u0,
					 Time&  t0,
					 Time&  dt,
					 Size&  n
					 ){ //, ObserverPtr& obs_ptr = std::nullptr);

		for(Size i = 0; i < n; ++i){

			do{

			} while(true);

		}


	}

	static void list(StepperPtr&  stp_ptr,
			         EquationPtr& eqn_ptr,
					 State& u0,
					 std::vector<Time>&  tlist,
					 ); //, ObserverPtr& obs_ptr = std::nullptr);

};


} // namespace odeint
} // namespace geoflow

#endif /* SRC_ODEINT_INTEGRATOR_HPP_ */
