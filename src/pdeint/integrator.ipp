/*
 * integrator.ipp
 *
 *  Created on: Nov 27, 2018
 *      Author: bflynt
 */

#include <algorithm>
#include <limits>

#include "integrator.hpp"
#include "tbox/assert.hpp"
#include "tbox/error_handler.hpp"

namespace geoflow {
namespace pdeint {

template<typename EquationType>
Integrator<EquationType>::Integrator(const EqnBasePtr& equation,
		                     const ObsBasePtr& observer,
                                     Grid&             grid,
		                     const Traits&     traits) :
	cycle_(0), traits_(traits), eqn_ptr_(equation), obs_ptr_(observer), grid_(&grid) {
	ASSERT(nullptr != eqn_ptr_);
	ASSERT(nullptr != obs_ptr_);
}

template<typename EquationType>
void Integrator<EquationType>::time_integrate( Time&       t,
                                               State&      ub,
		                               State&      u ){
	ASSERT(nullptr != eqn_ptr_);
	ASSERT(nullptr != obs_ptr_);
	using std::min;

        Time t0 = t;
        Size n  = traits_.cycle_end - cycle_ + 1;
        if      ( traits_.integ_type == INTEG_CYCLE ) {
          steps( t0 
               , traits_.dt
               , n
               , ub
               , u
               , t);
        } 
        else if ( traits_.integ_type == INTEG_TIME  ) {
          time ( t
               , traits_.time_end
               , traits_.dt
               , ub
               , u);
          t = traits_.time_end;
        }

}


template<typename EquationType>
void Integrator<EquationType>::time( const Time& t0,
		                     const Time& t1,
		                     const Time& dt,
		                     State&      ub,
		                     State&      u ){
	ASSERT(nullptr != eqn_ptr_);
	ASSERT(nullptr != obs_ptr_);
	using std::min;

	Time t = t0;
	do {

		// Limit dt if requested
		Time dt_eff = dt;
		this->init_dt(t,u,dt_eff);
		dt_eff = min(dt_eff, t1 - t);
		if(0 >= dt_eff) {
			EH_ERROR("Effective Time Step is 0");
		}

		// Take Step
		this->eqn_ptr_->step(t,u, ub, dt_eff);
		t += dt_eff;


		// Call observer on solution
                for ( auto j = 0 ; j < this->obs_ptr_->size(); j++ ) 
		  (*this->obs_ptr_)[j]->observe(t,u);
 
                ++cycle_;

	} while( t != t1 );


}

template<typename EquationType>
void Integrator<EquationType>::steps( const Time&  t0,
		                      const Time&  dt,
			              const Size&  n,
		                      State&       ub,
		                      State&       u,
			              Time&        t ){
	ASSERT(nullptr != eqn_ptr_);
	ASSERT(nullptr != obs_ptr_);

	t = t0;
	for(Size i = 0; i < n; ++i, ++cycle_){

		// Limit dt if requested
		Time dt_eff = dt;
		this->init_dt(t,u,dt_eff);
		if(0 >= dt_eff) {
			EH_ERROR("Effective Time Step is 0");
		}

		// Take Step
		this->eqn_ptr_->step(t,u, ub, dt_eff);
		t += dt_eff;

		// Call observer on solution at new t
                for ( auto j = 0 ; j < this->obs_ptr_->size(); j++ ) 
		  (*this->obs_ptr_)[j]->observe(t,u);
	}

}

template<typename EquationType>
void Integrator<EquationType>::list( const std::vector<Time>& tlist,
		                     State&                   ub,
	                             State&                   u ){
	ASSERT(nullptr != eqn_ptr_);
	ASSERT(nullptr != obs_ptr_);

	for(Size i = 0; i < tlist.size()-1; ++i){

		// Calculate dt
		Time dt = tlist[i+1] - tlist[i];

		// Step till next time
		this->time(tlist[i],tlist[i+1],dt,u);
	}

}


template<typename EquationType>
void Integrator<EquationType>::init_dt(const Time& t, State& u, Time& dt) const{
	using std::max;
	using std::min;

	// Get pointer to equation
	auto eqn_ptr = this->eqn_ptr_; // ->getEquationPtr();
	ASSERT(nullptr != eqn_ptr);

	// Make sure dt is between absolute limits
	dt = min(dt, traits_.dt_max);
	dt = max(dt, traits_.dt_min);

	// Make sure dt is between CFL limits
	if( eqn_ptr->has_dt() ){
		Time dt_cfl_1 = std::numeric_limits<Time>::max()/(1+traits_.cfl_max);
		eqn_ptr->dt(t,u,dt_cfl_1);
		dt = min(dt, dt_cfl_1*traits_.cfl_max);
		dt = max(dt, dt_cfl_1*traits_.cfl_min);
	}

}



} // namespace pdeint
} // namespace geoflow
