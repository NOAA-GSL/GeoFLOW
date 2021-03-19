/*
 * observer_base.hpp
 *
 *  Created on: Nov 27, 2018
 *      Author: bflynt
 */

#ifndef SRC_PDEINT_OBSERVER_BASE_HPP_
#define SRC_PDEINT_OBSERVER_BASE_HPP_

#include <memory>
#include <vector>
#include <string>
#include "pdeint/equation_base.hpp"

namespace geoflow {
namespace pdeint {

/**
 * Base class to provided an interface for all observations.
 *
 * The base class provides the interface using the strategy pattern
 * for all observation types.  Observers take the solution after
 * every time step and can extract any user requested information,
 * display, save states, etc.
 */
template<typename TypePack>
class ObserverBase {

public:
        enum ObsType     {OBS_CYCLE=0, OBS_TIME};
	using Types         = TypePack;
        using EqnBase       = EquationBase<TypePack>;
        using EqnBasePtr    = std::shared_ptr<EqnBase>;
	using State         = typename Types::State;
        using StateInfo     = typename Types::StateInfo;
	using Grid          = typename Types::Grid;
	using Value         = typename Types::Value;
	using Derivative    = typename Types::Derivative;
	using Time          = typename Types::Time;
	using Jacobian      = typename Types::Jacobian;
	using Size          = typename Types::Size;
//      using ObsTraits     = ObserverBase<Types>::Traits;

        /**
         * Data structure to describe quantities derived from state data
         */
        struct dqTraits {
                std::vector<int>    icomponents;      // which state components to use     
                std::vector<std::string>   
                                    snames;           // names for computed components
                std::string         agg_sname;        // accumulated derived quantity filename
                std::string         smath_op;         // which math op to use to compute
        };

        /**
         * Data structure to hold user selected parameters
         */
        struct Traits {
                bool      treat_as_1d    = false;     // treat obs in 1d sense
                int       imisc          = 0;         // obs output type
                int       itag1          = 0;         // integer tag
                int       itag2          = 0;         // integer tag
                int       itag3          = 0;         // integer tag
                ObsType   itype          = OBS_CYCLE; // obs output type
                std::vector<int>     
                          state_index;                // which state members to observe
                std::vector<std::string>
                          state_names;                // file/ref names for each state member
                std::string 
                          agg_state_name          ;   // accumulated state filename
                std::vector<std::string>
                          grid_names;                 // file/ref names for each grid comp
                std::string 
                          agg_grid_name          ;    // accumulated grid filename
                std::vector<dqTraits>
                          derived_quantities;         // derived types: [ [type, function, output name],..]
                std::vector<dqTraits>
                          state_derived_quantities;   // derived types: [ [type, function, output name],..]
                size_t    start_ocycle   = 0  ;       // starting output cycle 
//              size_t    start_cycle    = 0  ;       // starting evol cycle 
                size_t    cycle_interval = 10 ;       // cycle interval for observation
                double    start_time     = 0.0;       // starting evol time
                double    time_interval  = 1.0;       // time interval for observation
                double    freq_fact      = 1.0;       // ties output freq to another observer (e.g., restart)
                std::string
                          idir                ;       // input directory (for data e.g.)
                std::string
                          odir                ;       // input directory (for I/O, e.g.)
                std::string
                          stag1               ;       // string tag
                std::string
                          stag2               ;       // string tag
                std::string
                          stag3               ;       // string tag
        };

        ObserverBase() = default;
	ObserverBase(EqnBasePtr& equation, Grid& grid, Traits& traits) {eqn_ptr_=equation; traits_=traits; grid_= &grid; utmp_=nullptr;}
	ObserverBase(ObserverBase& obs) = default;
	virtual ~ObserverBase() = default;
	ObserverBase& operator=(const ObserverBase& obs) = default;

	/**
	 * Observe the current state at time
	 *
	 * @param[in] Traits structure
	 * @param[in] t Time of current state
	 * @param[in] u State at the current time
	 */
	void observe(const Time t, const State& u, const State& uf){
		this->observe_impl(t,u,uf);
	}

	/**
	 * Do initialization
	 *
	 */
	void init(StateInfo &info){
		init_impl(info);
        } 
	/**
	 * Set tmp space
	 *
	 * @param[in] State variable for tmp
	 */
	void set_tmp(State& utmp){
		utmp_ = &utmp;
        } 

        /**
         * Get traits.
         *
         */
        Traits &get_traits() {return traits_;}

protected:
        Traits       traits_;
        Grid        *grid_;
        State       *utmp_;
        EqnBasePtr   eqn_ptr_;
	/**
	 * Must be provided by implementation
	 */
	virtual void observe_impl(const Time& t, const State& u, const State& uf) = 0;
	virtual void init_impl   (StateInfo &) = 0;

};

} // namespace pdeint
} // namespace geoflow

#endif /* SRC_PDEINT_OBSERVER_BASE_HPP_ */
