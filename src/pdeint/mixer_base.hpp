/*
 * mixer_base.hpp
 *
 *  Created on: Mar 27, 2019
 *      Author: bflynt, d.rosemnerg
 */

#ifndef SRC_PDEINT_MIXER_BASE_HPP_
#define SRC_PDEINT_MIXER_BASE_HPP_

#include <memory>
#include <vector>
#include <string>
#include "equation_base.hpp"

namespace geoflow {
namespace pdeint {

/**
 * Base class to provided an interface for all 'mixrers'
 *
 * The base class provides the interface using the strategy pattern
 * for all mixrer types.  Mixers take the forcing after
 * every time step and determine if it should 'randomized', and
 * how.
 */
template<typename TypePack>
class MixerBase {

public:
        enum MixType     {STIR_CYCLE=0, STIR_TIME};
	using Types       = TypePack;
        using EqnBase     = EquationBase<TypePack>;
        using EqnBasePtr  = std::shared_ptr<EqnBase>;
        using State       = typename Types::State;
	using Grid        = typename Types::Grid;
	using Ftype       = typename Types::Ftype;
	using Derivative  = typename Types::Derivative;
	using Time        = typename Types::Time;
	using Jacobian    = typename Types::Jacobian;
	using Size        = typename Types::Size;

        /**
         * Data structure to hold user selected parameters
         */
        struct Traits {
                MixType  itype          = STIR_TIME; // mix cadence type
                std::vector<int>     
                          state_index;                // which state members to 'mix'
                size_t    corr_cycle     = 10;        // correlation cycle interval 
                double    corr_time      = 1.0;       // correlation time interval 
        };

        MixerBase() = default;
	MixerBase(Traits& traits, Grid& grid){traits_=traits; grid_= &grid;}
	MixerBase(const MixerBase& obs) = default;
	virtual ~MixerBase() = default;
	MixerBase& operator=(const MixerBase& obs) = default;

	/**
	 * Mix the forcing state at time
	 *
	 * @param[in] t Time of current state
	 * @param[in] u State at the current time
	 * @param[in] uf forcing tendency at the current time
	 */
	void mix(const Time t, const State& u, State& uf){
		this->mix_impl(t,u,uf);
	}

        /**
         * Set tmp space
         *
         * @param[in] State variable for tmp
         */
        void set_tmp(State& utmp){
                utmp_ = &utmp;
        }


protected:
        Traits traits_;
        Grid  *grid_;
        State *utmp_;
	/**
	 * Must be provided by implementation
	 */
	virtual void mix_impl(const Time& t, const State& u, State &uf) = 0;

};

} // namespace pdeint
} // namespace geoflow

#endif /* SRC_PDEINT_MIXER_BASE_HPP_ */
