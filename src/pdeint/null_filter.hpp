/*
 * null_filter.hpp
 *
 *  Created on: Sep. 14, 2020
 *      Author: d.rosenberg
 */

#ifndef SRC_PDEINT_NULL_FILTER_HPP_
#define SRC_PDEINT_NULL_FILTER_HPP_


#include "pdeint/filter_base.hpp"

namespace geoflow {
namespace pdeint {

/**
 * Do nothing filter implementation
 *
 * The NullFilter is the default filter when nothing else is
 * provided and does nothing when called.
 *
 * @see FilterBase
 */
template<typename TypePack>
class NullFilter : public FilterBase<TypePack> {

public:
        using Interface  = FilterBase<TypePack>;
        using State      = typename Interface::State;
        using StateComp  = typename Interface::StateComp;
        using Grid       = typename Interface::Grid;
        using Time       = typename Interface::Time;

        NullFilter() {}
        
#if 0
        NullFilter(const NullFilter& eb) = default;
       ~NullFilter() = default;
        NullFilter& operator=(const NullFilter& eb) = default;
#endif


protected:

	void apply_impl(const Time& /* t */, State& /* u */, State& /* utmp */, State& /* uo */){
		// Do nothing ...
	}

	void apply_impl(const Time& /* t */, State& /* u */, State& /* utmp */){
		// Do nothing ...
	}

};


} // namespace pdeint
} // namespace geoflow



#endif /* SRC_PDEINT_NULL_FILTER_HPP_ */
