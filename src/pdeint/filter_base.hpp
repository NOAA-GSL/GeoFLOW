/*
 * filter_base.hpp
 *
 *  Created on: Sep. 14, 2020
 *      Author: d.rosenberg
 */

#ifndef SRC_PDEINT_FILTER_BASE_HPP_
#define SRC_PDEINT_FILTER_BASE_HPP_

#include <functional>


namespace geoflow {
namespace pdeint {

/**
 * Base class to provided an interface for all (state) filter types.
 *
 * The base class provides the interface using the strategy pattern
 * for all equation implementations.
 */
template<typename TypePack>
class FilterBase {

public:
	using Types      = TypePack;
	using State      = typename Types::State;
	using StateComp  = typename Types::StateComp;
        using StateInfo  = typename Types::StateInfo;
	using Grid       = typename Types::Grid;
	using Ftype      = typename Types::Ftype;
	using Time       = typename Types::Time;
	using Size       = typename Types::Size;


        FilterBase() = default;
	FilterBase(const FilterBase& eb) = default;
	virtual ~FilterBase() = default;
	FilterBase& operator=(const FilterBase& eb) = default;

	/**
	 * Apply flter.
	 *
	 * Apply the filter to the state u.
	 *
	 * @param[in]      t     Current time of state u before taking step
	 * @param[in]      u     Is the state of the system of equations
	 * \param[in/out]  utmp  tmp state vectors
	 * @param[out]     uo    filtered state 
	 */
	void apply(const Time &t, State& u, State& utmp, State& uo){
		this->apply_impl(t,u,utmp,uo);
	}

	/**
	 * Apply filter, in-place.
	 *
	 * Apply the filter to the state u, in place.
	 *
	 * @param[in]     t     Current time of state u before taking step
	 * @param[in,out] u     Is the state of the system of equations
	 * \param[in]     utmp  tmp state vectors
	 */
	void apply(const Time &t, State& u, State& utmp){
		this->apply_impl(t,u,utmp);
	}


protected:

	/**
	 * Must be provided by implementation
	 */
	virtual void apply_impl(const Time& t, State& u, State& utmp, State& uo) = 0;
	virtual void apply_impl(const Time& t, State& u, State& utmp) = 0;

};


} // namespace pdeint
} // namespace geoflow


#endif /* SRC_PDEINT_FILTER_BASE_HPP_ */
