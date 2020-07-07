/*
 * null_update_bdy.hpp
 *
 *  Created on: June 25, 2020 
 *      Author: d.rosenberg
 */

#ifndef SRC_PDEINT_NULL_UPDATEBDY_HPP_
#define SRC_PDEINT_NULL_UPDATEBDY_HPP_


#include <memory>
#include <vector>
#include "update_bdy_base.hpp"

namespace geoflow {
namespace pdeint {

/**
 * Class to handle no-op for bdy updates with time
 *
 */
template<typename TypePack>
class NullUpdateBdy : public UpdateBdyBase<TypePack> {

public:
        using Types        = TypePack;
	using State        = typename Types::State;
	using StateInfo    = typename Types::StateInfo; // May contain time, time index, var name etc
	using Grid         = typename Types::Grid;
	using Ftype        = typename Types::Value;
        using Time         = typename Types::Time;
	using Size         = typename Types::Size;

      
	UpdateBdyBase() = delete;

	/**
	 * Constructor to initialize everything needed to do IO
	 *
	 */
	UpdateBdyBase();
	UpdateBdyBase(const UpdateBdyBase& I) = default;
	~UpdateBdyBase() = default;
	UpdateBdyBase& operator=(const UpdateBdyBase& I) = default;

	/**
	 * Update bdy conditions with state at t
	 *
	 * @param[in,out] ptree  Initial time at start, and final time
	 * @param[in,out] u  Current state values
	 */
	bool update (Grid      &grid, 
                     StateInfo &stinfo, 
                     Time      &time, 
                     State     &utmp, 
                     State     &u, 
                     State     &ub){
                        return this->update_impl(grid, stinfo, time, utmp, u, ub);
                     }

protected:
	bool update_impl (
                     Grid      &grid, 
                     StateInfo &stinfo, 
                     Time      &time, 
                     State     &utmp, 
                     State     &u, 
                     State     &ub) {
                // Do nothing ...
                return true;
              }

};


} // namespace pdeint
} // namespace geoflow



#endif /* SRC_PDEINT_NULL_UPDATEBDY_HPP_ */
