/*
 * io_base.hpp
 *
 *  Created on: June 25, 2020 
 *      Author: d.rosenberg
 */

#ifndef SRC_PDEINT_UPDATEBDYBASE_HPP_
#define SRC_PDEINT_UPDATEBDYBASE_HPP_


#include <memory>
#include <vector>

namespace geoflow {
namespace pdeint {

/**
 * Base to handle time updates of bdy information 
 * (or possibly state due to bdy) 
 *
 */
template<typename TypePack>
class UpdateBdyBase {

public:
        enum GIOType        {GIO_POSIX=0, GIO_COLL}; // POSIX or collective
        using Types        = TypePack;
	using State        = typename Types::State;
	using StateInfo    = typename Types::StateInfo; // May contain time, time index, var name etc
	using Grid         = typename Types::Grid;
	using Ftype        = typename Types::Ftype;
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
                     State     &ub) = 0;

};


} // namespace pdeint
} // namespace geoflow



#endif /* SRC_PDEINT_UPDATEBDYBASE_HPP_ */
