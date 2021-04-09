/*
 * update_bdy_base.hpp
 *
 *  Created on: June 25, 2020 
 *      Author: d.rosenberg
 */

#ifndef SRC_PDEINT_UPDATEBDYBASE_HPP_
#define SRC_PDEINT_UPDATEBDYBASE_HPP_


#include <memory>
#include <vector>
#include "equation_base.hpp"

namespace geoflow {
namespace pdeint {

/**
 * Base to handle specification and time updates of bdy information 
 * (or possibly state due to bdy) 
 *
 */
template<typename TypePack>
class UpdateBdyBase {

public:

        using Types      = TypePack;
        using State      = typename Types::State;
        using StateComp  = typename Types::StateComp;
        using StateInfo  = typename Types::StateInfo;
        using Grid       = typename Types::Grid;
        using Mass       = typename Types::Mass;
        using Ftype      = typename Types::Ftype;
        using Derivative = typename Types::Derivative;
        using Time       = typename Types::Time;
        using CompDesc   = typename Types::CompDesc;
        using Jacobian   = typename Types::Jacobian;
	using IBdyVol    = typename Types::IBdyVol;
	using TBdyVol    = typename Types::TBdyVol;
        using Size       = typename Types::Size;
        using EqnBase    = typename Types::EqnBase;
        using EqnBasePtr = typename Types::EqnBasePtr;
      
	UpdateBdyBase() = default;
	UpdateBdyBase(const UpdateBdyBase& I) = default;
	~UpdateBdyBase() = default;
	UpdateBdyBase& operator=(const UpdateBdyBase& I) = default;

#if 0
	/**
	 * spec bdy conditions 
	 *
	 * @param[in]     sptree : property tree block for specification
	 * @param[in]     grid   : grid object
	 * @param[in]     bdyid  : bdy id
	 * @param[in]     ibdy   : indirection array into volmume indicating bdy indices
	 * @param[in,out] tbdy   : bdy types for each ibdy element
	 */
	bool spec   (PropertyTree &sptree,
                     Grid         &grid, 
                     int           bdyid, 
                     IBdyVol      &ibdy, 
                     TBdyvol      &tbdy){
                        return this->spec_impl(sptree, grid, bdyid, ibdy, tbdy);
                     }
#endif


	/**
	 * Update bdy conditions with state at t
	 *
	 * @param[in,out] eqn    : equation base (not the shared pointer)
	 * @param[in,out] grid   : initial time at start, and final time
	 * @param[in,out] time   : current time
	 * @param[in,out] utmp   : tmp arrays
	 * @param[in,out] u      : current state array
	 */
	bool update (EqnBasePtr &eqn,
                     Grid       &grid, 
                     Time       &time, 
                     State      &utmp, 
                     State      &u){ 
                        return this->update_impl(eqn, grid, time, utmp, u);
                     }

protected:

#if 0
	virtual bool spec   (PropertyTree &sptree,
                     Grid         &grid, 
                     int           bdyid, 
                     IBdyVol      &ibdy, 
                     TBdyvol      &tbdy) { return true;}
#endif
                     
	virtual bool update_impl (
                     EqnBasePtr &eqn,
                     Grid       &grid, 
                     Time       &time, 
                     State      &utmp, 
                     State      &u) = 0; 

};


} // namespace pdeint
} // namespace geoflow



#endif /* SRC_PDEINT_UPDATEBDYBASE_HPP_ */
