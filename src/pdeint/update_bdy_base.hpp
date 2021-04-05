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
        using Types        = TypePack;
	using State        = typename Types::State;
        using EqnBase      = EquationBase<TypePack>;
        using EqnBasePtr   = std::shared_ptr<EqnBase>;
	using Grid         = typename Types::Grid;
	using Ftype        = typename Types::Ftype;
        using Time         = typename Types::Time;
	using IBdyVol      = typename Types::IBdyVol;
	using TBdyVol      = typename Types::TBdyVol;

      
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
	 * @param[in,out] eqn    : equation pointer
	 * @param[in,out] grid   : initial time at start, and final time
	 * @param[in,out] time   : current time
	 * @param[in,out] utmp   : tmp arrays
	 * @param[in,out] u      : current state array
	 */
	bool update (EqnBasePtr &eqn
                     Grid       &grid, 
                     Time       &time, 
                     State      &utmp, 
                     State      &u, ){ 
                        return this->update_impl(eqn, grid, stinfo, time, utmp, u);
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
