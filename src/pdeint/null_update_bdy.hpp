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
        using EqnBase      = EquationBase<Types>;
        using EqnBasePtr   = std::shared_ptr<EqnBase>;
	using Grid         = typename Types::Grid;
	using Ftype        = typename Types::Ftype;
        using Time         = typename Types::Time;

      
	NullUpdateBdy() = default;
	NullUpdateBdy(const NullUpdateBdy& I) = default;
	~NullUpdateBdy() = default;
	NullUpdateBdy& operator=(const NullUpdateBdy& I) = default;

	/**
	 * Update bdy conditions with state at t
	 *
	 * @param[in,out] ptree  Initial time at start, and final time
	 * @param[in,out] u  Current state values
	 */
	bool update (EqnBasePtr &eqn,
                     Grid       &grid, 
                     Time       &time, 
                     State      &utmp, 
                     State      &u){
                        return this->update_impl(eqn, grid, time, utmp, u);
                     }

protected:
	bool update_impl (
                     EqnBasePtr &eqn,
                     Grid       &grid, 
                     Time       &time, 
                     State      &utmp, 
                     State      &u) {
                // Do nothing ...
                return true;
              }

        std::vector<int>&
             get_istate_impl(){ return idummy_; }

private:
        std::vector<int> idummy_;

};


} // namespace pdeint
} // namespace geoflow



#endif /* SRC_PDEINT_NULL_UPDATEBDY_HPP_ */
