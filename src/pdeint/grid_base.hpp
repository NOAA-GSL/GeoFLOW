/*
 * grid_base.hpp
 *
 *  Created on: April 5, 2021 
 *      Author: d.rosenberg
 */

#ifndef SRC_PDEINT_GRIDBASE_HPP_
#define SRC_PDEINT_GRIDBASE_HPP_


#include <memory>
#include <vector>
#include "equation_base.hpp"
#include "update_bdy_base.hpp"

namespace geoflow {
namespace pdeint {

/**
 * Base to use for abstracting grids
 *
 */
template<typename TypePack>
class GridBase {

public:

        using Types          = TypePack;
        using State          = typename Types::State;
        using Grid           = typename Types::Grid;
        using Ftype          = typename Types::Ftype;
        using Time           = typename Types::Time;
        using EqnBase        = typename Types::EqnBase;
        using EqnBasePtr     = typename Types::EqnBasePtr;
        using UpdateBase     = UpdateBdyBase<Types>;
        using UpdateBasePtr  = std::shared_ptr<UpdateBase>;
        using BdyUpdateList  = std::vector<std::vector<UpdateBasePtr>>;

        // NOTE: BdyUpdateList is a list of boundary condition updates
        //       for each canonical or other boundary

      
	GridBase() = default;
	GridBase(const GridBase& I) = default;
	~GridBase() = default;
	GridBase& operator=(const GridBase& I) = default;


	/**
	 * Return bdy condition/update list 
	 *
	 */
	BdyUpdateList& bdy_update_list () {
                       return bdylist_;
                     }

protected:

                     
         BdyUpdateList bdylist_;

};


} // namespace pdeint
} // namespace geoflow



#endif /* SRC_PDEINT_GRIDBASE_HPP_ */
