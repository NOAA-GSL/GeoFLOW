//==================================================================================
// Module       : gbdy_factory.hpp
// Date         : 7/11/19 (DLR)
// Description  : GeoFLOW config and initialization factory object. 
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GBDY_FACTORY_HPP)
#define _GBDY_FACTORY_HPP 

#include "tbox/property_tree.hpp"
#include "gcomm.hpp"
#include "gtvector.hpp"
#include "ggrid.hpp"
#include "gbdy_impl.hpp"


namespace geoflow {
namespace pdeint {


class GBdyFactory
{
  public:
        using Equation      = EquationType;
        using EqnBase       = EquationBase<EquationType>;
        using EqnBasePtr    = std::shared_ptr<EqnBase>;
        using State         = typename Equation::State;
        using Grid          = typename Equation::Grid;
        using Value         = typename Equation::Value;
        using Time          = typename Equation::Time;


	static void init(const geoflow::tbox::PropertyTree& ptree, GGrid &grid,  Time &time, State &utmp, State &u, State &ub);
      
         GBOOL      is_time_dep() { return btime_dep_; }
         std::function<void(GGrid &grid, const Time &t, State &utmp,
                     State &u   , State &ub)>
                    get_update_callback() { return update_bdy_callback_; }
  private:
  GBOOL                                       btime_dep_;
  std::function<void(GGrid &grid, const Time &t, State &utmp, 
                     State &u   , State &ub)> 
                  update_bdy_callback_;
}; // class GBdyFactory




} // namespace pdeint
} // namespace geoflow


#include "gbdy_factory.ipp"

#endif 
