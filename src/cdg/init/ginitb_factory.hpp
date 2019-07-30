//==================================================================================
// Module       : gspecb_factory.hpp
// Date         : 7/11/19 (DLR)
// Description  : GeoFLOW boundary specification factory object. 
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GINITB_FACTORY_HPP)
#define _GINITB_FACTORY_HPP 

#include "tbox/property_tree.hpp"
#include "gcomm.hpp"
#include "gtvector.hpp"
#include "ggrid.hpp"
#include "gbdyspec_impl.hpp"

namespace geoflow {
namespace pdeint {


template<typename EquationType>
class GInitBFactory
{
  public:
        using Equation      = EquationType;
        using EqnBase       = EquationBase<EquationType>;
        using EqnBasePtr    = std::shared_ptr<EqnBase>;
        using State         = typename Equation::State;
        using Grid          = typename Equation::Grid;
        using Value         = typename Equation::Value;
        using Time          = typename Equation::Time;


	static GBOOL init   (const geoflow::tbox::PropertyTree& ptree, GGrid &grid,  Time &time, State &utmp, State &u, State &ub);
	static std::function<void(PropertyTree &ptree, GGrid &grid, const Time &t,
                           State &utmp, State &u   , State &ub)>
                     initfcn(const geoflow::tbox::PropertyTree& ptree, GGrid &grid,  Time &time, State &utmp, State &ub, State &u);

  private:
         GBOOL       doinits(const geoflow::tbox::PropertyTree& ptree, GGrid &grid, Time &time, State &utmp, State &u, State &ub);

}; // class GInitBFactory




} // namespace pdeint
} // namespace geoflow


#include "ginits_factory.ipp"

#endif 
