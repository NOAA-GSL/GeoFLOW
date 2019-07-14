//==================================================================================
// Module       : ginit_factory.hpp
// Date         : 7/11/19 (DLR)
// Description  : GeoFLOW state initialization factory object. 
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GINIT_FACTORY_HPP)
#define _GINIT_FACTORY_HPP 

#include "tbox/property_tree.hpp"
#include "gcomm.hpp"
#include "gtvector.hpp"
#include "ggrid.hpp"

namespace geoflow {
namespace pdeint {


class GInitFactory
{
  public:
        using Equation      = EquationType;
        using EqnBase       = EquationBase<EquationType>;
        using EqnBasePtr    = std::shared_ptr<EqnBase>;
        using State         = typename Equation::State;
        using Grid          = typename Equation::Grid;
        using Value         = typename Equation::Value;
        using Time          = typename Equation::Time;


	static void init(const geoflow::tbox::PropertyTree& ptree, EqnBasePtr &eqn_ptr,  Time &time, State *ub, State &u);

  private:
}; // class GInitFactory




} // namespace pdeint
} // namespace geoflow


#include "ginit_factory.ipp"

#endif 
