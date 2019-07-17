//==================================================================================
// Module       : ginits_factory.hpp
// Date         : 7/11/19 (DLR)
// Description  : GeoFLOW state variable initialization factory object. 
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GINITS_FACTORY_HPP)
#define _GINITS_FACTORY_HPP 

#include "tbox/property_tree.hpp"
#include "gcomm.hpp"
#include "gtvector.hpp"
#include "ggrid.hpp"

namespace geoflow {
namespace pdeint {


class GInitSFactory
{
  public:
        using Equation      = EquationType;
        using EqnBase       = EquationBase<EquationType>;
        using EqnBasePtr    = std::shared_ptr<EqnBase>;
        using State         = typename Equation::State;
        using Grid          = typename Equation::Grid;
        using Value         = typename Equation::Value;
        using Time          = typename Equation::Time;


	static void init(const geoflow::tbox::PropertyTree& ptree, GGrid &grid,  Time &time, State &utmp, State &ub, State &u);

  private:
}; // class GInitSFactory




} // namespace pdeint
} // namespace geoflow


#include "ginit_factory.ipp"

#endif 