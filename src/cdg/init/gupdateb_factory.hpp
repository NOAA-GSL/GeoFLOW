//==================================================================================
// Module       : gupdateb_factory.hpp
// Date         : 7/11/19 (DLR)
// Description  : GeoFLOW boundary update factory object. 
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GUPDATEB_FACTORY_HPP)
#define _GUPDATEB_FACTORY_HPP 

#include "tbox/property_tree.hpp"
#include "gcomm.hpp"
#include "gtvector.hpp"
#include "ggrid.hpp"
#include "gupdateb.hpp"

namespace geoflow {
namespace pdeint {


template<typename EquationType>
class GUpdateBFactory
{
  public:
        using Equation      = EquationType;
        using EqnBase       = EquationBase<EquationType>;
        using EqnBasePtr    = std::shared_ptr<EqnBase>;
        using State         = typename Equation::State;
        using Grid          = typename Equation::Grid;
        using Value         = typename Equation::Value;
        using Time          = typename Equation::Time;


	static void update(const geoflow::tbox::PropertyTree& ptree, GGrid &grid,  Time &time, State &utmp, State &ub, State &u);

  private:
}; // class GUpdateBFactory




} // namespace pdeint
} // namespace geoflow


#include "ginits_factory.ipp"

#endif 
