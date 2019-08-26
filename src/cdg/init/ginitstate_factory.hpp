//==================================================================================
// Module       : ginitstate_factory.hpp
// Date         : 7/11/19 (DLR)
// Description  : GeoFLOW state variable initialization factory object. 
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GINITSTATE_FACTORY_HPP)
#define _GINITSTATE_FACTORY_HPP 

#include "tbox/property_tree.hpp"
#include "pdeint/equation_base.hpp"
#include "gcomm.hpp"
#include "gtvector.hpp"
#include "ggrid.hpp"
#include "ginitstate_direct_user.hpp"
#include "ginitstate_comp.h"

using namespace geoflow::pdeint;
using namespace std;

template<typename EquationType>
class GInitStateFactory
{
  public:
        using Equation      = EquationType;
        using EqnBase       = EquationBase<EquationType>;
        using EqnBasePtr    = std::shared_ptr<EqnBase>;
        using State         = typename Equation::State;
        using Grid          = typename Equation::Grid;
        using CompDesc      = typename Equation::CompDesc;
        using Value         = typename Equation::Value;
        using Time          = typename Equation::Time;


	static GBOOL init(const geoflow::tbox::PropertyTree& ptree, GGrid &grid, EqnBasePtr &peqn, Time &time, State &utmp, State &ub, State &u);

  private:
	static GBOOL set_by_direct(const PropertyTree& ptree, GGrid &grid, EqnBasePtr &peqn,  Time &time, State &utmp, State &ub, State &u);
	static GBOOL set_by_comp  (const PropertyTree& ptree, GGrid &grid, EqnBasePtr &peqn,  Time &time, State &utmp, State &ub, State &u);

        static GBOOL doinitv      (const PropertyTree &vtree, GGrid &grid, Time &time, State &utmp, State &ub, State &u);
        static GBOOL doinitb      (const PropertyTree &vtree, GGrid &grid, Time &time, State &utmp, State &ub, State &u);
        static GBOOL doinits      (const PropertyTree &vtree, GGrid &grid, Time &time, State &utmp, State &ub, State &u);
        static GBOOL doinitps     (const PropertyTree &vtree, GGrid &grid, Time &time, State &utmp, State &ub, State &u);
        static GBOOL doinitc      (const PropertyTree &vtree, GGrid &grid, Time &time, State &utmp, State &ub, State &u);

}; // end, class GInitStateFactory

#include "ginitstate_factory.ipp"

#endif 
