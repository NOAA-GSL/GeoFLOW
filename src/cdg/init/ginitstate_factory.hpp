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
#include "gcomm.hpp"
#include "gtvector.hpp"
#include "ggrid.hpp"
#include "ginitstate_user.hpp"

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


	static GBOOL init(const geoflow::tbox::PropertyTree& ptree, GGrid &grid, EqnBasePtr &peqn,  Time &time, State &utmp, State &ub, State &u);

  private:
	GBOOL set_by_direct(const PropertyTree& ptree, GGrid &grid, EqnBasePtr &peqn,  Time &time, State &utmp, State &ub, State &u);
	GBOOL set_by_comp  (const PropertyTree& ptree, GGrid &grid, EqnBasePtr &peqn,  Time &time, State &utmp, State &ub, State &u);
}; // class GInitSFactory
        GBOOL doinitv      (const PropertyTree &vtree, GGrid &grid, EqnBasePtr &peqn,  Time &time, State &utmp, State &ub, State &u);
        GBOOL doinitb      (const PropertyTree &vtree, GGrid &grid, EqnBasePtr &peqn,  Time &time, State &utmp, State &ub, State &u);
        GBOOL doinits      (const PropertyTree &vtree, GGrid &grid, EqnBasePtr &peqn,  Time &time, State &utmp, State &ub, State &u);
        GBOOL doinitps     (const PropertyTree &vtree, GGrid &grid, EqnBasePtr &peqn,  Time &time, State &utmp, State &ub, State &u);
        GBOOL doinitc      (const PropertyTree &vtree, GGrid &grid, EqnBasePtr &peqn,  Time &time, State &utmp, State &ub, State &u);



#include "ginitstate_factory.ipp"

#endif 
