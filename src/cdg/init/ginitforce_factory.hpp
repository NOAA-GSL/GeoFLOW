//==================================================================================
// Module       : ginitforce_factory.hpp
// Date         : 7/11/19 (DLR)
// Description  : GeoFLOW forcing variable initialization factory object. 
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GINITFORCE_FACTORY_HPP)
#define _GINITFORCE_FACTORY_HPP 

#include "tbox/property_tree.hpp"
#include "pdeint/equation_base.hpp"
#include "gcomm.hpp"
#include "gtvector.hpp"
#include "ggrid.hpp"
#include "ginitforce_direct_user.hpp"
#include "ginitstate_direct_user.hpp"
#include "ginitforce_comp.h"

using namespace geoflow::pdeint;
using namespace std;

template<typename EquationType>
class GInitForceFactory
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
	GBOOL set_by_comp (const PropertyTree& ptree, GGrid &grid, EqnBasePtr &peqn,  Time &time, State &utmp, State &ub, State &u);

        GBOOL doinitfv (const PropertyTree &vtree, GGrid &grid, Time &time, State &utmp, State &ub, State &u);
        GBOOL doinitfb (const PropertyTree &vtree, GGrid &grid, Time &time, State &utmp, State &ub, State &u);
        GBOOL doinitfs (const PropertyTree &vtree, GGrid &grid, Time &time, State &utmp, State &ub, State &u);
        GBOOL doinitfps(const PropertyTree &vtree, GGrid &grid, Time &time, State &utmp, State &ub, State &u);

}; // end, class GInitForceFactory


#include "ginitforce_factory.ipp"

#endif 
