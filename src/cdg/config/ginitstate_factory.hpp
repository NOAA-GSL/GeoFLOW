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
#include "ginitstate_direct_user.hpp"
#include "ginitstate_comp.h"
#include "ggrid.hpp"
#include "ggrid_icos.hpp"
#include "ggrid_box.hpp"
#include "ginitb.hpp"
#include "ginitc.hpp"
#include "ginitps.hpp"
#include "ginits.hpp"
#include "ginitv.hpp"

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
        using StateInfo     = typename Equation::StateInfo;
        using Grid          = typename Equation::Grid;
        using CompDesc      = typename Equation::CompDesc;
        using Value         = typename Equation::Value;
        using Time          = typename Equation::Time;


	static GBOOL init(const geoflow::tbox::PropertyTree& ptree, Grid &grid, StateInfo &stinfo, Time &time, State &utmp, State &ub, State &u);

  private:
	static GBOOL set_by_direct(const PropertyTree& ptree, Grid &grid, StateInfo &stinfo,  Time &time, State &utmp, State &ub, State &u);
	static GBOOL set_by_comp  (const PropertyTree& ptree, Grid &grid, StateInfo &stinfo,  Time &time, State &utmp, State &ub, State &u);

        static GBOOL doinitv      (const PropertyTree &ptree, GString &sconfig,  Grid &grid, StateInfo &stinfo,  Time &time, State &utmp, State &ub, State &u);
        static GBOOL doinitb      (const PropertyTree &ptree, GString &sconfig,  Grid &grid, StateInfo &stinfo, Time &time, State &utmp, State &ub, State &u);
        static GBOOL doinitdt     (const PropertyTree &ptree, GString &sconfig,  Grid &grid, StateInfo &stinfo, Time &time, State &utmp, State &ub, State &u);
        static GBOOL doinitmfrac  (const PropertyTree &ptree, GString &sconfig,  Grid &grid, StateInfo &stinfo, Time &time, State &utmp, State &ub, State &u);
        static GBOOL doinitenergy (const PropertyTree &ptree, GString &sconfig,  Grid &grid, StateInfo &stinfo, Time &time, State &utmp, State &ub, State &u);
        static GBOOL doinittemp   (const PropertyTree &ptree, GString &sconfig,  Grid &grid, StateInfo &stinfo, Time &time, State &utmp, State &ub, State &u);
        static GBOOL doinitc      (const PropertyTree &ptree, GString &sconfig,  Grid &grid, StateInfo &stinfo,  Time &time, State &utmp, State &ub, State &u);

}; // end, class GInitStateFactory

#include "ginitstate_factory.ipp"

#endif 
