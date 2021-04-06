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

using namespace geoflow::pdeint;
using namespace std;

template<typename TypePack>
class GInitStateFactory
{
  public:
        using Types         = TypePack;
        using EqnBase       = EquationBase<TypePack>;
        using EqnBasePtr    = std::shared_ptr<EqnBase>;
        using State         = typename Types::State;
        using Grid          = typename Types::Grid;
        using CompDesc      = typename Types::CompDesc;
        using Time          = typename Types::Time;


	static GBOOL init(const geoflow::tbox::PropertyTree& ptree, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u);

  private:
	static GBOOL set_by_direct(const PropertyTree& ptree, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u);
	static GBOOL set_by_comp  (const PropertyTree& ptree, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u);

        static GBOOL doinitv      (const PropertyTree &ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u);
        static GBOOL doinitdt     (const PropertyTree &ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u);
        static GBOOL doinitmfrac  (const PropertyTree &ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u);
        static GBOOL doinitenergy (const PropertyTree &ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u);
        static GBOOL doinittemp   (const PropertyTree &ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u);
        static GBOOL doinitc      (const PropertyTree &ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u);

}; // end, class GInitStateFactory

#include "ginitstate_factory.ipp"

#endif 
