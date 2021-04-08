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
#include "ginitforce_direct_user.hpp"
#include "ginitstate_direct_user.hpp"
#include "ginitforce_comp.h"

using namespace geoflow::pdeint;
using namespace std;

template<typename TypePack>
class GInitForceFactory
{
  public:
        using Types         = TypePack;
        using EqnBase       = EquationBase<TypePack>;
        using EqnBasePtr    = std::shared_ptr<EqnBase>;
        using State         = typename Types::State;
        using Grid          = typename Types::Grid;
        using GridBox       = typename Types::GridBox;
        using GridIcos      = typename Types::GridIcos;
        using CompDesc      = typename Types::CompDesc;
        using Ftype         = typename Types::Ftype;
        using Time          = typename Types::Time;


	static GBOOL init(const geoflow::tbox::PropertyTree& ptree, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u, State &uf);

  private:
	static GBOOL set_by_direct(const PropertyTree& ptree, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u, State &uf);
	static GBOOL set_by_comp  (const PropertyTree& ptree, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u, State &uf);

        static GBOOL doinitfv     (const PropertyTree &ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u, State &uf);
        static GBOOL doinitftemp  (const PropertyTree &ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u, State &uf);

}; // end, class GInitForceFactory


#include "ginitforce_factory.ipp"

#endif 
