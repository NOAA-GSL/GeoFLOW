//==================================================================================
// Module       : gupdatebdy_factory.hpp
// Date         : 7/11/19 (DLR)
// Description  : GeoFLOW bdy update object factory. 
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GUPDATEBDY_FACTORY_HPP)
#define _GUPDATEBDY_FACTORY_HPP 

#include "tbox/property_tree.hpp"
#include "gcomm.hpp"
#include "gtvector.hpp"
#include "ggrid.hpp"
#include "gfrominit_bdy.hpp"
#include "gsimple_outflow_bdy.hpp"
#include "gsponge_bdy.hpp"

using namespace geoflow;
using namespace geoflow::tbox;
using namespace std;

template<typename TypePack>
class GUpdateBdyFactory
{
  public:
        using Types         = TypePack;
        using State         = typename Equation::State;
        using StateInfo     = typename Equation::StateInfo;
        using Grid          = typename Equation::Grid;
        using Value         = typename Equation::Value;
        using Time          = typename Equation::Time;
        using UpdateBase    = UpdateBdyBase<Types>
        using UpdateBasePtr = std::shared_prt<Types>


	static UpdateBdyBasePtr build(const PropertyTree& ptree, Grid &grid, StateInfo &stinfo);

  private:
        static void  set_bdy_from_state(const PropertyTree& ptree, GString &sconfig, Grid &grid, StateInfo &stinfo, Time &time, State &utmp, State &u, State &ub);

}; // class GUpdateBdyFactory


#include "gupdatebdy_factory.ipp"

#endif 