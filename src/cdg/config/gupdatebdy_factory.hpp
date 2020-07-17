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
#include "gdirichlet_bdy.hpp"
#include "ginflow_bdy.hpp"
#include "gnoslip_bdy.hpp"
#include "gsimple_outflow_bdy.hpp"
#include "gsponge_bdy.hpp"
#include "gns_inflow_user.hpp"
#include "gutils.hpp"
#include "gspecbdy_user.hpp"

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
        using UpdateBase    = UpdateBdyBase<Types>;
        using UpdateBasePtr = std::shared_prt<UpdateBase>;
        using UpdateBaseList= std::vector<UpdateBasePtr>;
        using CallbackPtr   = std::function<void(
                                Grid       &grid,
                                StateInfo  &stinfo,
                                Time       &time,
                                State      &utmp,
                                State      &u,
                                State      &ub)>;


	static UpdateBdyBasePtr build(const PropertyTree& sptree, GString &supdate, Grid &grid, const GINT id, GBdyType bdytype, GTVector<GINT> &istate, GTVector<GSIZET> &ibdy);

	static UpdateBdyBasePtr  get_bdy_class (const PropertyTree& sptree, Grid &grid, const GINT id, const GBdyType bdytype, GTVector<GINT> &istate, GTVector<GSIZET> &ibdy);

  private:
               CallbackPtr       get_inflow_callback(const GString& sname, const GINT id);

}; // class GUpdateBdyFactory


#include "gupdatebdy_factory.ipp"

#endif 
