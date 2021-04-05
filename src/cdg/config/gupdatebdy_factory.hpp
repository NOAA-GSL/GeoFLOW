//==================================================================================
// Module       : gupdatebdy_factory.hpp
// Date         : 7/11/19 (DLR)
// Description  : GeoFLOW bdy update object factory. 
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GUPDATEBDY_FACTORY_HPP)
#define _GUPDATEBDY_FACTORY_HPP 

#include "gcomm.hpp"
#include "gtvector.hpp"
#include "ggrid.hpp"
#include "ggrid_box.hpp"
#include "ggrid_icos.hpp"
#include "g0flux_bdy.hpp"
#include "gdirichlet_bdy.hpp"
#include "ginflow_bdy.hpp"
#include "gnoslip_bdy.hpp"
#include "goutflow_bdy.hpp"
#include "gsponge_bdy.hpp"
#include "gns_inflow_user.hpp"
#include "gutils.hpp"
#include "gspecbdy_user.hpp"
#include "ginitstate_factory.hpp"
#include "pdeint/null_update_bdy.hpp"
#include "tbox/property_tree.hpp"

using namespace geoflow;
using namespace geoflow::tbox;
using namespace std;

struct stBdyBlock {
  GBOOL            use_init=FALSE;    // use initialisation method for setting bdy conditions
  GBOOL            compute_once=FALSE;// compute data once (not time-dep)
  GINT             idir=-1;           // coord direction of bdy surface
  GINT             bdyid=-1;          // bdy id 
  GBdyType         tbdy=GBDY_NONE;    // bdy type
  GFTYPE           xstart=0.0;        // start of sponge surf in idir direction
  GFTYPE           xmax=0.0;          // max coord in direction idir
  vector<GINT>     istate;            // vector of state indices
  vector<GFTYPE>   value;             // vector of Dirichlet values for each istate
  vector<GFTYPE>   farfield;          // vector of far field bdys for SPONGE bcs for each istate
  vector<GFTYPE>   falloff;           // vector of fall-off rates for SPONGE bcs for each istate
  vector<GFTYPE>   diffusion;         // vector of diffusion factors for SPONGE bcs for each istate
  GString          bdyclass;          // bdy ('uniform', 'mixed')
  GString          smethod;           // name of method providing bdy values (e.g. for inflow)
};


template<typename TypePack>
class GUpdateBdyFactory
{
  public:
        using Types            = TypePack;
        using EqnBase          = EquationBase<TypePack>;
        using EqnBasePtr       = std::shared_ptr<EqnBase>;
        using State            = typename Types::State;
        using Grid             = typename Types::Grid;
        using Ftype            = typename Types::Ftype;
        using Time             = typename Types::Time;
        using UpdateBdyBasePtr = shared_ptr<UpdateBdyBase<Types>>;
        using CallbackPtr      = std::function<GBOOL(
                                EqnBasePtr &eqn,
                                Grid       &grid,
                                Time       &time,
                                const GINT  id,
                                State      &utmp,
                                State      &u)>;


	static UpdateBdyBasePtr build(const PropertyTree &ptree, const GString &sbdy, Grid &grid, stBdyBlock &bcblock, GTVector<GSIZET> &ibdy, GSIZET igbdy_start);

	static UpdateBdyBasePtr  get_bdy_class (const PropertyTree& ptree, Grid &grid, stBdyBlock &bcblock,  GTVector<GSIZET> &ibdy, GSIZET igbdy_start);

        static GBOOL             get_bdy_block(const geoflow::tbox::PropertyTree &sptree, GString &sbdy, GINT ibc, stBdyBlock &stblock);

        static GINT              bdy_block_conform_per(const geoflow::tbox::PropertyTree &sptree);

  private:
        static  CallbackPtr       get_inflow_callback(const GString& sname, const GINT id);

}; // class GUpdateBdyFactory


#include "gupdatebdy_factory.ipp"

#endif 
