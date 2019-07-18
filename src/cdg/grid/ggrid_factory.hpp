//==================================================================================
// Module       : ggrid_factory
// Date         : 2/1/19 (DLR)
// Description  : GeoFLOW grid factory object. 
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#if !defined(_GGRID_FACTORY_HPP)
#define _GGRID_FACTORY_HPP 

#include "tbox/property_tree.hpp"
#include "gcomm.hpp"
#include "gtvector.hpp"
#include "ggrid.hpp"


class GGridFactory
{
  public:

	static GGrid *build(const geoflow::tbox::PropertyTree& ptree, GTVector<GNBasis<GCTYPE,GFTYPE>*> gbasis, GC_COMM &comm);

         GBOOL      is_time_dep_bdy() { return btime_dep_bdy_; }
         std::function<void(GGrid &grid, const Time &t, State &utmp,
                     State &u   , State &ub)>
                    get_bdy_update_callback() { return update_bdy_callback_; }


  private:
        static void   read_grid(const geoflow::tbox::PropertyTree& ptree, GC_COMM comm,  GTMatrix<GINT> &p, GTVector<GTVector<GFTYPE>> &xnodes);

         GBOOL   btime_dep_bdy_;
  std::function<void(GGrid &grid, const Time &t, State &utmp,
                     State &u   , State &ub)>
                  update_bdy_callback_;

}; // class GGridFactory

//#include "ggrid_factory.ipp"

#endif 
