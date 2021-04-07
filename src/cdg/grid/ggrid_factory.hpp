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
#include "pdeint/io_base.hpp"
#include "pdeint/observer_base.hpp"
#include "ggrid_icos.hpp"
#include "ggrid_box.hpp"


//typedef GFTYPE                      Time;
//typedef GTVector<GTVector<GFTYPE>*> State;

template<typename TypePack>
class GGridFactory
{
  public:

        using Types       = TypePack;
        using State       = typename Types::State;
        using StateInfo   = typename Types::StateInfo;
        using Grid        = typename Types::Grid;
        using Time        = typename Types::Time;
        using IOBaseType  = IOBase<Types>;
        using IOBasePtr   = std::shared_ptr<IOBaseType>;
        using ObsTraits   = typename ObserverBase<Types>::Traits;
        using BdyUpdatePtr= std::vector<std::shared_ptr<UpdateBdyBase<Types>>>;



	static GGrid<Types>  *build(const geoflow::tbox::PropertyTree& ptree, GTVector<GNBasis<GCTYPE,GFTYPE>*> gbasis, IOBasePtr pIO, ObsTraits &obstraits, GC_COMM &comm);


  private:
        static void   read_grid(const geoflow::tbox::PropertyTree& ptree, GTMatrix<GINT> &p, GTVector<GTVector<GFTYPE>> &xnodes, IOBasePtr pIO, ObsTraits &obstraits, GC_COMM &comm);


}; // class GGridFactory

#include "ggrid_factory.ipp"

#endif 
