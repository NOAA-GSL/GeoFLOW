//==================================================================================
// Module       : ggrid_factory
// Date         : 2/1/19 (DLR)
// Description  : GeoFLOW grid factory object. 
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#if !defined(_GGRID_FACTORY_HPP)
#define _GGRID_FACTORY_HPP 

#include "gcomm.hpp"
#include "gtvector.hpp"
#include "gnbasis.hpp"
#include "ggrid_box.hpp"
#include "ggrid_icos.hpp"
#include "pdeint/update_bdy_base.hpp"
#include "pdeint/io_base.hpp"
#include "pdeint/observer_base.hpp"
#include "pdeint/equation_base.hpp"
#include "tbox/property_tree.hpp"
#include "tbox/tracer.hpp"



template<typename TypePack>
class GGridFactory
{
  public:

        using Types       = EquationBase<TypePack>;;
        using State       = typename Types::State;
        using StateInfo   = typename Types::StateInfo;
        using Time        = typename Types::Time;
        using Ftype       = typename Types::Ftype;
        using IOBaseType  = IOBase<Types>;
        using IOBasePtr   = std::shared_ptr<IOBaseType>;
        using ObsTraits   = typename ObserverBase<Types>::Traits;
        using BdyUpdatePtr= std::vector<std::shared_ptr<UpdateBdyBase<Types>>>;



	static GGrid  *build(const geoflow::tbox::PropertyTree& ptree, GTVector<GNBasis<GCTYPE,GFTYPE>*> gbasis, IOBasePtr pIO, ObsTraits &obstraits, GC_COMM &comm);


        static void   read_grid(const geoflow::tbox::PropertyTree& ptree, GTMatrix<GINT> &p, GTVector<GTVector<GFTYPE>> &xnodes, IOBasePtr pIO, ObsTraits &obstraits, GC_COMM &comm);


}; // class GGridFactory


#include "ggrid_factory.ipp"


#endif 
