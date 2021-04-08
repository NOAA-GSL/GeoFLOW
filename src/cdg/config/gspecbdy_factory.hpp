//==================================================================================
// Module       : gspecbdy_factory.hpp
// Date         : 7/11/19 (DLR)
// Description  : GeoFLOW boundary specification factory object. 
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GSPECBDY_FACTORY_HPP)
#define _GSPECBDY_FACTORY_HPP 

#include "gtypes.h"
#include "gcomm.hpp"
#include "gtvector.hpp"
#include "ggrid.hpp"
#include "gspecbdy_user.hpp" // include user-provided spec functions
#include "tbox/property_tree.hpp"
#include "pdeint/equation_base.hpp"


using namespace geoflow;
using namespace geoflow::tbox;
using namespace std;

template<typename TypePack>
class GSpecBdyFactory
{
        using Types         = TypePack;
        using EqnBase       = EquationBase<TypePack>;
        using EqnBasePtr    = std::shared_ptr<EqnBase>;
        using State         = typename Types::State;
        using Grid          = typename Types::Grid;
        using Ftype         = typename Types::Ftype;
        using Time          = typename Types::Ftype;


  public:


        static GBOOL dospec(const geoflow::tbox::PropertyTree& ptree, Grid &grid, const GINT id,  GTVector<GSIZET> &ibdy);


}; // class GSpecBdyFactory


#include "gspecbdy_factory.ipp"


#endif 
