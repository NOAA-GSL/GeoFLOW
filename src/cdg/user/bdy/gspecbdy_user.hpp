//==================================================================================
// Module       : gspecbdy_user.hpp
// Date         : 7/11/19 (DLR)
// Description  : User-specified boundary specification function implementations 
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GSPECBDY_USER_HPP)
#define _GSPECBDY_USER_HPP

#include "gtypes.h"
#include "gtvector.hpp"
#include "gutils.hpp"
#include "tbox/property_tree.hpp"
#include "pdeint/equation_base.hpp"


using namespace geoflow;
using namespace geoflow::tbox;



template<typename TypePack>
struct gspecbdy
{
        using Types         = TypePack;
        using EqnBase       = EquationBase<TypePack>;
        using EqnBasePtr    = std::shared_ptr<EqnBase>;
        using State         = typename Types::State;
        using Grid          = typename Types::Grid;
        using Ftype         = typename Types::Ftype;
        using Time          = typename Types::Ftype;



        static GBOOL impl_my_mixed_bdy(const PropertyTree &stree, Grid &grid,  const GINT id, GTVector<GSIZET> &ibdy);

};


#include "gspecbdy_user.ipp"


#endif
