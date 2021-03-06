//==================================================================================
// Module       : gterrainSpecSph.hpp
// Date         : 7/11/19 (DLR)
// Description  : Terrain specification function implementations provided
//                by user for spherical grids (2d and 3d).
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GTERRAINSPECSPH_USER_HPP)
#define _GTERRAINSPECSPH_USER_HPP

#include "tbox/property_tree.hpp"
#include "gtypes.h"
#include "gtvector.hpp"
#include "ggrid.hpp"


using namespace geoflow;
using namespace geoflow::tbox;

template<typename TypePack>
struct gterrainSpecSph
{       
        using Types      = TypePack;
        using State      = typename Types::State;
        using StateComp  = typename Types::StateComp;
        using EqnBase    = EquationBase<Types>;
        using EqnBasePtr = std::shared_ptr<EqnBase>;
        using Grid       = typename Types::Grid;
        using Time       = typename Types::Time;
        using Ftype      = typename Types::Ftype;


static GBOOL impl_gauss_range     (const PropertyTree& ptree, GString sblk, Grid &grid,  State &utmp, State &xb);
static GBOOL impl_poly_range      (const PropertyTree& ptree, GString sblk, Grid &grid,  State &utmp, State &xb);

};


#include "gterrain_specsph_user.ipp"


#endif
