//==================================================================================
// Module       : gupdateb_user.hpp
// Date         : 7/11/19 (DLR)
// Description  : Boundary update function implementations specified by
//                user
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GUPDATEB_USER_HPP)
#define _GUPDATEB_USER_HPP

#include "tbox/property_tree.hpp"
#include "gtypes.h"
#include "gtvector.hpp"
#include "ggrid.hpp"


using namespace geoflow;
using namespace geoflow::tbox;

namespace gupdateb
{

GBOOL impl_mybdyupdate          (const PropertyTree& stree, GGrid &grid,  Time &time, State &utmp, State &u, State &ub);

};

#include "gupdateb_user.ipp"

#endif
