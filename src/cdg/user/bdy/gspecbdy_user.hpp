//==================================================================================
// Module       : gspecbdy_user.hpp
// Date         : 7/11/19 (DLR)
// Description  : Boundary specification function implementations specified by
//                user
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GSPECBDY_USER_HPP)
#define _GSPECBDY_USER_HPP

#include "tbox/property_tree.hpp"
#include "gtypes.h"
#include "gtvector.hpp"
#include "ggrid.hpp"


using namespace geoflow;
using namespace geoflow::tbox;

namespace gspecbdy
{

GBOOL impl_mybdyspec          (const PropertyTree& stree, GGrid &grid,  const GINT id, GTVector<GSIZET> &ibdy, GTVector<GBdyType> &tbdy);

};

#include "gspecbdy_user.ipp"

#endif
