//==================================================================================
// Module       : gspecbdy.hpp
// Date         : 7/11/19 (DLR)
// Description  : Boundary specification function implementations provided by
//                system
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GSPECBDY_HPP)
#define _GSPECBDY_HPP

#include "tbox/property_tree.hpp"
#include "gtypes.h"
#include "gtvector.hpp"
#include "ggrid.hpp"
#include "geoflow.hpp"


using namespace geoflow;
using namespace geoflow::tbox;

typedef GFTYPE                      Time;
typedef GTVector<GTVector<GFTYPE>*> State;


namespace gspecbdy
{

GBOOL impl_uniform        (const PropertyTree& stree, GGrid &grid,  const GINT id, GTVector<GSIZET> &ibdy, GTVector<GBdyType> &tbdy);

};

#include "gspecbdy.ipp"

#endif
