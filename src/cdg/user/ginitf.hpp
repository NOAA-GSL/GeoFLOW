//==================================================================================
// Module       : ginitf.hpp
// Date         : 7/10/19 (DLR)
// Description  : Force initialization function implementations
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GINITF_HPP)
#define _GINITF_HPP

#include "tbox/property_tree.hpp"
#include "gtypes.h"
#include "gtvector.hpp"
#include "ggrid.hpp"
#include "gcomm.hpp"


namespace ginitf
{

GBOOL initf_impl_null    (const PropertyTree &ftree, const Time &t, State &u, State &ub);
GBOOL initf_impl_rand    (const PropertyTree &ftree, const Time &t, State &u, State &ub);

};

#include "ginitf_null.ipp"
#include "ginitf_rand.ipp"

#endif
