//==================================================================================
// Module       : ginitforce_user.hpp
// Date         : 7/10/19 (DLR)
// Description  : User force initialization function implementations
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GINITFORCE_USER_HPP)
#define _GINITFORCE_USER_HPP

#include "tbox/property_tree.hpp"
#include "gtypes.h"
#include "gtvector.hpp"
#include "ggrid.hpp"
#include "gcomm.hpp"


namespace ginitforce
{

GBOOL impl_null    (const PropertyTree &ftree, const Time &t, State &u, State &uf);
GBOOL impl_rand    (const PropertyTree &ftree, const Time &t, State &u, State &uf);

};

#include "ginitforce_null.ipp"
#include "ginitforce_rand.ipp"

#endif
