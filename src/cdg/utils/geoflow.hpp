//==================================================================================
// Module       : geoflow.hpp
// Date         : 1/31/19 (DLR)
// Description  : GeoFLOW utilities namespace
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#if !defined(_GEOFLOW_HPP)
#define _GEOFLOW_HPP

#include "gtypes.h"
#include <assert.h>


namespace geoflow
{

GBdyType str2bdytype(const GString &stype);

} // end, namespace

#endif
