//==================================================================================
// Module       : gtools.h
// Date         : 5/5/19 (DLR)
// Description  : GeoFLOW initialization tools
// Copyright    : Copyright 2019. Colorado State University. All rights reserved.
// Derived From : 
//==================================================================================
#if !defined(_GTOOLS_H)
#define _GTOOLS_H

#include "gtypes.h"
#include <cstdio>
#include <cmath>
#include <unistd.h>
#include <iostream>
#include <memory>
#include <cstdlib>
#include <cassert>
#include <random>
#include <typeinfo>
#include "gcomm.hpp"
#include "ggfx.hpp"
#include "gllbasis.hpp"
#include "gmorton_keygen.hpp"
#include "ggrid.hpp"
#include "ggrid_box.hpp"
#include "ggrid_icos.hpp"
#include "tbox/property_tree.hpp"

#if defined(GEOFLOW_USE_GPTL)
    #include "gptl.h"
#endif

using namespace geoflow::tbox;
using namespace std;

void init_ggfx(PropertyTree &ptree, GGrid &grid, GGFX<GFTYPE> &ggfx);


#endif
