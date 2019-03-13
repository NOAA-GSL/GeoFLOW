//==================================================================================
// Module       : gio.h
// Date         : 3/2/19 (DLR)
// Description  : Simple GeoFLOW IO
// Copyright    : Copyright 2019. Colorado State University. All rights reserved.
// Derived From : 
//==================================================================================

#include "gtypes.h"
#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include "gcomm.hpp"
#include "ggrid.hpp"
#include "gtvector.hpp"
#include "tbox/property_tree.hpp"
#include "tbox/mpixx.hpp"

using namespace geoflow::tbox;
using namespace std;

typedef GTVector<GTVector<GFTYPE>*> State;
typedef GTVector<GFTYPE> StateElem;

int gio(GGrid &grid, State &u, GINT nu, GSIZET tindex, GTVector<GString> &svars, GC_COMM comm, GBOOL &bgridprinted);
