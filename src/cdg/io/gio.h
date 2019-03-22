//==================================================================================
// Module       : gio.h
// Date         : 3/2/19 (DLR)
// Description  : Simple GeoFLOW IO
// Copyright    : Copyright 2019. Colorado State University. All rights reserved.
// Derived From : 
//==================================================================================
#if !defined(_GIO_H)

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


void gio(GGrid &grid, const State &u, const GTVector<GINT> &nu, const GSIZET tindex, const GFTYPE time, const GTVector<GString> &svars, const GString &sdir, GC_COMM comm, const GBOOL bprgrid);


#endif
