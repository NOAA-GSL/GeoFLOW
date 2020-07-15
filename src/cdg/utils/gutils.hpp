//==================================================================================
// Module       : gutils.hpp
// Date         : 1/31/19 (DLR)
// Description  : GeoFLOW utilities namespace
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#if !defined(_GUTILS_HPP)
#define _GUTILS_HPP

#include <cassert>
#include "gtypes.h"
#include "ggrid.hpp"
#include "tbox/property_tree.hpp"

using namespace geoflow::pdeint;
using namespace std;


namespace geoflow
{

GBdyType       str2bdytype (const GString &stype);
GStateCompType str2comptype(const GString &stype);
GBOOL          file_empty(GString filename);

template<typename T>
void           unique(GSIZET ibeg, GSIZET iend, GTVector<GSIZET> &iunique);
template<typename T>
void           smooth(GGrid &grid, GGFX_OP op,  GTVector<T> &tmp, GTVector<T> &v);
template<typename T>
void           coord_dims(const PropertyTree& ptree, GTVector<T> &xmin, GTVector<T> &xmax);


} // end, namespace

#include "gutils.ipp"

#endif
