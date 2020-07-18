//==================================================================================
// Module       : gutils.hpp
// Date         : 1/31/19 (DLR)
// Description  : GeoFLOW utilities namespace
// Copyright    : Copyright 2019. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GUTILS_HPP)
#define _GUTILS_HPP

#include <cassert>
#include "gtypes.h"
#include "ggfx.hpp"
#include "gmass.hpp"
#include "gmtk.hpp"
//#include "ggrid_icos.hpp"
//#include "ggrid_box.hpp"
#include "tbox/property_tree.hpp"

using namespace geoflow::pdeint;
using namespace std;

class GGrid;
class GGridBox;
class GGridIcos;

struct stBdyBlock {
  GTVector<GTVector<GINT>>  istate; // vector of staet index vectors
  GTVector<GBdyType>        tbdy;   // vector of bdy types; one for each istate vector
  GString                   config_method; // name of bdy node config method
  GString                   inflow_mehod;  // name of inflow method, if any
};

namespace geoflow
{

GBdyType       str2bdytype (const GString &stype);
GStateCompType str2comptype(const GString &stype);
GBOOL          file_empty(GString filename);
GBOOL          get_bdy_block(const geoflow::tbox::PropertyTree &sptree, stBdyBlock &stblock);
template<typename T>
void           append(GTVector<T> &base, GTVector<T> &add);
template<typename T>
void           unique(GSIZET ibeg, GSIZET iend, GTVector<GSIZET> &iunique);
template<typename T>
void           smooth(GGrid &grid, GGFX_OP op,  GTVector<T> &tmp, GTVector<T> &v);
template<typename T>
void           coord_dims(const geoflow::tbox::PropertyTree& ptree, GTVector<T> &xmin, GTVector<T> &xmax);


} // end, namespace

#include "gutils.ipp"

#endif
