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
#include "gtpoint.hpp"
#include "gmtk.hpp"
#include "tbox/property_tree.hpp"

using namespace std;



namespace geoflow
{

GBdyType       str2bdytype (const GString &stype);
GStateCompType str2comptype(const GString &stype);
GBOOL          file_empty(GString filename);

template<typename T>
void           append(GTVector<T> &base, GTVector<T> &add);
template<typename T>
void           unique(GSIZET ibeg, GSIZET iend, GTVector<GSIZET> &iunique);
template<typename T>
void           coord_dims(const geoflow::tbox::PropertyTree& ptree, GTVector<T> &xmin, GTVector<T> &xmax);

template<typename T>
void compute_temp(const GTVector<T> &e, const GTVector<T> &d, const GTVector<T> &cv,  GTVector<T> &temp);
template<typename T>
void compute_p(const GTVector<T> &Temp, const GTVector<T> &d, const GTVector<T> &q, GFTYPE R, GTVector<T> &p);
template<typename T>
void compute_p(const GTVector<T> &e, const GTVector<T> &q, GFTYPE R, const GTVector<T> &cv, GTVector<T> &p);

template<typename T>
GSIZET in_seg(const GTVector<GTPoint<T>> &vertices, const GTVector<GTVector<T>> &x, T eps, GTVector<GSIZET> &ind);
template<typename T>
GSIZET in_seg(const GTVector<GTPoint<T>> &vertices, const GTVector<GTVector<T>> &x, const GTVector<GSIZET> &ix, T eps, GTVector<GSIZET> &ind);
template<typename T>
GSIZET in_poly(const GTVector<GTPoint<T>> &vertices, const GTVector<GTVector<T>> &x, T eps, GTVector<GSIZET> &ind);
template<typename T>
GSIZET in_poly(const GTVector<GTPoint<T>> &vertices, const GTVector<GTVector<T>> &x, const GTVector<GSIZET> &ix, T eps, GTVector<GSIZET> &ind);

template<typename T>
GSIZET fuzzyeq(const GTVector<T> &v, T vcomp, T eps, GTVector<GSIZET> &ind);
template<typename T>
GSIZET fuzzyeq(const GTVector<T> &v, const GTVector<GSIZET> &iv, T vcomp, T eps, GTVector<GSIZET> &ind);

template<typename TOLD, typename TNEW>
void convert(const GTVector<TOLD> &v, GTVector<TNEW> &vnew);

} // end, namespace geoflow

#include "gutils.ipp"

#endif
