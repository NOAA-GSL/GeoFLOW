//==================================================================================
// Module       : gexec.h
// Date         : 7/1/18 (DLR)
// Description  : Header to be included in exececutive code only
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#if !defined(_GEXEC_HPP)
#define _GEXEC_HPP

#include <cstddef>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include "gtvector.hpp"
#include "gtmatrix.hpp"


#if 0   //!defined(_G_PRE_INSTANTIATION)
#define _G_PRE_INSTANTIATION
template class GTVector <GFLOAT>;
template class GTVector<GDOUBLE>;
template class GTVector  <GQUAD>;
template class GTVector   <GINT>;
template class GTVector  <GUINT>;
template class GTVector<GUSHORT>;
template class GTVector <GSHORT>;
template class GTVector  <GLONG>;
template class GTVector <GULONG>;
template class GTVector <GLLONG>;
template class GTVector <GSIZET>;

template class GTMatrix <GFLOAT>;
template class GTMatrix<GDOUBLE>;
template class GTMatrix  <GQUAD>;
template class GTMatrix   <GINT>;
template class GTMatrix  <GUINT>;
template class GTMatrix<GUSHORT>;
template class GTMatrix <GSHORT>;
template class GTMatrix  <GLONG>;
template class GTMatrix <GULONG>;
template class GTMatrix <GLLONG>;
template class GTMatrix <GSIZET>;
#endif

#endif // _GEXEC_HPP
