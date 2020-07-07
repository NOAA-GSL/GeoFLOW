//==================================================================================
// Module       : gstateinfo.hpp
// Date         : 6/3/20 (DLR)
// Description  : Header defining StateInfo structure
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GSTATEINFO_HPP)
#define _GSTATEINFO_HPP

#include "gtypes.h"
#include <stdio.h>
#include <cstdlib>
#include "gtmatrix.hpp"


struct GStateInfo {
  GINT        sttype  = 1;       // state type index (grid=0 or state=1)
  GINT        gtype   = 0;       // check src/cdg/include/gtypes.h
  GINT        nevolved= 0;       // number state comps evolved
  GINT        npresc  = 0;       // number state comps prescribed
  GSIZET      index   = 0;       // time index
  GSIZET      nelems  = 0;       // num elems
  GSIZET      cycle   = 0;       // continuous time cycle
  GFTYPE      time    = 0.0;     // state time
  std::vector<GString>
              svars;             // names of state members
  GTVector<GStateCompType>
              icomptype;         // encoding of state component types    
  GTMatrix<GINT>
              porder;            // if ivers=0, is 1 X GDIM; else nelems X GDIM;
  GString     idir;              // input directory
  GString     odir;              // output directory
};




#endif
