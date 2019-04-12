//==================================================================================
// Module       : gio.h
// Date         : 3/2/19 (DLR)
// Description  : Simple GeoFLOW IO
// Copyright    : Copyright 2019. Colorado State University. All rights reserved.
// Derived From : 
//==================================================================================
#if !defined(_GIO_H)
#define _GIO_H

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

extern GString  fname_;
extern char    *cfname_;
extern GINT     nfname_;

#define GIO_VERSION 0

struct GIOTraits {
  GBOOL        prgrid  = FALSE;
  GINT         wtime   = 6;
  GINT         wtask   = 5;
  GINT         wfile   = 2048;
  GINT         ivers   = 0;
  GINT         dim     = GDIM;
  GINT         gtype   = GE_REGULAR;
  GSIZET       index   = 0;
  GSIZET       nelems  = 0;
  GSIZET       cycle   = 0;
  GFTYPE       time    = 0.0;
  GString      dir     = ".";
  GTMatrix<GINT> porder;
};

void gio_write_state(GIOTraits &, GGrid &grid, 
                     const GTVector<GTVector<GFTYPE>*> &u, 
                     GTVector<GINT> &iu, 
                     const GTVector<GString> &svars, 
                     GC_COMM comm);

void gio_write_grid (GIOTraits &, GGrid &grid, 
                     const GTVector<GString> &svars, 
                     GC_COMM comm);

void  gio_restart(const geoflow::tbox::PropertyTree& ptree, GINT igrid, 
                  GTVector<GTVector<GFTYPE>*> &u, GTMatrix<GINT> &p,
                  GSIZET &icycle, GFTYPE &time, GC_COMM comm);

GSIZET gio_read_header(GIOTraits&, GString filename);

GSIZET gio_write(GIOTraits&, GString filename, const GTVector<GFTYPE> &u);

GSIZET gio_read (GIOTraits&, GString filename, GTVector<GFTYPE> &u);

void   gio_resize(GINT n);



#endif
