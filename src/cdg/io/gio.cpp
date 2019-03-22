//==================================================================================
// Module       : gio.cpp
// Date         : 3/2/19 (DLR)
// Description  : Simple GeoFLOW IO
// Copyright    : Copyright 2019. Colorado State University. All rights reserved.
// Derived From : 
//==================================================================================
#include "gio.h"


//**********************************************************************************
//**********************************************************************************
// METHOD : gout
// DESC   : Do simple GIO POSIX output
// ARGS   : grid  : GGrid object
//          u     : state
//          iu    : indirection indices, point to u members to print
//          tindex: time output index
//          time  : state evol. time
//          svars : variable names for output
//          sdir  : output directory
//          comm  : communicator
//          bprgrid: flag to print grid, which is not tagged by time index.
// RETURNS: none
//**********************************************************************************
void gio(GGrid &grid, const State &u, const GTVector<GINT> &iu, const GSIZET tindex, const GFTYPE time, const GTVector<GString> &svars, const GString &sdir, GC_COMM comm, const GBOOL bprgrid)
{

    GString serr ="gio: ";
    char    fname[2048];
    char    fn2;
    
    GINT myrank = GComm::WorldRank(comm);

    assert(svars.size() >=  iu.size()
        && "Insufficicent number of state variable names specified");

    GINT           dim    = GDIM;
    GSIZET         nelems = grid.nelems();
    GTVector<GINT> porder(GDIM);
    GTVector<GTVector<GFTYPE>>
                  *xnodes = &grid.xNodes();
    GElemList     *elems  = &grid.elems();
    for ( auto j=0; j<xnodes->size(); j++ ) porder[j] = (*elems)[0]->order(j);
   

    // Print convergence data to file:
    FILE     *fp;
    GString   sx[3] = {"xgrid", "ygrid", "zgrid"};
    GElemType gtype = grid.gtype();

    // Print grid, if not printed already:
    if ( bprgrid ) {
      for ( auto j=0; j<xnodes->size(); j++ ) {
        sprintf(fname, "%s/%s.%04d.out", sdir.c_str(), sx[j].c_str(),  myrank);
        fp = fopen(fname,"wb");

        assert(fp != NULLPTR && "gio.cpp: error opening file");
        // write header: dim, numelems, poly_order:
        fwrite(&dim         , sizeof(GINT)  ,    1, fp); // problem dimension
        fwrite(&nelems      , sizeof(GSIZET),    1, fp); // # elems
        fwrite(porder.data(), sizeof(GINT)  , GDIM, fp); // expansion order in each dir
        fwrite(&gtype       , sizeof(GINT)  ,    1, fp); // grid type
        fwrite(&time        , sizeof(GFTYPE),    1, fp); // time stamp
        // write this tasks field data:
        fwrite((*xnodes)[j].data(), sizeof(GFTYPE), (*xnodes)[j].size(), fp);
        fclose(fp);
      }
    }

    // Print field data:
    GINT j;
    GINT n = iu.size() == 0 ? u.size() : iu.size();
    for ( auto jj=0; jj<n; jj++ ) {
      j = iu.size() == 0 ? iu[jj] : jj;
      sprintf(fname, "%s/%s.%06d.%04d.out", sdir.c_str(), svars[j].c_str(), tindex, myrank);
      fp = fopen(fname,"wb");
      assert(fp != NULLPTR && "gio.cpp: error opening file");

      // write header: dim, numelems, poly_order:
      fwrite(&dim         , sizeof(GINT)  ,    1, fp);
      fwrite(&nelems      , sizeof(GSIZET),    1, fp);
      fwrite(porder.data(), sizeof(GINT)  , GDIM, fp);
      fwrite(&gtype       , sizeof(GINT)  ,    1, fp); // grid type
      fwrite(&time        , sizeof(GFTYPE),    1, fp); // time stamp
      // write this tasks field data:
      fwrite((*u[j]).data(), sizeof(GFTYPE), (*u[j]).size(), fp);
      fclose(fp);
    }

} // end, gio

