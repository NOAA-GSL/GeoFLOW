//==================================================================================
// Module       : gio.cpp
// Date         : 3/2/19 (DLR)
// Description  : Simple GeoFLOW IO
// Copyright    : Copyright 2019. Colorado State University. All rights reserved.
// Derived From : 
//==================================================================================
#include "gio.h"

using namespace std;

//**********************************************************************************
//**********************************************************************************
// METHOD : gout
// DESC   : Do simple GIO POSIX output
// ARGS   : grid  : GGrid object
//          u     : state
//          nu    : indirection indices
//          tindex: time output index
//          time  : state evol. time
//          svars : variable names for output
//          comm  : communicator
//          bprgrid: flag to print grid (not tagged by time index).
// RETURNS: none
//**********************************************************************************
GINT gio(GGrid &grid, State &u, GTVector<GINT> &nu, GSIZET tindex, GFTYPE time, GTVector<GString> &svars, GC_COMM comm, GBOOL &bprgrid)
{

    GString serr ="gio: ";
    char    fname[1024];
    
    GINT myrank = GComm::WorldRank(comm);

    assert(svars.size() >=  nu
        && "Insufficicent number of state variable names specified");

    GINT           dim = GDIM;
    GSIZET         nelems = grid.nelems();
    GTVector<GINT> porder(GDIM);
    GTVector<GTVector<GFTYPE>> *xnodes = &grid.xNodes();
    GElemList     *elems = &grid.elems();
    for ( auto j=0; j<xnodes->size(); j++ ) porder[j] = (*elems)[0]->order(j);
   

    // Print convergence data to file:
    FILE     *fp;
    GString   sx[3] = {"xgrid", "ygrid", "zgrid"};
    GElemType gtype = grid.gtype();

    // Print grid, if not printed already:
    if ( bprgrid ) {
      for ( auto j=0; j<xnodes->size(); j++ ) {
        sprintf(fname, "%s.%04d.out", sx[j].c_str(),  myrank);
        fp = fopen(fname,"wb");
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
      bprgrid = FALSE;
    }

    // Print field data:
    GINT j;
    GINT n = nu.size() == 0 ? nu.size() : u.size();
    for ( auto jj=0; jj<n; jj++ ) {
      j = nu.size() == 0 ? nu[jj] : jj;
      sprintf(fname, "%s.%06d.%04d.out", svars[j].c_str(), tindex, myrank);
      fp = fopen(fname,"wb");
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
    return(0);

} // end, gio

