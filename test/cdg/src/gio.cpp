//==================================================================================
// Module       : gio.cpp
// Date         : 3/2/19 (DLR)
// Description  : Simple GeoFLOW IO
// Copyright    : Copyright 2019. Colorado State University. All rights reserved.
// Derived From : 
//==================================================================================
#include "gio.h"

int gio(GGrid &grid, State &u, GTVector<GString> &svars, GComm comm, GBOOL &bgridprinted)
{

    GString serr ="gio: ";
    char    fname[1024];
    
    GINT myrank = GComm::WorldRank(comm);

    assert(svars.size() >= u.size() 
        && "Insufficicent number of state variable names specified");

    GINT           dim = GDIM;
    GSIZET         nelems = grid.nelems();
    GTVector<GINT> porder(GDIM);
    GElemList     *elems = &grid.elems();
    for ( auto j=0; j<xnodes->size(); j++ ) porder[j] = (*elems)[0].order(j);
   

    // Print convergence data to file:
    FILE *fp;
    GTVector<GTVector<GFTYPE>> *xnodes = &grid.xnodes();
    GString sx(3) = {"xgrid", "ygrid", "zgrid"};

    // Print grid, if not printed already:
    if ( !bgridprinted ) {
      for ( auto j=0; j<xnodes->size(); j++ ) {
        sprintf(fname, "%s.%04d.out", sx[j], myrank);
        fp = fopen(fname,"wb");
        // write header: dim, numelems, poly_order:
        fwrite(&dim         , sizeof(GINT)  ,    1, fp);
        fwrite(&nelems      , sizeof(GSIZET),    1, fp);
        fwrite(porder.data(), sizeof(GINT)  , GDIM, fp);
        // write this tasks field data:
        fwrite((*xnodes)[j].data(), sizeof(GFTYPE), (*xnodes)[j].size(), fp);
        fclose(fp);
      }
      bgridprinted = TRUE;
    }

    for ( auto j=0; j<xnodes->size(); j++ ) {
      sprintf(fname, "%s.%04d.out", svars[j], myrank);
      fp = fopen(fname,"wb");
      // write header: dim, numelems, poly_order:
      fwrite(&dim         , sizeof(GINT)  ,    1, fp);
      fwrite(&nelems      , sizeof(GSIZET),    1, fp);
      fwrite(porder.data(), sizeof(GINT)  , GDIM, fp);
      // write this tasks field data:
      fwrite((*u[j]).data(), sizeof(GFTYPE), (*u[j]).size(), fp);
      fclose(fp);
    }

    return(0);

} // end, gio

