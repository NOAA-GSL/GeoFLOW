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
// METHOD : gio_write
// DESC   : Do simple GIO POSIX output
// ARGS   : grid  : GGrid object
//          u     : state
//          iu    : indirection indices, point to u members to print
//          tindex: time output index
//          time  : state evol. time
//          svars : variable names for output
//          sdir  : output directory
//          ivers : version number
//          comm  : communicator
//          bprgrid: flag to print grid, which is not tagged by time index.
// RETURNS: none
//**********************************************************************************
template <typename T>
void gio_write(GGrid &grid, const GTVector<GTVector<T>*> &u, const GTVector<GINT> &iu, const GSIZET tindex, const GFTYPE time, const GTVector<GString> &svars, const GString &sdir, GINT ivers, GC_COMM comm, const GBOOL bprgrid)
{

    GString serr ="gio_write: ";
    char    fname[2048];
    char    fn2;
    
    GINT myrank = GComm::WorldRank(comm);

    assert(svars.size() >=  iu.size()
        && "Insufficient number of state variable names specified");

    GINT           dim    = GDIM;
    GSIZET         nelems = grid.nelems();
    GTMatrix<GINT> porder;
    GTVector<GTVector<GFTYPE>>
                  *xnodes = &grid.xNodes();
    GElemList     *elems  = &grid.elems();

    porder.resize(1,ivers == 0 ? GDIM : nelems);
    for ( auto i=0; i<porder.size(1); i++ )  { // for each element
      for ( auto j=0; j<porder->size(2); j++ ) porder(i,j) = (*elems)[i]->order(j);
    }
   
    // Print convergence data to file:
    FILE     *fp;
    GString   sx[3] = {"xgrid", "ygrid", "zgrid"};
    GElemType gtype = grid.gtype();

    // Print grid, if not printed already:
    if ( bprgrid ) {
      for ( auto j=0; j<xnodes->size(); j++ ) {
        sprintf(fname, "%s/%s.%05d.out", sdir.c_str(), sx[j].c_str(),  myrank);
        fp = fopen(fname,"wb");
        if ( nb != nd ) {
          cout << serr << "Error opening file: " << fname << endl;
          exit(1);
        }
        
        // write header: dim, numelems, poly_order:
        fwrite(&ivers       , sizeof(GINT)  ,    1, fp); // GIO version number
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
      sprintf(fname, "%s/%s.%06d.%05d.out", sdir.c_str(), svars[j].c_str(), tindex, myrank);
      fp = fopen(fname,"wb");
       if ( nb != nd ) {
         cout << serr << "Error opening file: " << fname << endl;
         exit(1);
       }

      // write header: dim, numelems, poly_order:
      fwrite(&ivers       , sizeof(GINT)  ,    1, fp); // GIO version number
      fwrite(&dim         , sizeof(GINT)  ,    1, fp);
      fwrite(&nelems      , sizeof(GSIZET),    1, fp);
      fwrite(porder.data(), sizeof(GINT)  , GDIM, fp);
      fwrite(&gtype       , sizeof(GINT)  ,    1, fp); // grid type
      fwrite(&time        , sizeof(GFTYPE),    1, fp); // time stamp
      // write this tasks field data:
      fwrite((*u[j]).data(), sizeof(GFTYPE), (*u[j]).size(), fp);
      fclose(fp);
    }

} // end, gio_write


//**********************************************************************************
//**********************************************************************************
// METHOD : gio_read
// DESC   : Do simple GIO POSIX input
// ARGS   : fname : file name (fully resolved)
//          comm  : communicator
//          u     : state member vector, returned, resized to exact required size,
//                  if necessary
// RETURNS: none
//**********************************************************************************
template <typename T>
GSIZET gio_read(const GString fname, GC_COMM comm, GTVector<T> &u)
{

    GString serr ="gio_read: ";
    GINT           dim, ivers;
    GElemType      gtype;
    GSIZET         nb, nd, nh, nelems, nt;
    GFTYPE         time;
    GTMatrix<GINT> porder(GDIM);
    
    nh  = gio_read_header(fname, comm, ivers, GINT dim, nelems, porder, gtype, time); 

    assert(dim == GDIM && "File dimension incompatible with GDIM");

    // Print convergence data to file:
    FILE *fp;
    fp = fopen(fname.c_str(),"wb");
    if ( nb != nd ) {
      cout << serr << "Error opening file: " << fname << endl;
      exit(1);
    }

    // Seek to first byte after header:
    fseek(fp, nh, SEEK_SET); 

    // Compute data size from header data:
    nd = 0;
    if ( ivers == 0 ) {
      nt = 1; 
      for ( GSIZET j=0; j<GDIM; j++ ) nt *= (porder(0,j) + 1);
      nd += nt * nelems;
    }
    else {
      for ( GSIZET i=0; i<porder.size(1); i++ ) {
        nt = 1; 
        for ( GSIZET j=0; j<GDIM; j++ ) nt *= (porder(i,j) + 1);
        nd += nt;
      }
    }

    u.resize(nd);
    nb = fread((*u[j]).data(), sizeof(T), ns, fp);

    fclose(fp);

    if ( nb != nd ) {
      cout << serr << "Incorrect amount of data read from file: " << fname << endl;
      exit(1);
    }

   return nb;

} // end, gio_read


//**********************************************************************************
//**********************************************************************************
// METHOD : gio_read_header
// DESC   : Read GIO POSIX file header
// ARGS   : fname : file name (fully resolved)
//          comm  : communicator
//          ivers : version number
//          dim   : problem dimensionality
//          nelems: number elements represented
//          porder: matrix of dimensions nelem X dim containing expansion order
//                  for each element, if ivers > 0. If ivers == 0, then exp.
//                  order is assumed to be constant, and matrix size is just
//                  1 X GDIM.
//          gtype : GElemType grid type
//          time  : GFTYPE time
// RETURNS: no header bytes read
//**********************************************************************************
template <typename T>
GSIZET gio_read_header(const GString fname, GC_COMM comm, GINT &ivers, GINT &dim, GSIZET &nelems, 
                GTMatrix<GINT> porder, GElemType &gtype, GFTYPE &time)
{

    GString serr ="gio_read_header: ";
    GSIZET nb, nd, nh, numr;
  
    // Read field data:
    nb = 0;
    fp = fopen(fname.c_str(),"wb");
    assert(fp != NULLPTR && "gio.cpp: error opening file");
  
    // Read header: dim, numelems, poly_order:
    nh = fread(&ivers       , sizeof(GINT)  ,    1, fp); nb += nd;
    nh = fread(&dim         , sizeof(GINT)  ,    1, fp); nb += nd;
    nh = fread(&nelems      , sizeof(GSIZET),    1, fp); nb += nd;
    numr = ivers == 0 ? dim : nelems*dim;;
    nh = fread(porder.data(), sizeof(GINT)  , numr, fp); nb += nd;
    nh = fread(&gtype       , sizeof(GINT)  ,    1, fp); nb += nd;
    nh = fread(&time        , sizeof(GFTYPE),    1, fp); nb += nd;
  
    fclose(fp);
  
    // Get no. bytest that should have been read:
    nd = (numr+3)*sizeof(GINT) + sizeof(GSIZET) + sizeof(GFTYPE);
  
    if ( nb != nd ) {
      cout << serr << "Incorrect amount of data read from file: " << fname << endl;
      exit(1);
    }

    return nb;

} // end, gio_read_header

