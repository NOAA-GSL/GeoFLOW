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
// METHOD : gio_write_state
// DESC   : Do simple GIO POSIX output of state
// ARGS   : traits: GIOTraits structure
//          grid  : GGrid object
//          u     : state
//          iu    : indirection indices to fields to write
//          svars : variable names for state members
//          comm  : communicator
// RETURNS: none
//**********************************************************************************
void gio_write_state(GIOTraits &traits, GGrid &grid, 
                     const GTVector<GTVector<GFTYPE>*> &u, 
                     GTVector<GINT>              &iu, 
                     const GTVector<GString> &svars, GC_COMM comm) 
{
    GString serr = "gio_write_state: ";
    FILE          *fp;
    GINT           myrank = GComm::WorldRank(comm);
    GString        fname;
    GSIZET         nb, nc, nd;
    GTVector<GTVector<GFTYPE>>
                  *xnodes = &grid.xNodes();
    GElemList     *elems  = &grid.elems();
    std::stringstream format;

    // Required number of coord vectors:
    nc = grid.gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;

    // Do some checks:
    assert(svars.size() >=  iu.size()
        && "Insufficient number of state variable names specified");
    assert(traits.gtype == grid.gtype()
        && "Incompatible grid type");

    traits.dim    = GDIM;
    traits.nelems = grid.nelems();
    traits.gtype  = grid.gtype();

    assert(nc ==  xnodes->size()
        && "Incompatible grid or coord  dimensions");

    // Build format string:
    fname.resize(traits.wfile);
    format    << "\%s/\%s.\%0" << traits.wtime << "\%0" << traits.wtask << "d.out";

    // Set porder vector depending on version:
    traits.porder.resize(traits.ivers == 0 ? 1: traits.nelems, GDIM);
    for ( auto i=0; i<traits.porder.size(1); i++ )  { // for each element
      for ( auto j=0; j<traits.porder.size(2); j++ ) traits.porder(i,j) = (*elems)[i]->order(j);
    }

    // Cycle over all fields, and write:
    for ( auto j=0; j<iu.size(); j++ ) {
      printf(fname.c_str(), format.str().c_str(), traits.dir.c_str(), 
             svars[j].c_str(), traits.index,  myrank);
      fp = fopen(fname.c_str(),"wb");
      if ( fp == NULL ) {
        cout << serr << "Error opening file: " << fname << endl;
        exit(1);
      }
      gio_write(traits, fname, *u[iu[j]]);
      fclose(fp);
    }

} // end, gio_write_state


//**********************************************************************************
//**********************************************************************************
// METHOD : gio_write_grid
// DESC   : Do simple GIO POSIX output of grid
// ARGS   : traits: GIOTraits structure
//          grid  : GGrid object
//          svars : variable names for grid components
//          comm  : communicator
// RETURNS: none
//**********************************************************************************
void gio_write_grid(GIOTraits &traits, GGrid &grid, const GTVector<GString> &svars, GC_COMM comm) 
{
    GString serr ="gio_write_grid: ";
    FILE          *fp;
    GINT           myrank = GComm::WorldRank(comm);
    GString        fname;
    GSIZET         nb, nc, nd;
    GTVector<GTVector<GFTYPE>>
                  *xnodes = &grid.xNodes();
    GElemList     *elems  = &grid.elems();
    std::stringstream format;

    // Required number of coord vectors:
    nc = grid.gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;

    // Do some checks:
    assert(svars.size() >=  xnodes->size()
        && "Insufficient number of grid variable names specified");
    assert(nc ==  xnodes->size()
        && "Incompatible grid or coord  dimensions");

    traits.dim    = GDIM;
    traits.nelems = grid.nelems();
    traits.gtype  = grid.gtype();


    // Build format string:
    fname.resize(traits.wfile);
    format    << "\%s/\%s.\%0" << traits.wtask << "d.out";

    // Set porder vector depending on version:
    traits.porder.resize(traits.ivers == 0 ? 1: traits.nelems, GDIM);
    for ( auto i=0; i<traits.porder.size(1); i++ )  { // for each element
      for ( auto j=0; j<traits.porder.size(2); j++ ) traits.porder(i,j) = (*elems)[i]->order(j);
    }

    // Cycle over all coords, and write:
    for ( auto j=0; j<xnodes->size(); j++ ) {
      printf(fname.c_str(), format.str().c_str(), traits.dir.c_str(), svars[j].c_str(),  myrank);
      fp = fopen(fname.c_str(),"wb");
      if ( fp == NULL ) {
        cout << serr << "Error opening file: " << fname << endl;
        exit(1);
      }
      gio_write(traits, fname, (*xnodes)[j]);
      fclose(fp);
    }

} // end, gio_write_grid


//**********************************************************************************
//**********************************************************************************
// METHOD : gio_restart
// DESC   : Read restart files for entire state, based on prop tree configuration.
//          
// ARGS   : ptree : main property tree
//          igrid : if > 0, uses grid format; else uses state format
//          tindex: time index to restart
//          u     : state object
//          p     : matrix of polynomial orders of size nelems X GDIM if 
//                  order isn't constant, and of size x X GDIM if it is
//          icycle: time cycle of input
//          t     : time stamp of input
//          comm  : GC_COMM object
// RETURNS: none
//**********************************************************************************
void gio_restart(const geoflow::tbox::PropertyTree& ptree, GINT igrid, 
                 GTVector<GTVector<GFTYPE>*> &u, GTMatrix<GINT> &p,
                 GSIZET &icycle, GFTYPE &time, GC_COMM comm)
{
  GINT                 myrank = GComm::WorldRank(comm);
  GINT                 ivers, nc, nr, nt;
  GElemType            igtype;
  GSIZET               itindex;
  GString              fname;
  GString              def;
  GString              stmp;
  GString              sgtype;
  GTVector<GString>    gobslist;
  GTVector<GString>    ggnames;
  std::vector<GString> defgnames = {"xgrid", "ygrid", "zgrid"};
  std::vector<GString> stdlist;
  std::stringstream    format;
  geoflow::tbox::PropertyTree inputptree;
  GIOTraits            traits;

  // Get restart index:
  itindex = ptree.getValue<GSIZET>   ("restart_index", 0);

  // Get observer list:
  stdlist = ptree.getArray<GString>("observer_list");
  gobslist = stdlist;

  // Get version number:
  def    = "constant";
  stmp   = ptree.getValue<GString>("exp_order_type", def);
  ivers  = stmp == "constant" ? 0 : 1;

  // Get configured grid type:
  def    = "grid_box";
  sgtype = ptree.getValue<GString>("grid_type", def);
  igtype = GE_REGULAR;
  if ( sgtype == "grid_icos" && GDIM == 2 ) igtype = GE_2DEMBEDDED;
  if ( sgtype == "grid_icos" && GDIM == 3 ) igtype = GE_DEFORMED;
  nc = GDIM;
  if ( igtype == GE_2DEMBEDDED ) nc = GDIM + 1;


  if ( gobslist.contains("posixio_observer") ) {
    inputptree    = ptree.getPropertyTree       ("posixio_observer");
    fname.resize(inputptree.getValue<GINT>("filename_size",2048));
    stdlist       = inputptree.getArray<GString>("grid_names", defgnames);
    ggnames       = stdlist;
    traits.wfile  = inputptree.getValue   <GINT>("filename_size",2048);
    traits.wtask  = inputptree.getValue   <GINT>("task_field_width", 5);
    traits.wtime  = inputptree.getValue   <GINT>("time_field_width", 6);
    def = ".";
    traits.dir    = inputptree.getValue<GString>("indirectory",def);

    // Get data:
    if ( igrid )  // use grid format
      format    << "\%s/\%s.\%0" << traits.wtask << "d.out";
    else          // use state format
      format    << "\%s/\%s.\%0" << traits.wtime << "d\%0" << traits.wtask << "d.out";
    
    for ( GSIZET j=0; j<nc; j++ ) { // Retrieve all grid coord vectors
      if ( igrid )  // use grid format
        printf(fname.c_str(), format.str().c_str(), traits.dir.c_str(), myrank);
      else          // use state format
        printf(fname.c_str(), format.str().c_str(), traits.dir.c_str(), itindex, myrank);
  
      nr = gio_read(traits, fname, *u[j]);
    }
    p.resize(traits.porder.size(1),traits.porder.size(2));
    p = traits.porder;

    if ( traits.ivers > 0 ) traits.ivers = 1;
    assert(traits.ivers == ivers  && "Incompatible version identifier");
    assert(traits.dim   == GDIM   && "Incompatible problem dimension");
    assert(traits.gtype == igtype && "Incompatible grid type");
  }
  else {
    assert( FALSE && "Configuration does not exist for grid files");
  }

} // end, gio_restart method


//**********************************************************************************
//**********************************************************************************
// METHOD : gio_write
// DESC   : Do simple (lowes-level) GIO POSIX output of specified field
// ARGS   : 
//          traits  : GIOTraits structure
//          fname   : filename, fully resolved
//          u       : field to output
//          icycle  : evol. time cycle
//          time    : evol. time
// RETURNS: none
//**********************************************************************************
GSIZET gio_write(GIOTraits &traits, GString fname, const GTVector<GFTYPE> &u) 
{

    GString  serr ="gio_write: ";
    FILE     *fp;
    GSIZET    i, j, nb, nd;

   
    fp = fopen(fname.c_str(),"wb");
    if ( fp == NULL ) {
      cout << serr << "Error opening file: " << fname << endl;
      exit(1);
    }

    if ( traits.ivers > 0 && traits.porder.size(1) <= 1 ) {
      cout << serr << " porder of insufficient size for version: " << traits.ivers
                   << ". Error writing file " << fname << endl;
      exit(1);
    }

    // Write header: dim, numelems, poly_order:
    fwrite(&traits.ivers       , sizeof(GINT)  ,    1, fp); // GIO version number
    fwrite(&traits.dim         , sizeof(GINT)  ,    1, fp); // dimension
    fwrite(&traits.nelems      , sizeof(GSIZET),    1, fp); // num elements
    nd = traits.porder.size(1) * traits.porder.size(2);
    fwrite(traits.porder.data().data()
                               , sizeof(GINT)  ,   nd, fp); // exp order
    fwrite(&traits.gtype       , sizeof(GINT)  ,    1, fp); // grid type
    fwrite(&traits.cycle       , sizeof(GSIZET),    1, fp); // time cycle stamp
    fwrite(&traits.time        , sizeof(GFTYPE),    1, fp); // time stamp

    // Write field data:
    nb = fwrite(u.data(), sizeof(GFTYPE), u.size(), fp);
    fclose(fp);

    return nb;

} // end, gio_write


//**********************************************************************************
//**********************************************************************************
// METHOD : gio_read
// DESC   : Do simple GIO POSIX read
// ARGS   : 
//          traits  : GIOTraits structure, filled here
//          fname   : filename, fully resolved
//          u       : field to output; resized if required
// RETURNS: none
//**********************************************************************************
GSIZET gio_read(GIOTraits &traits, GString fname, GTVector<GFTYPE> &u)
{

    GString serr = "gio_read: ";
    FILE     *fp;
    GSIZET    nb, nd, nh, nt;
    
    nh  = gio_read_header(traits, fname);

    assert(traits.dim == GDIM && "File dimension incompatible with GDIM");

    fp = fopen(fname.c_str(),"rb");
    if ( fp == NULL ) {
      cout << serr << "Error opening file: " << fname << endl;
      exit(1);
    }

    // Seek to first byte after header:
    fseek(fp, nh, SEEK_SET); 

    // Compute field data size from header data:
    nd = 0;
    if ( traits.ivers == 0 ) { // expansion order is constant
      nt = 1; 
      for ( GSIZET j=0; j<GDIM; j++ ) nt *= (traits.porder(0,j) + 1);
      nd += nt * traits.nelems;
    }
    else {                     // expansion order varies
      for ( GSIZET i=0; i<traits.porder.size(1); i++ ) {
        nt = 1; 
        for ( GSIZET j=0; j<GDIM; j++ ) nt *= (traits.porder(i,j) + 1);
        nd += nt;
      }
    }

    u.resize(nd);
    nb = fread(u.data(), sizeof(GFTYPE), nd, fp);

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
// ARGS   : traits: GIOTraits structure, filled with what header provides
//          fname : file name (fully resolved)
// RETURNS: no. header bytes read
//**********************************************************************************
GSIZET gio_read_header(GIOTraits &traits, GString fname)
{

    GString serr ="gio_read_header: ";
    GSIZET nb, nd, nh, numr;
  
    // Read field data:
    FILE *fp;
    nb = 0;
    fp = fopen(fname.c_str(),"rb");
    assert(fp != NULLPTR && "gio.cpp: error opening file");
  
    // Read header: 
    nh = fread(&traits.ivers       , sizeof(GINT)  ,    1, fp); nb += nh;
    nh = fread(&traits.dim         , sizeof(GINT)  ,    1, fp); nb += nh;
    nh = fread(&traits.nelems      , sizeof(GSIZET),    1, fp); nb += nh;
    numr = traits.ivers == 0 ? traits.dim : traits.nelems * traits.dim;
    nh = fread(traits.porder.data().data()
                                   , sizeof(GINT)  , numr, fp); nb += nh;
    nh = fread(&traits.gtype       , sizeof(GINT)  ,    1, fp); nb += nh;
    nh = fread(&traits.cycle       , sizeof(GSIZET),    1, fp); nb += nh;
    nh = fread(&traits.time        , sizeof(GFTYPE),    1, fp); nb += nh;
  
    fclose(fp);
  
    // Get no. bytest that should have been read:
    nd = (numr+3)*sizeof(GINT) + 2*sizeof(GSIZET) + sizeof(GFTYPE);
  
    if ( nb != nd ) {
      cout << serr << "Incorrect amount of data read from file: " << fname << endl;
      exit(1);
    }

    return nb;

} // end, gio_read_header

