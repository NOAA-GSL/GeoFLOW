//==================================================================================
// Module       : ggrid_factory
// Date         : 2/1/19 (DLR)
// Description  : GeoFLOW grid factory object. 
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#include "ggrid_factory.hpp"
#include "ggrid_icos.hpp"
#include "ggrid_box.hpp"


//**********************************************************************************
//**********************************************************************************
// METHOD : build
// DESC   : Do build and return of GGrid object
// ARGS   : ptree : main property tree
//          gbasis: basis object
//          comm  : GC_Comm object
// RETURNS: GGrid object ptr
//**********************************************************************************
GGrid *GGridFactory::build(const geoflow::tbox::PropertyTree& ptree, GTVector<GNBasis<GCTYPE,GFTYPE>*> gbasis, GC_COMM &comm)
{
  GINT    itindex = ptree.getValue<GINT>   ("restart_index", 0);
  GString sdef    = "grid_box",
  GString gname   = ptree.getValue<GString>("grid_type", sdef);
  sdef            = "constant",
  GString ptype   = ptree.getValue<GString>("exp_order_type", sdef);
  geoflow::tbox::PropertyTree& gridptree = ptree.getPropertyTree(gname);
  GGrid *grid;
  GMatrix<GINT> p;
  GTVector<GTVector<GFTYPE>> xnodes;

  comm_ = comm;


  // NOTE: If doing a restart, and if exp_order_type = constant, then we build 
  //       grid as usual, based on prop tree. If doing a restart and
  //       exp_order_type = variable, then we read old grid from data file
  //       and build grid from that:

  if ( itindex == 0 
    || "constant" == ptype ) { // not doing a restart, or p doesn't change

    // In this case, gbasis is assumed to contain the basis
    // functions for all elements; these are assumed to be 
    // constant:
    if      ( "grid_icos"   == gname   // 2d or 3d Icos grid
        ||    "grid_sphere" == gname ) {
      grid = new GGridIcos(gridptree, gbasis, comm_);
      grid->grid_init();
    }
    else if ( "grid_box"    ==  gname) { // 2d or 3D Cart grid
      grid = new GGridBox(gridptree, gbasis, comm_);
      grid->grid_init();
    }
    else {
      assert(FALSE && "Invalid PropertyTree grid specification");
    }

  }
  else {                       // doing restart w/ variable p

    // In this case, gbasis is interpreted as a 'pool' of 
    // basis functions with various orders. It is an error
    // if correct order is not found on restart:
    read_grid(ptree, p, xnodes);
    if      ( "grid_icos"   == gname   // 2d or 3d Icos grid
        ||    "grid_sphere" == gname ) {
      grid = new GGridIcos(gridptree, gbasis, comm_);
      grid->grid_init(p, xnodes);
    }
    else if ( "grid_box"    ==  gname) { // 2d or 3D Cart grid
      grid = new GGridBox(gridptree, gbasis, comm_);
      grid->grid_init(p, xnodes);
    }
    else {
      assert(FALSE && "Invalid PropertyTree grid specification");
    }

  }

  return grid;
} // end, build method


//**********************************************************************************
//**********************************************************************************
// METHOD : read_grid
// DESC   : Read node positions from info provided in ptree
// ARGS   : ptree : main property tree
//          p     : matrix of poly exp. order in each direction for each element.
//                  Returned quantity has dimensions NElems X GDIM
//          xnodes: node positions from file. Returned quantity is a vector
//                  of 2 or 3 vectors representing x, y, ... coordinates. 
//                  Is resized according to the input data.
// RETURNS: GGrid object ptr
//**********************************************************************************
void GGridFactory::read_grid(const geoflow::tbox::PropertyTree& ptree, GTMatrix<GFTYPE> &p, 
                             GTVector<GTVector<GFTYPE>> &xnodes)
{
  GINT                 myrank = GCOMM::WorldRank(comm_);
  GINT                 dim, ivers, nd;
  GElemType            gtype;
  GTMatrix<GINT>       porder;
  char                 sfname[2048];
  GFTYPE               time;
  GString              def = ".";
  GString              sfname;
  GString              sin;
  GTVector<GString>    gobslist;
  std::vector<GString> stdobslist;
  geoflow::tbox::PropertyTree& inputptree;
  stdobslist = ptree.getArray<GString>("observer_list");
  gobslist = stdobslist;

  // Verify size of xnodes:
  printf(fname, "%s/%s.%05d.out", sin.c_str(), "xgrid",  myrank);
  gio_read_header(fname, comm_, ivers, dim, nelems, porder, gtype, time);
  nd = gtype == GE_2DEMBEDDED ? dim + 1 : dim;
  xnodes.resize(nd);

  GString   sx[3] = {"xgrid", "ygrid", "zgrid"};
  if ( gobslist.contains("posixio_observer") ) { 
    inputptree = ptree.getPropertyTree("posixio_observer");
    sin        = inputptree.getValue<GString>("indirectory",def);
    for ( j=0; j<nd; j++ ) { // Retrieve all grid coord vectors
      sfname      = sin + "/" + sx[j];
      printf(fname, "%s/%s.%05d.out", sin.c_str(), sfname.c_str(),  myrank);
      sfname.assign(fname);
      gio_read(sfname, comm_, xnodes[j]);
    }
  }
  else {
    assert( FALSE && "Configuration does not exist for grid files");
  }

} // end, read_grid method
