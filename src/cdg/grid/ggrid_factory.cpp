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
#include "gio.h"


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
  GSIZET  itindex = ptree.getValue<GSIZET>   ("restart_index", 0);
  GString sdef    = "grid_box";
  GString gname   = ptree.getValue<GString>("grid_type", sdef);
  sdef            = "constant";
  GString ptype   = ptree.getValue<GString>("exp_order_type", sdef);
  geoflow::tbox::PropertyTree gridptree = ptree.getPropertyTree(gname);
  GGrid *grid;
  GTMatrix<GINT> p;
  GTVector<GTVector<GFTYPE>> xnodes;



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
      grid = new GGridIcos(gridptree, gbasis, comm);
      grid->grid_init();
    }
    else if ( "grid_box"    ==  gname) { // 2d or 3D Cart grid
      grid = new GGridBox(gridptree, gbasis, comm);
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
    read_grid(ptree, comm, p, xnodes);
    if      ( "grid_icos"   == gname   // 2d or 3d Icos grid
        ||    "grid_sphere" == gname ) {
      grid = new GGridIcos(gridptree, gbasis, comm);
      grid->grid_init(p, xnodes);
    }
    else if ( "grid_box"    ==  gname) { // 2d or 3D Cart grid
      grid = new GGridBox(gridptree, gbasis, comm);
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
void GGridFactory::read_grid(const geoflow::tbox::PropertyTree& ptree, GC_COMM comm,  GTMatrix<GINT> &p, 
                             GTVector<GTVector<GFTYPE>> &xnodes)
{
  GINT                        nc=GDIM;
  GElemType                   igtype;
  GSIZET                      cycle;
  GFTYPE                      time;
  GString                     sgtype;
  GTVector<GTVector<GFTYPE>*> u;


  // Resize xnodes based on configuration file:
  sgtype = ptree.getValue<GString>("grid_type");
  igtype = GE_REGULAR;
  if ( sgtype == "grid_icos" && GDIM == 2 ) igtype = GE_2DEMBEDDED;
  if ( sgtype == "grid_icos" && GDIM == 3 ) igtype = GE_DEFORMED;
  if ( igtype == GE_2DEMBEDDED ) nc = GDIM + 1;
  xnodes.resize(nc);

  u.resize(nc);
  for ( GSIZET j=0; j<nc; j++ ) u[j] = &xnodes[j];

  // Read restart data; time and  cycle are dummy here:
  gio_restart(ptree, 1, u, p, cycle, time, comm);

} // end, read_grid method
