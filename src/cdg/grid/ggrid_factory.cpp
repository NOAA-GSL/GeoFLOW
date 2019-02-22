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
// ARGS   : ptree : property tree
//          gbasis: basis object
//          comm  : GC_Comm object
// RETURNS: GGrid object ptr
//**********************************************************************************
GGrid *GGridFactory::build(const geoflow::tbox::PropertyTree& ptree, GTVector<GNBasis<GCTYPE,GFTYPE>*> gbasis, GC_COMM &comm)
{

  GString gname = ptree.getValue<GString>("grid_name");

  GGrid *grid;
  if      ( gname.compare  ("grid_icos") == 0   // 2d or 3d Icos grid
      ||    gname.compare("grid_sphere") == 0 ) {
    grid = new GGridIcos(ptree, gbasis, comm);
    grid->grid_init();
  }
  else if ( gname.compare("grid_box") == 0 ) { // 2d or 3D Cart grid
    grid = new GGridBox(ptree, gbasis, comm);
    grid->grid_init();
  }
  else {
    assert(FALSE && "Invalid PropertyTree grid specification");
  }

  return grid;
}


