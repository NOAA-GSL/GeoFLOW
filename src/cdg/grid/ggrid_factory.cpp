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
GGrid *GGridFactory::build(const geoflow::tbox::PropertyTree& ptree, GTVector<GNBasis<GCTYPE,GFTYPE>*> gbasis, GC_COMM comm)
{

  GINT myrank  = GComm::WorldRank(comm);
  GINT nprocs  = GComm::WorldSize(comm);
        
  GString gname = ptree.getValue<GString>("grid_name");
  
  GGrid  *grid = new GGrid(comm);

  GTVector<GBdyType> bdytype;

  if      ( gname.compare("grid_icos") == 0 ) { // 2d ICOS grid
    assert(GDIM == 2 && "GDIM must be 2");
    GGridIcos gen_icos(ptree, gbasis, comm);
    gen_icos.do_grid(*grid, myrank);
  }
  else if ( gname.compare("grid_sphere") == 0 ) { // 3D grid build on ICOIS
    assert(GDIM == 3 && "GDIM must be 3");
    GGridIcos gen_icos(gbasis, comm);
    gen_icos.do_grid(*grid, myrank);
  }
  else if ( gname.compare("grid_box") == 0 ) { // 2- 3D Cart grid
    GGridBox gen_box(ptree, gbasis, comm);
    gen_box.do_grid(*grid, myrank);
  }
  else {
    assert(FALSE && "Invalid PropertyTree grid specification");
  }

  return grid;
}


