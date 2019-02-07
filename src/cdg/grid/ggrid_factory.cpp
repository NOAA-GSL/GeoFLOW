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
GGrid *GGridFactory::build(const tbox::PropertyTree& ptree, GTVector<GNBasis<GCTYPE,GFTYPE>*> gbasis, GC_COMM comm)
{

  GINT myrank  = GComm::WorldRank(comm);
  GINT nprocs  = GComm::WorldSize(comm);
        
  GString gname = ptree("grid_name");
  GGrid  *grid = new GGrid(comm);

  if      ( gname.compare("grid_icos") { // 2d ICOS grid
    assert(GDIM == 2 && "GDIM must be 3");
    GFTYPE radius = ptree.getValue<GFTYPE>("radius");
    GFTYPE ilevel = ptree.getValue<GINT>  ("ilevel");
    GGridIcos gen_icos(radiusi, ilevel, gbasis, nprocs);
    gen_icos.do_grid(*grid, myrank);
  }
  else if ( gname.compare("grid_sphere") { // 3D grid build on ICOIS
    assert(GDIM == 3 && "GDIM must be 3");
    GFTYPE radiusi = ptree.getValue<GFTYPE>("radiusi");
    GFTYPE radiuso = ptree.getValue<GFTYPE>("radiuso");
    std::vector<GINT> stdne = ptree.getArray("num_elems",0);
    GTVector<GINT> ne = stdne;
    GGridIcos gen_icos(radiusi, radiuso, ne, gbasis, nprocs);
    gen_icos.do_grid(*grid, myrank);
  }
  else if ( gname.compare("grid_box") { // 2- 3D Cart grid
    GFTYPE radiusi = ptree.getValue<GFTYPE>("radiusi");
    std::vector<GFTYPE> xyz0   = ptree.getArray("xyz0",0);
    std::vector<GFTYPE> delxyz = ptree.getArray("delxyz",0);
    GTPoint<GFTYPE> P0(3);  // start point, lower LHS diagonal
    GTPoint<GFTYPE> P1(3);  // end point, upper RHS diagonal 
    for ( GSIZET j=0; j<3; j++ ) {
      P0[j] = xyz0[j];
      P1[j] = p0[j] + delxyz[j];
    }
    std::vector<GINT> stdne = ptree.getArray("num_elems",0);
    GTVector<GINT> ne = stdne;
    GGridBox gen_box(P0, P1, ne, gbasis, nprocs);
    gen_box.do_grid(*grid, myrank);
  }
  else {
    assert(FALSE && "Invalid PropertyTree grid specification");
  }

  return grid;
}

