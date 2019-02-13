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
    assert(GDIM == 2 && "GDIM must be 3");
    bdytype.resize(2); bdytype = GBDY_NONE;
    GGridIcos::Traits icostraits;
    icostraits.ilevel  = ptree.getValue<GINT>("ilevel");
    icostraits.radiusi = ptree.getValue<GFTYPE>("radius");
    icostraits.radiuso = icostraits.radiusi;
    icostraits.bdyTypes= bdytype;
    GGridIcos gen_icos(icostraits, gbasis, nprocs);
    gen_icos.do_grid(*grid, myrank);
  }
  else if ( gname.compare("grid_sphere") == 0 ) { // 3D grid build on ICOIS
    assert(GDIM == 3 && "GDIM must be 3");
    bdytype.resize(2);
    bdytype[0] = str2bdytype(ptree.getValue<GString>("bdy_inner"));
    bdytype[1] = str2bdytype(ptree.getValue<GString>("bdy_outer"));
    GGridIcos::Traits sphtraits;
    sphtraits.radiusi = ptree.getValue<GFTYPE>("radiusi");
    sphtraits.radiuso = ptree.getValue<GFTYPE>("radiuso");
    sphtraits.bdyTypes= bdytype;
    std::vector<GINT> stdne = ptree.getArray<GINT>("num_elems");
    GTVector<GINT> ne(stdne.size());
    ne = stdne;
    GGridIcos gen_icos(sphtraits, ne, gbasis, nprocs);
    gen_icos.do_grid(*grid, myrank);
  }
  else if ( gname.compare("grid_box") == 0 ) { // 2- 3D Cart grid
    bdytype.resize(2*GDIM);
    bdytype[0] = str2bdytype(ptree.getValue<GString>("bdy_y_0"));
    bdytype[1] = str2bdytype(ptree.getValue<GString>("bdy_x_1"));
    bdytype[2] = str2bdytype(ptree.getValue<GString>("bdy_y_1"));
    bdytype[3] = str2bdytype(ptree.getValue<GString>("bdy_x_0"));
    if ( GDIM == 3 ) {
    bdytype[4] = str2bdytype(ptree.getValue<GString>("bdy_z_0"));
    bdytype[5] = str2bdytype(ptree.getValue<GString>("bdy_z_1"));
    }
    std::vector<GFTYPE> xyz0   = ptree.getArray<GFTYPE>("xyz0");
    std::vector<GFTYPE> delxyz = ptree.getArray<GFTYPE>("delxyz");
    std::vector  <GINT> stdne  = ptree.getArray<GINT>("num_elems");
    GTVector<GINT> ne(stdne.size());
    ne     = stdne;
    GGridBox::Traits boxtraits;
    GTPoint<GFTYPE> P0(3);  // start point, lower LHS diagonal
    GTPoint<GFTYPE> P1(3);  // end point, upper RHS diagonal 
    for ( GSIZET j=0; j<3; j++ ) {
      P0[j] = xyz0[j];
      P1[j] = P0[j] + delxyz[j];
    }
    boxtraits.P0 = P0;
    boxtraits.P1 = P1;
    boxtraits.bdyTypes= bdytype;
    GGridBox gen_box(boxtraits, ne, gbasis, nprocs);
    gen_box.do_grid(*grid, myrank);
  }
  else {
    assert(FALSE && "Invalid PropertyTree grid specification");
  }

  return grid;
}


//**********************************************************************************
//**********************************************************************************
// METHOD : str2bdytype
// DESC   : Convert string to GBdyType
// ARGS   : stype: string type
// RETURNS: GBdyType
//**********************************************************************************
GBdyType GGridFactory::str2bdytype(const GString &stype)
{
  GString s0;
  for ( auto j=0; j<GBDY_NONE; j++ ) {
    s0 = sGBdyType[j];
    if ( stype.compare(s0) == 0 ) return static_cast<GBdyType>(j);
  }
  assert(FALSE && "Invalid boundary type");
} // end, str2bdytype
