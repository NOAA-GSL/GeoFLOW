//==================================================================================
// Module       : ginit_factory
// Date         : 7/11/19 (DLR)
// Description  : GeoFLOW state initialization factory
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#include "ginitv_factory.hpp"
#include "ginitv_impl.hpp"

namespace geoflow {
namespace pdeint {


//**********************************************************************************
//**********************************************************************************
// METHOD : init
// DESC   : Do init of state components
// ARGS   : ptree  : main property tree
//          eqn_ptr: pointer to equation
//          time   : initialization time
//          ub     : boundary state (also initialized here)
//          u      : state to be initialized. Each component must be 
//                   labelled in EqnBasePtr::icomptype_.
// RETURNS: none.
//**********************************************************************************
void GInitFactory::static void init(const geoflow::tbox::PropertyTree& ptree, EqnBasePtr &eqn_ptr,  Time &time, State *ub, State &u)
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
} // end, init method



} // namespace pdeint
} // namespace geoflow

