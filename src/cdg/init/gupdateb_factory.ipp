//==================================================================================
// Module       : gupdateb_factory
// Date         : 7/11/19 (DLR)
// Description  : GeoFLOW state initialization factory
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================

//**********************************************************************************
//**********************************************************************************
// METHOD : update
// DESC   : Do bdy update
// ARGS   : ptree  : main property tree
//          time   : initialization time
//          utmp   : tmp arrays
//          u      : state to be initialized. 
//          ub     : boundary state 
// RETURNS: none.
//**********************************************************************************
template<typename EquationType>
void GUpdateBFactory<EquationType>::update(const geoflow::tbox::PropertyTree& ptree, GGrid &grid, Time &time, State &utmp, State &u, State &ub)
{
  GBOOL         bret = FALSE, use_inits;
  Time          tt;
  GString       sgrid, supdate;
  PropertyTree  gtree;

  sgrid    = ptree.getValue<GString>("grid_type");
  gtree    = ptree.getPropertyTree(sgrid);
  supdate  = gtree.getValue<GStrig>("bdy_update_method","none");
  use_inits= gtree.getValue<GStrig>("use_state_init",FALSE);
  tt       = time;

  if ( "updateb_none" == supdate
    || "none"         == supdate
    || ""             == supdate ) {
    bret = TRUE;
  }
  else if ( use_inits ) {
    bret = gupdateb::impl_bystateinit(ptree, grid, tt, utmp, u, ub);
  }
  else if ( "simple_outflow" == supdate ) {
    bret = gupdateb::impl_simple_output  (ptree, grid, tt, utmp, u, ub);
  }

  else if ( "mybdyupdate"    == supdate ) {
    bret = gupdateb::impl_mybdyupdate   (ptree, grid, tt, utmp, u, ub);
  }
  else                                        {
    assert(FALSE && "Specified bdy update method unknown");
  }

  return bret;
} // end, init method update

