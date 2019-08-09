//==================================================================================
// Module       : ginitb_factory
// Date         : 7/11/19 (DLR)
// Description  : GeoFLOW state initialization factory
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================

namespace geoflow {
namespace pdeint {


//**********************************************************************************
//**********************************************************************************
// METHOD : init
// DESC   : Make call to to do bdy initialization
// ARGS   : ptree  : main property tree
//          time   : initialization time
//          utmp   : tmp arrays
//          u      : state to be initialized. 
//          ub     : boundary state 
// RETURNS: none.
//**********************************************************************************
template<typename EquationType>
GBOOL GInitBFactory<EquationType>::init(const geoflow::tbox::PropertyTree& ptree, GGrid &grid, Time &time, State &utmp, State &u, State &ub)
{
  GBOOL         bret=FALSE;
  GBOOL         use_inits; // use state init method to set bdy?
  GFTYPE        tt;
  GString       sinitb    = ptree.getValue<GString>("initb_block");
  GString       sinits    = ptree.getValue<GString>("inits_block");
  PropertyTree  vtreeb    = ptree.getPropertyTree(sinitb);
  PropertyTree  vtrees    = ptree.getPropertyTree(sinits);

  use_inits = vtree.getValue<GBOOL>("use_state_init",FALSE);

  tt = time;

  if ( "initb_none" == sinitb
    || "none"       == sinitb 
    || ""           == sinitb ) {
    bret = TRUE;
  }
  else if ( use_inits ) { // initialize from state initialization
    bret = ginitb::impl_bystateinit(ptree, grid, tt, utmp, u, ub);
  }
  else if ( "mybdyinit" == sinit ) { // initialize from user-specified fcn
    bret = ginitb::impl_mybdyinit  (ptree, grid, tt, utmp, u, ub);
  }
  else                                        {
    assert(FALSE && "Specified bdy initialization method unknown");
  }

  return bret;
} // end, init method init





} // namespace pdeint
} // namespace geoflow

