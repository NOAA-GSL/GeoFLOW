//==================================================================================
// Module       : gspecb_factory
// Date         : 7/11/19 (DLR)
// Description  : GeoFLOW boundary specification factory
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================

namespace geoflow {
namespace pdeint {


//**********************************************************************************
//**********************************************************************************
// METHOD : init
// DESC   : Do bdy initialization
// ARGS   : ptree  : main property tree
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state 
//          u      : state to be initialized. 
// RETURNS: none.
//**********************************************************************************
void GSpecBFactory::static void init(const geoflow::tbox::PropertyTree& ptree, GGrid &grid, Time &time, State &utmp, State &ub, State &u)
{
  GString       sinit   = ptree.getValue<GString>("specb_block");
  PropertyTree  vtree   = ptree.getPropertyTree(sinit);

  if ( "specb_none" == sinit
    || "none"       == sinit
    || ""           == init ) {
    return;
  }
  else if ( "specb_icosgaussburgers"  == sinit ) {
    gspecb::impl_icosgauss       (vtree, eqn_ptr, grid, time, utmp, ub, u);
  }
  else if ( "specb_nwave"              == sinit ) {
    gspecb::impl_nwave           (vtree, eqn_ptr, grid, time, utmp, ub, u);
  }
  else                                        {
    assert(FALSE && "Specified state initialization method unknown");
  }

} // end, init method



} // namespace pdeint
} // namespace geoflow

