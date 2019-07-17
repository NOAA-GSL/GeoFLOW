//==================================================================================
// Module       : ginits_factory
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
// DESC   : Do init of state components
// ARGS   : ptree  : main property tree
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          u      : state to be initialized. 
// RETURNS: none.
//**********************************************************************************
void GInitSFactory::static void init(const geoflow::tbox::PropertyTree& ptree, GGrid &grid, Time &time, State &utmp, State &ub, State &u)
{
  GString       sinit   = ptree.getValue<GString>("inits_block");
  PropertyTree  vtree   = ptree.getPropertyTree(sinit);

  if      ( "inits_icosgaussburgers"  == sinit ) {
    ginits::impl_icosgauss       (vtree, eqn_ptr, grid, time, utmp, ub, u);
  }
  else if ( "inits_boxdirgauss"        == sinit ) {
    ginits::impl_boxdirgauss     (vtree, eqn_ptr, grid, time, utmp, ub, u);
  }
  else if ( "inits_boxpergauss"        == sinit ) {
    ginits::impl_boxpergauss     (vtree, eqn_ptr, grid, time, utmp, ub, u);
  }
  else if ( "inits_nwave"              == sinit ) {
    ginits::impl_nwave           (vtree, eqn_ptr, grid, time, utmp, ub, u);
  }
  else                                        {
    assert(FALSE & "Specified state initialization unknown");
  }

} // end, init method



} // namespace pdeint
} // namespace geoflow

