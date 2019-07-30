//==================================================================================
// Module       : gbdy_factory
// Date         : 7/11/19 (DLR)
// Description  : GeoFLOW bdy config & initialization factory
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================

namespace geoflow {
namespace pdeint {


//**********************************************************************************
//**********************************************************************************
// METHOD : init
// DESC   : Do init of forcing components
// ARGS   : ptree  : main property tree
//          eqn_ptr: pointer to equation
//          time   : initialization time
//          utmp   : tmp arrays
//          u      : current state vector
//          ub     : bdy components set from call
// RETURNS: none.
//**********************************************************************************
void GBdyFactory::static void init(const geoflow::tbox::PropertyTree& ptree, GGrid &grid, Time &time, State &utmp, State &u, State &ub)
btime_dep_           (FALSE),
update_bdy_callback_ (NULLPTR)
{
  GBOOL         bforced = ptree.getValue<GString>("use_forcing", FALSE);
  GString       sinit   = ptree.getValue<GString>("initf_block");
  PropertyTree  ftree   = ptree.getPropertyTree(sinit);


  if      ( "initf_null"        == sinit ) {
    ginitf::impl_null     (ftree, eqn_ptr, grid, time, utmp, u, ub);
  else if ( "initf_rand"        == sinit ) {
    ginitf::impl_rand     (ftree, eqn_ptr, grid, time, utmp, u, ub);
  }
  else                                        {
    assert(FALSE & "Specified forcing initialization unknown");
  }

} // end, init method



} // namespace pdeint
} // namespace geoflow
