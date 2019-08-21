//==================================================================================
// Module       : ginitf_factory
// Date         : 7/11/19 (DLR)
// Description  : GeoFLOW forcing initialization factory
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================

//**********************************************************************************
//**********************************************************************************
// METHOD : init
// DESC   : Do init of forcing components
// ARGS   : ptree  : main property tree
//          grid   : grid object
//          time   : initialization time
//          utmp   : tmp arrays
//          u      : current state vector
//          uf     : forcing components set from call
// RETURNS: none.
//**********************************************************************************
template<typename EquationType>
GBOOL GInitFFactory<EquationType>::init(const geoflow::tbox::PropertyTree& ptree, GGrid &grid, Time &time, State &utmp, State &u, State &uf)
{
  GBOOL         bret=FALSE;
  GBOOL         bforced = ptree.getValue<GBOOL>  ("use_forcing", FALSE);
  GString       sinit   = ptree.getValue<GString>("initf_block");
  PropertyTree  ftree   = ptree.getPropertyTree(sinit);

  if ( !bforced ) return TRUE;

  if ( "initf_none" == sinit
    || "none"       == sinit
    || ""           == sinit ) {
    bret = TRUE;
  }
  else if ( "initf_null"        == sinit ) {
    bret = ginitf::impl_null     (ftree, grid, time, utmp, u, uf);
  }
  else if ( "initf_rand"        == sinit ) {
    bret = ginitf::impl_rand     (ftree, grid, time, utmp, u, uf);
  }
  else {                                        {
    assert(FALSE && "Specified forcing initialization unknown");
  }

  return bret;

} // end, init method


