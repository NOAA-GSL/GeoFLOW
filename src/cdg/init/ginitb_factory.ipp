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
// DESC   : Do bdy specification 
// ARGS   : ptree  : main property tree
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state 
//          u      : state to be initialized. 
// RETURNS: none.
//**********************************************************************************
template<typename EquationType>
GBOOL GInitBFactory<EquationType>::init(const geoflow::tbox::PropertyTree& ptree, GGrid &grid, Time &time, State &utmp, State &ub, State &u)
{
  GBOOL         bret=FALSE;
  GString       sinit   = ptree.getValue<GString>("initb_block");
  PropertyTree  vtree   = ptree.getPropertyTree(sinit);

  if ( "initb_none" == sinit
    || "none"       == sinit 
    || ""           == sinit ) {
    return TRUE;;
  }
  else                                        {
    assert(FALSE && "Specified bdy initialization method unknown");
  }

  return bret;
} // end, init method



} // namespace pdeint
} // namespace geoflow

