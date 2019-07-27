//==================================================================================
// Module       : gupdateb_factory
// Date         : 7/11/19 (DLR)
// Description  : GeoFLOW state initialization factory
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================

namespace geoflow {
namespace pdeint {


//**********************************************************************************
//**********************************************************************************
// METHOD : update
// DESC   : Do bdy update
// ARGS   : ptree  : main property tree
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state 
//          u      : state to be initialized. 
// RETURNS: none.
//**********************************************************************************
template<typename EquationType>
void GInitBFactory<EquationType>::update(const geoflow::tbox::PropertyTree& ptree, GGrid &grid, Time &time, State &utmp, State &ub, State &u)
{
  GString       sinit   = ptree.getValue<GString>("updateb_block");
  PropertyTree  vtree   = ptree.getPropertyTree(sinit);

  if ( "updateb_none" == sinit
    || "none"       == sinit 
    || ""           == sinit ) {
    return;
  }
  else                                        {
    assert(FALSE && "Specified bdy update method unknown");
  }

} // end, init method



} // namespace pdeint
} // namespace geoflow

