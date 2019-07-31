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
//          u      : state to be initialized. 
//          ub     : boundary state 
// RETURNS: none.
//**********************************************************************************
template<typename EquationType>
void GUpdateBFactory<EquationType>::update(const geoflow::tbox::PropertyTree& ptree, GGrid &grid, Time &time, State &utmp, State &u, State &ub)
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

} // end, init method update



} // namespace pdeint
} // namespace geoflow

