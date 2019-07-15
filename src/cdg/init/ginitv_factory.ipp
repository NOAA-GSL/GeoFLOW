//==================================================================================
// Module       : ginitv_factory
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
void GInitVFactory::static void init(const geoflow::tbox::PropertyTree& ptree, EqnBasePtr &eqn_ptr, GGrid &grid, Time &time, State &ub, State &u)
{
  GString       sinit   = ptree.getValue<GString>("initv_block");
  PropertyTree  vtree   = ptree.getPropertyTree(sinit);

  if      ( "initv_icosgauss"      == sinit ) {
    ginitv::initv_impl_icosgauss(vtree, eqn_ptr, grid, time, ub, u);
  }
  else if ( "initv_lump"           == sinit ) {
    ginitv::initv_impl_lump     (vtree, eqn_ptr, grid, time, ub, u);
  }
  else if ( "initv_nwave"          == sinit ) {
    ginitv::initv_impl_nwave    (vtree, eqn_ptr, grid, time, ub, u);
  }
  else                                        {
    assert(FALSET & "Specified state initialization unknown");
  }

} // end, init method



} // namespace pdeint
} // namespace geoflow

