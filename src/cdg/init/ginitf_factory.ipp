//==================================================================================
// Module       : ginitf_factory
// Date         : 7/11/19 (DLR)
// Description  : GeoFLOW forcing initialization factory
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#include "ginitf_factory.hpp"
#include "ginitf_impl.hpp"

namespace geoflow {
namespace pdeint {


//**********************************************************************************
//**********************************************************************************
// METHOD : init
// DESC   : Do init of forcing components
// ARGS   : ptree  : main property tree
//          eqn_ptr: pointer to equation
//          time   : initialization time
//          u      : current state vector
//          uf     : forcing components set from call
// RETURNS: none.
//**********************************************************************************
void GInitFFactory::static void init(const geoflow::tbox::PropertyTree& ptree, EqnBasePtr &eqn_ptr, GGrid &grid, Time &time, State &u, State &uf)
{
  GBOOL         bforced = ptree.getValue<GString>("use_forcing", FALSE);
  GString       sinit   = ptree.getValue<GString>("initf_block");
  PropertyTree  vtree   = ptree.getPropertyTree(sinit);

  if ( !bforced ) return;

  if      ( "initf_null"        == sinit ) {
    ginitf::initf_impl_null     (vtree, eqn_ptr, grid, time, ub, u);
  else if ( "initf_rand"        == sinit ) {
    ginitf::initf_impl_icosgauss(vtree, eqn_ptr, grid, time, ub, u);
  }
  else                                        {
    assert(FALSET & "Specified forcing initialization unknown");
  }

} // end, init method



} // namespace pdeint
} // namespace geoflow

