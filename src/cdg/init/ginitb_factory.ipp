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
// DESC   : Get method to do bdy initialization. This function should be called
//          prior to calling the state initialization.
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
  std::function<void(PropertyTree &ptree, GGrid &grid, const Time &t,
                           State &utmp, State &u   , State &ub)>
          initb;

  initb = initfcn(ptree, grid, time, utmp, u, ub);

  return bret;
} // end, init method init


//**********************************************************************************
//**********************************************************************************
// METHOD : initfcn
// DESC   : Get function pointer to do bdy initialization
// ARGS   : ptree  : main property tree
//          time   : initialization time
//          utmp   : tmp arrays
//          u      : state to be initialized. 
//          ub     : boundary state 
// RETURNS: none.
//**********************************************************************************
template<typename EquationType>
std::function<void(PropertyTree &ptree, GGrid &grid, const Time &t,
                           State &utmp, State &u   , State &ub)>
GInitBFactory<EquationType>::initfcn(const geoflow::tbox::PropertyTree& ptree, GGrid &grid, Time &time, State &utmp, State &u, State &ub)
{
  GBOOL         bret=FALSE;
  GBOOL         use_inits; // use state init method to set bdy?
  GString       sinitb    = ptree.getValue<GString>("initb_block");
  GString       sinits    = ptree.getValue<GString>("inits_block");
  PropertyTree  vtreeb    = ptree.getPropertyTree(sinitb);
  PropertyTree  vtrees    = ptree.getPropertyTree(sinits);

  use_inits = vtree.getValue<GBOOL>("use_state_init",FALSE);

  static std::function<void(PropertyTree &ptree, GGrid &grid, const Time &t,
                           State &utmp, State &u   , State &ub)> fret =
                           this->doinits;

  if ( "initb_none" == sinit
    || "none"       == sinit 
    || ""           == sinit ) {
    return NULLPTR;
  }
  else if ( use_inits ) { // initialize from state initialization
    return fret;
  }
  else                                        {
    assert(FALSE && "Specified bdy initialization method unknown");
  }

  return NULLPTR;
} // end, init method initfcn


//**********************************************************************************
//**********************************************************************************
// METHOD : doinits
// DESC   : method to do bdy initialization from state initialization
// ARGS   : ptree  : main property tree
//          time   : initialization time
//          utmp   : tmp arrays
//          u      : state to be initialized. 
//          ub     : boundary state 
// RETURNS: none.
//**********************************************************************************
template<typename EquationType>
std::function<void(PropertyTree &ptree, GGrid &grid, const Time &t,
                           State &utmp, State &u   , State &ub)>
GInitBFactory<EquationType>::doinits(const geoflow::tbox::PropertyTree& ptree, GGrid &grid, Time &time, State &utmp, State &u, State &ub)
{
  GBOOL         bret=FALSE;

  return bret;

} // end, method doinits;




} // namespace pdeint
} // namespace geoflow

