//==================================================================================
// Module       : ginitstate_factory
// Date         : 7/11/19 (DLR)
// Description  : GeoFLOW state initialization factory
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================

//**********************************************************************************
//**********************************************************************************
// METHOD : init
// DESC   : Do init of state components
// ARGS   : ptree  : main property tree
//          grid   : grid object
//          peqn   : ptr to EqnBase 
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          u      : state to be initialized. 
// RETURNS: none.
//**********************************************************************************
template<typename EquationType>
GBOOL GInitStateFactory<EquationType>::init(const geoflow::tbox::PropertyTree& ptree, GGrid &grid, EqnBasePtr &peqn, Time &time, State &utmp, State &ub, State &u)
{
  GBOOL         bret    = FALSE;
  GString       stype ;  
  GString       sinit   = ptree.getValue<GString>("inits_block");
  PropertyTree  vtree   = ptree.getPropertyTree(sinit);

  // Get type of initialization: direct or by-var:
  stype = vtree.getValue<GString>("init_type","");
  if ( "direct"   == stype 
    || ""         == stype ) {
    bret = set_direct(ptree, grid, peqn, time, utmp, ub, u);
  }
  else if ( "by-var" == stype ) {
    bret = set_by_grp(ptree, grid, peqn, time, utmp, ub, u);
  }
  else {
    assert(FALSE && "Invalid state initialization type");
  }

  return bret;

} // end, init method


//**********************************************************************************
//**********************************************************************************
// METHOD : set_direct
// DESC   : Do init of state components by calling init method to set
//          entire state at once. E.g., one might classify initialization
//          schemes by PDE-type, and user is responsible to ensuring 
//          all state members are initialized.
// ARGS   : ptree  : main property tree
//          grid   : grid object
//          peqn   : ptr to EqnBase 
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          u      : state to be initialized. 
// RETURNS: none.
//**********************************************************************************
template<typename EquationType>
GBOOL GInitSFactory<EquationType>::set_direct(const geoflow::tbox::PropertyTree& ptree, GGrid &grid, EqnBasePtr &peqn, Time &time, State &utmp, State &ub, State &u)
{
  GBOOL         bret    = FALSE;
  GString       sinit   = ptree.getValue<GString>("inits_block");
  PropertyTree  vtree   = ptree.getPropertyTree(sinit);
  if      ( "initstate_icosgaussburgers"  == sinit ) {
    bret = ginitstate::impl_icosgauss       (vtree, eqn_ptr, grid, time, utmp, ub, u);
  }
  else if ( "initstate_boxdirgauss"        == sinit ) {
    bret = ginitstate::impl_boxdirgauss     (vtree, eqn_ptr, grid, time, utmp, ub, u);
  }
  else if ( "initstate_boxpergauss"        == sinit ) {
    bret = ginitstate::impl_boxpergauss     (vtree, eqn_ptr, grid, time, utmp, ub, u);
  }
  else if ( "initstate_nwave"              == sinit ) {
    bret = ginitstate::impl_nwave           (vtree, eqn_ptr, grid, time, utmp, ub, u);
  }
  else                                        {
    assert(FALSE & "Specified state initialization method unknown");
  }

  return bret;

} // end, set_direct method


//**********************************************************************************
//**********************************************************************************
// METHOD : set_by_grp
// DESC   : Do init of state components by specifying individual
//          variable group types, and initializing each group
//          in the state separately. This method uses the CompDesc data
//          in the EqnBase pointer to locate variable groups.
// ARGS   : ptree  : main property tree
//          grid   : grid object
//          peqn   : ptr to EqnBase 
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          u      : state to be initialized. 
// RETURNS: none.
//**********************************************************************************
template<typename EquationType>
GBOOL GInitSFactory<EquationType>::set_by_grp(const geoflow::tbox::PropertyTree& ptree, GGrid &grid, EqnBasePtr &peqn, Time &time, State &utmp, State &ub, State &u)
{
  GBOOL           bret    = FALSE;
  GSISET          ncomp;
  GString         sinit   = ptree.getValue<GString>("inits_block");
  GStateCompType *tcomp;
  PropertyTree    vtree   = ptree.getPropertyTree(sinit);
  CompDesc       *icomptype = &peqn.compdesc();

  ncomp = 0;
  tcomp = NULLPTR;

  icomptype->distinct(tcomp, ncomp);
  

  if      ( "initv"  == sinit ) {
    bret = ginitstate::impl_init_null  (vtree, eqn_ptr, grid, time, utmp, ub, u);
  }

  if ( tcomp  != NULLPTR ) delete [] tcomp;

  return bret;

} // end, set_by_grp method

