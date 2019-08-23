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
// RETURNS: TRUE on success; else FALSE
//**********************************************************************************
template<typename EquationType>
GBOOL GInitForceFactory<EquationType>::init(const PropertyTree& ptree, GGrid &grid, EqnBasePtr &peqn, Time &time, State &utmp, State &ub, State &u)
{
  GBOOL         bret    = FALSE;
  GString       stype ;  
  GString       sinit   = ptree.getValue<GString>("initforce_block");
  PropertyTree  vtree   = ptree.getPropertyTree(sinit);

  // Get type of initialization: direct or by-var:
  stype = vtree.getValue<GString>("init_type","");
  if ( "name"   == stype 
    || ""       == stype ) {
    bret = set_by_direct(ptree, grid, peqn, time, utmp, ub, u);
  }
  else if ( "block" == stype ) {
    bret = set_by_comp  (ptree, grid, peqn, time, utmp, ub, u);
  }
  else {
    assert(FALSE && "Invalid state initialization type");
  }

  return bret;

} // end, init method


//**********************************************************************************
//**********************************************************************************
// METHOD : set_by_direct
// DESC   : Do init of state components by calling initstate_block by name,
//          E.g., one might classify initialization
//          schemes by PDE-type, and user is responsible to ensuring 
//          all state members are initialized.
// ARGS   : ptree  : main property tree
//          grid   : grid object
//          peqn   : ptr to EqnBase 
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          u      : state to be initialized. 
// RETURNS: TRUE on success; else FALSE
//**********************************************************************************
template<typename EquationType>
GBOOL GInitForceFactory<EquationType>::set_by_direct(const PropertyTree& ptree, GGrid &grid, EqnBasePtr &peqn, Time &time, State &utmp, State &ub, State &uf)
{
  GBOOL         bret    = FALSE;
  GString       sinit   = ptree.getValue<GString>("initforce_block");
  PropertyTree  vtree   = ptree.getPropertyTree(sinit);

  if      ( "zero"                        == sinit ) {
    for ( GINT i=0; i<uf.size(); i++ ) {
      if ( uf[i] != NULLPTR ) *uf[i] = 0.0;
    }
  }
//else if ( "initforce_myinit"            == sinit ) {
//  bret = ginitstate::impl_mhyforce      (vtree, grid, time, utmp, ub, uf);
//}
  else                                        {
    assert(FALSE && "Specified state initialization method unknown");
  }

  return bret;

} // end, set_by_direct method


//**********************************************************************************
//**********************************************************************************
// METHOD : set_by_comp
// DESC   : Do init of state components by specifying block name for
//          initstate_block, and initializing each group
//          in the state separately. This method uses the CompDesc data
//          in the EqnBase pointer to locate variable groups.
// ARGS   : ptree  : main property tree
//          grid   : grid object
//          peqn   : ptr to EqnBase 
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          uf     : force to be initialized. May be NULLPTR. 
// RETURNS: TRUE on success; else FALSE
//**********************************************************************************
template<typename EquationType>
GBOOL GInitForceFactory<EquationType>::set_by_comp(const PropertyTree& ptree, GGrid &grid, EqnBasePtr &peqn, Time &time, State &utmp, State &ub, State &uf)
{
  GBOOL           bret    = TRUE;
  GSIZET          ndistcomp, mvar, nvar;
  GString         sblk    = ptree.getValue<GString>("initforce_block");
  GString         sinit;
  GStateCompType *distcomp, *ivar;
  PropertyTree    vtree   = ptree.getPropertyTree(sblk);
  State           comp;
  CompDesc       *icomptype = &peqn.comptype();

  ndistcomp = 0; // # distinct comp types
  distcomp  = NULLPTR; // list of ids for  distinct comp types

  mvar      = 0; // # max of specific comp types
  ivar      = NULLPTR;  // list of ids for specific comp types

  // Get distinct component types:
  icomptype->distinct(distcomp, ndistcomp);

  // Cycle over all types required, get components of that
  // type, and initialize all of them. There should be a member
  // function for each GStateCompType:
  for ( GSIZET j=0; j<ndistcomp && bret; j++ ) {
    
    switch ( distcomp[j] ) {
      
      case GSC_KINETIC:
        sinit = vtree.getValue<GString>("initfv");
        nvar = icomptype->contains(distcomp[j], ivar, mvar);
        comp.resize(nvar);
        for ( GINT i=0; i<nvar; i++ ) comp[i] = uf[ivar[i]];
        bret = doinitfv(vtree, grid, time, utmp, ub, comp);
        break;
      case GSC_MAGNETIC:
        sinit = vtree.getValue<GString>("initfb");
        nvar = icomptype->contains(distcomp[j], ivar, mvar);
        comp.resize(nvar);
        for ( GINT i=0; i<nvar; i++ ) comp[i] = uf[ivar[i]];
        bret = doinitfb(vtree, grid, time, utmp, ub, comp);
        break;
      case GSC_ACTIVE_SCALAR:
        sinit = vtree.getValue<GString>("initfs");
        nvar = icomptype->contains(distcomp[j], ivar, mvar);
        comp.resize(nvar);
        for ( GINT i=0; i<nvar; i++ ) comp[i] = uf[ivar[i]];
        bret = doinitfs(vtree, grid, time, utmp, ub, comp);
        break;
      case GSC_PASSIVE_SCALAR:
        sinit = vtree.getValue<GString>("initfps");
        nvar = icomptype->contains(distcomp[j], ivar, mvar);
        comp.resize(nvar);
        for ( GINT i=0; i<nvar; i++ ) comp[i] = uf[ivar[i]];
        bret = doinitfps(vtree, grid, time, utmp, ub, comp);
        break;
      case GSC_PRESCRIBED:
      case GSC_NONE:
        break;
      default:
        assert(FALSE && "Invalid component type");
    } // end, switch

  } // end, j loop

  if ( distcomp  != NULLPTR ) delete [] distcomp;
  if ( ivar   != NULLPTR ) delete [] ivar ;

  return bret;

} // end, set_by_comp method


//**********************************************************************************
//**********************************************************************************
// METHOD : doinitfv
// DESC   : Do init of kinetic force components. Full list of available
//          kinetic initializations are contained here. Only
//          kinetic components are passed in.
// ARGS   : vtree  : initial condition property tree
//          grid   : grid object
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          uf     : force to be initialized. 
// RETURNS: TRUE on success; else FALSE
//**********************************************************************************
template<typename EquationType>
GBOOL GInitForceFactory<EquationType>::doinitfv(const PropertyTree &vtree, GGrid &grid, Time &time, State &utmp, State &ub, State &uf)
{
  GBOOL           bret    = TRUE;
  GString         sinit = vtree.getValue<GString>("name");

  if      ( "null"   == sinit
       ||   ""             == sinit ) {
    for ( GINT i=0; i<uf.size(); i++ ) *uf[i] = 0.0;
    bret = TRUE;
  }
  else if ( "zer0" == sinit ) {
    for ( GINT i=0; i<uf.size(); i++ ) *uf[i] = 0.0;
  }
//else if ( "random" == sinit ) {
//  bret = ginitfv::random(vtree, grid, time, utmp, ub, uf);
//} 
  else {
    assert(FALSE && "Unknown velocity initialization method");
  }

  return bret;

} // end, doinitfv method


//**********************************************************************************
//**********************************************************************************
// METHOD : doinitfb
// DESC   : Do init of magnetic force components. Full list of available
//          magnetic initializations are contained here. Only
//          magnetic components are passed in.
// ARGS   : vtree  : initial condition property tree
//          grid   : grid object
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          uf     : force to be initialized. 
// RETURNS: TRUE on success; else FALSE
//**********************************************************************************
template<typename EquationType>
GBOOL GInitForceFactory<EquationType>::doinitfb(const PropertyTree &vtree, GGrid &grid, Time &time, State &utmp, State &ub, State &uf)
{
  GBOOL           bret    = FALSE;
  GString         sinit = vtree.getValue<GString>("name");

  if      ( "null"   == sinit
       ||   ""             == sinit ) {
    bret = TRUE;
  }
  else if ( "zero" == sinit ) {
    for ( GINT i=0; i<uf.size(); i++ ) *uf[i] = 0.0;
    bret = TRUE;
  }
//else if ( "random" == sinit ) {
//  bret = ginitfb::random(vtree, grid, time, utmp, ub, uf);
//} 
  else {
    assert(FALSE && "Unknown b-field initialization method");
  }

  return bret;

} // end, doinitfb method



//**********************************************************************************
//**********************************************************************************
// METHOD : doinitfs
// DESC   : Do init of active scalar force components. Full list of available
//          scalar (passive & active) initializations are contained here.
//          Only scalar components are passed in.
// ARGS   : vtree  : initial condition property tree
//          grid   : grid object
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          uf     : force to be initialized. 
// RETURNS: TRUE on success; else FALSE
//**********************************************************************************
template<typename EquationType>
GBOOL GInitForceFactory<EquationType>::doinitfs(const PropertyTree &vtree, GGrid &grid,  Time &time, State &utmp, State &ub, State &uf)
{
  GBOOL           bret    = TRUE;
  GString         sinit = vtree.getValue<GString>("name");

  if      ( "null"   == sinit
       ||   ""             == sinit ) {
    bret = TRUE;
  }
  else if ( "zero" == sinit ) {
    for ( GINT i=0; i<uf.size(); i++ ) *uf[i] = 0.0;
    bret = TRUE;
  }
//else if ( "random" == sinit ) {
//  bret = ginitfs::random(vtree, grid, time, utmp, ub, uf);
//} 
  else {
    assert(FALSE && "Unknown b-field initialization method");
  }

  return bret;
} // end, doinitfs method


//**********************************************************************************
//**********************************************************************************
// METHOD : doinitfps
// DESC   : Do init of passive scalar force components. Full list of available
//          scalar (passive & active) initializations are contained here.
//          Only scalar components are passed in.
// ARGS   : vtree  : initial condition property tree
//          grid   : grid object
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          uf     : state to be initialized. 
// RETURNS: TRUE on success; else FALSE
//**********************************************************************************
template<typename EquationType>
GBOOL GInitForceFactory<EquationType>::doinitfps(const PropertyTree &vtree, GGrid &grid,  Time &time, State &utmp, State &ub, State &uf)
{
  GBOOL           bret    = TRUE;
  GString         sinit = vtree.getValue<GString>("name");

  if      ( "null"   == sinit
       ||   ""             == sinit ) {
    bret = TRUE;
  }
  else if ( "zero" == sinit ) {
    for ( GINT i=0; i<uf.size(); i++ ) *uf[i] = 0.0;
    bret = TRUE;
  }
//else if ( "random" == sinit ) {
//  bret = ginitfps::random(vtree, grid, time, utmp, ub, uf);
//} 
  else {
    assert(FALSE && "Unknown b-field initialization method");
  }

  return bret;
} // end, doinitfps method


