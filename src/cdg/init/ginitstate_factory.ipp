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
GBOOL GInitStateFactory<EquationType>::init(const PropertyTree& ptree, GGrid &grid, EqnBasePtr &peqn, Time &time, State &utmp, State &ub, State &u)
{
  GBOOL         bret    = FALSE;
  GString       stype;  

  // Get type of initialization: by-name or by-block:
  stype = ptree.getValue<GString>("initstate_type","");
  if ( "direct"   == stype 
    || ""         == stype ) {
    bret = set_by_direct(ptree, grid, peqn, time, utmp, ub, u);
  }
  else if ( "component" == stype ) {
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
// RETURNS: TRUE on success; else FALSE
//**********************************************************************************
template<typename EquationType>
GBOOL GInitStateFactory<EquationType>::set_by_direct(const PropertyTree& ptree, GGrid &grid, EqnBasePtr &peqn, Time &time, State &utmp, State &ub, State &u)
{
  GBOOL         bret    = FALSE;
  GString       sinit   = ptree.getValue<GString>("inits_block");
  PropertyTree  vtree   = ptree.getPropertyTree(sinit);

  if      ( "zero"                        == sinit ) {
    for ( GINT i=0; i<u.size(); i++ ) {
      if ( u[i] != NULLPTR ) *u[i] = 0.0;
    } 
  }
  else if ( "initstate_icosgaussburgers"  == sinit ) {
    bret = ginitstate::impl_icosgauss       (vtree, grid, time, utmp, ub, u);
  }
  else if ( "initstate_boxdirgauss"       == sinit ) {
    bret = ginitstate::impl_boxdirgauss     (vtree, grid, time, utmp, ub, u);
  }
  else if ( "initstate_boxpergauss"       == sinit ) {
    bret = ginitstate::impl_boxpergauss     (vtree, grid, time, utmp, ub, u);
  }
  else if ( "initstate_nwave"             == sinit ) {
    bret = ginitstate::impl_boxnwaveburgers (vtree, grid, time, utmp, ub, u);
  }
  else                                        {
    assert(FALSE && "Specified state initialization method unknown");
  }

  return bret;

} // end, set_by_direct method


//**********************************************************************************
//**********************************************************************************
// METHOD : set_by_comp
// DESC   : Do init of state components by specifying individual
//          variable group types within JSON block
//          in the state separately. This method uses the CompDesc data
//          in the EqnBase pointer to locate variable groups.
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
GBOOL GInitStateFactory<EquationType>::set_by_comp(const PropertyTree& ptree, GGrid &grid, EqnBasePtr &peqn, Time &time, State &utmp, State &ub, State &u)
{
  GBOOL           bret    = TRUE;
  GSIZET          ndist, *ivar, mvar, nvar;
  GSIZET         *idist;
  GStateCompType  itype;
  GString         scomp    = ptree.getValue<GString>("initstate_block");
  GString         sinit;
  PropertyTree    vtree     = ptree.getPropertyTree(scomp);
  State           comp;
  CompDesc       *icomptype = &peqn->comptype();

  ndist = 0; // # distinct comp types
  idist = NULLPTR; // list of ids for  distinct comp types

  mvar      = 0; // # max of specific comp types
  ivar      = NULLPTR;  // list of ids for specific comp types

  // Get distinct component types:
  icomptype->distinct(idist, ndist);

  // Cycle over all types required, get components of that
  // type, and initialize all of them. There should be a member
  // function for each GStateCompType:
  for ( GSIZET j=0; j<ndist && bret; j++ ) {
    
    itype = (*icomptype)[idist[j]]; 
    switch ( itype ) {

      case GSC_KINETIC:
        sinit = vtree.getValue<GString>("initv");
        nvar  = icomptype->contains(itype, ivar, mvar);
        comp.resize(nvar);
        for ( GINT i=0; i<nvar; i++ ) comp[i] = u[ivar[i]];
        bret = doinitv(vtree, grid, time, utmp, ub, comp);
        break;
      case GSC_MAGNETIC:
        sinit = vtree.getValue<GString>("initb");
        nvar  = icomptype->contains(itype, ivar, mvar);
        comp.resize(nvar);
        for ( GINT i=0; i<nvar; i++ ) comp[i] = u[ivar[i]];
        bret = doinitb(vtree, grid, time, utmp, ub, comp);
        break;
      case GSC_ACTIVE_SCALAR:
        sinit = vtree.getValue<GString>("inits");
        nvar  = icomptype->contains(itype, ivar, mvar);
        comp.resize(nvar);
        for ( GINT i=0; i<nvar; i++ ) comp[i] = u[ivar[i]];
        bret = doinits(vtree, grid, time, utmp, ub, comp);
        break;
      case GSC_PASSIVE_SCALAR:
        sinit = vtree.getValue<GString>("initp");
        nvar  = icomptype->contains(itype, ivar, mvar);
        comp.resize(nvar);
        for ( GINT i=0; i<nvar; i++ ) comp[i] = u[ivar[i]];
        bret = doinits(vtree, grid, time, utmp, ub, comp);
        break;
      case GSC_PRESCRIBED:
        sinit = vtree.getValue<GString>("initc");
        nvar  = icomptype->contains(itype, ivar, mvar);
        comp.resize(nvar);
        for ( GINT i=0; i<nvar; i++ ) comp[i] = u[ivar[i]];
        bret = doinitc(vtree, grid, time, utmp, ub, comp);
        break;
      case GSC_NONE:
        break;
      default:
        assert(FALSE && "Invalid component type");
    } // end, switch

  } // end, j loop

  if ( idist  != NULLPTR ) delete [] idist;
  if ( ivar   != NULLPTR ) delete [] ivar ;

  return bret;

} // end, set_by_comp method


//**********************************************************************************
//**********************************************************************************
// METHOD : doinitv
// DESC   : Do init of kinetic components. Full list of available
//          kinetic initializations are contained here. Only
//          kinetic components are passed in.
// ARGS   : vtree  : initial condition property tree
//          grid   : grid object
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          u      : state to be initialized. 
// RETURNS: TRUE on success; else FALSE
//**********************************************************************************
template<typename EquationType>
GBOOL GInitStateFactory<EquationType>::doinitv(const PropertyTree &vtree, GGrid &grid, Time &time, State &utmp, State &ub, State &u)
{
  GBOOL           bret    = TRUE;
  GString         sinit = vtree.getValue<GString>("name");

  if      ( "null"   == sinit
       ||   ""       == sinit ) {
    for ( GINT i=0; i<u.size(); i++ ) *u[i] = 0.0;
    bret = TRUE;
  }
  else if ( "zer0" == sinit ) {
    for ( GINT i=0; i<u.size(); i++ ) *u[i] = 0.0;
  }
//else if ( "random" == sinit ) {
//  bret = ginitv::random(vtree, grid, time, utmp, ub, u);
//} 
  else {
    assert(FALSE && "Unknown velocity initialization method");
  }

  return bret;

} // end, doinitv method


//**********************************************************************************
//**********************************************************************************
// METHOD : doinitb
// DESC   : Do init of magnetic components. Full list of available
//          magnetic initializations are contained here. Only
//          magnetic components are passed in.
// ARGS   : vtree  : initial condition property tree
//          grid   : grid object
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          u      : state to be initialized. 
// RETURNS: TRUE on success; else FALSE
//**********************************************************************************
template<typename EquationType>
GBOOL GInitStateFactory<EquationType>::doinitb(const PropertyTree &vtree, GGrid &grid, Time &time, State &utmp, State &ub, State &u)
{
  GBOOL           bret    = FALSE;
  GString         sinit = vtree.getValue<GString>("name");

  if      ( "null"   == sinit
       ||   ""             == sinit ) {
    bret = TRUE;
  }
  else if ( "zero" == sinit ) {
    for ( GINT i=0; i<u.size(); i++ ) *u[i] = 0.0;
    bret = TRUE;
  }
//else if ( "random" == sinit ) {
//  bret = ginitb::random(vtree, grid, time, utmp, ub, u);
//} 
  else {
    assert(FALSE && "Unknown b-field initialization method");
  }

  return bret;

} // end, doinitb method



//**********************************************************************************
//**********************************************************************************
// METHOD : doinits
// DESC   : Do init of active scalar components. Full list of available
//          scalar (passive & active) initializations are contained here.
//          Only scalar components are passed in.
// ARGS   : vtree  : initial condition property tree
//          grid   : grid object
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          u      : state to be initialized. 
// RETURNS: TRUE on success; else FALSE
//**********************************************************************************
template<typename EquationType>
GBOOL GInitStateFactory<EquationType>::doinits(const PropertyTree &vtree, GGrid &grid,  Time &time, State &utmp, State &ub, State &u)
{
  GBOOL           bret    = TRUE;
  GString         sinit = vtree.getValue<GString>("name");

  if      ( "null"   == sinit
       ||   ""             == sinit ) {
    bret = TRUE;
  }
  else if ( "zero" == sinit ) {
    for ( GINT i=0; i<u.size(); i++ ) *u[i] = 0.0;
    bret = TRUE;
  }
//else if ( "random" == sinit ) {
//  bret = ginits::random(vtree, grid, time, utmp, ub, u);
//} 
  else {
    assert(FALSE && "Unknown b-field initialization method");
  }

  return bret;
} // end, doinits method


//**********************************************************************************
//**********************************************************************************
// METHOD : doinitps
// DESC   : Do init of passive scalar components. Full list of available
//          scalar (passive & active) initializations are contained here.
//          Only scalar components are passed in.
// ARGS   : vtree  : initial condition property tree
//          grid   : grid object
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          u      : state to be initialized. 
// RETURNS: TRUE on success; else FALSE
//**********************************************************************************
template<typename EquationType>
GBOOL GInitStateFactory<EquationType>::doinitps(const PropertyTree &vtree, GGrid &grid,  Time &time, State &utmp, State &ub, State &u)
{
  GBOOL           bret = FALSE;
  GString         sinit = vtree.getValue<GString>("name");

  if      ( "null"   == sinit
       ||   ""             == sinit ) {
    bret = TRUE;
  }
  else if ( "zero" == sinit ) {
    for ( GINT i=0; i<u.size(); i++ ) *u[i] = 0.0;
    bret = TRUE;
  }
//else if ( "random" == sinit ) {
//  bret = ginits::random(vtree, grid, time, utmp, ub, u);
//} 
  else {
    assert(FALSE && "Unknown b-field initialization method");
  }

  return bret;
} // end, doinitps method



//**********************************************************************************
//**********************************************************************************
// METHOD : doinitc
// DESC   : Do init of prescribed components. Full list of available
//          prescribed initializations are contained here. Only
//          prescribed components are passed in.
// ARGS   : vtree  : initial condition property tree
//          grid   : grid object
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          u      : state to be initialized. 
// RETURNS: TRUE on success; else FALSE
//**********************************************************************************
template<typename EquationType>
GBOOL GInitStateFactory<EquationType>::doinitc(const PropertyTree &vtree, GGrid &grid, Time &time, State &utmp, State &ub, State &u)
{
  GBOOL           bret  = FALSE;
  GString         sinit = vtree.getValue<GString>("name");

  if      ( "null"   == sinit
       ||   ""             == sinit ) {
    bret = TRUE;
  }
  else if ( "zero" == sinit ) {
    for ( GINT i=0; i<u.size(); i++ ) *u[i] = 0.0;
    bret = TRUE;
  }
//else if ( "random" == sinit ) {
//  bret = ginitc::random(vtree, grid, time, utmp, ub, u);
//} 
  else {
    assert(FALSE && "Unknown b-field initialization method");
  }

  return bret;
} // end, doinitc method

