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
GBOOL GInitStateFactory<EquationType>::set_direct(const geoflow::tbox::PropertyTree& ptree, GGrid &grid, EqnBasePtr &peqn, Time &time, State &utmp, State &ub, State &u)
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
GBOOL GInitStateFactory<EquationType>::set_by_grp(const geoflow::tbox::PropertyTree& ptree, GGrid &grid, EqnBasePtr &peqn, Time &time, State &utmp, State &ub, State &u)
{
  GBOOL           bret    = TRUE;
  GSISET          ndistcomp, mvar, nvar;
  GString         sblk    = ptree.getValue<GString>("inits_block");
  GString         sinit;
  GStateCompType *distcomp, *ivar;
  PropertyTree    vtree   = ptree.getPropertyTree(sblk);
  State           comp;
  CompDesc       *icomptype = &peqn.compdesc();

  ndistcomp = 0; // # distinct comp types
  distcomp  = NULLPTR; // list of ids for  distinct comp types

  mvar      = 0; // # max of specific comp types
  ivar      = NULLPTR;  // list of ids for specific comp types

  // Get distinct component types:
  icomptype->distinct(distcomp, ndistcomp);

  // Cycle over all types required, get components of that
  // type, and initialize all of them. There should be a member
  // function for each GStateCompType:
  for ( GSIZET j=0; j<ncomp && bret; j++ ) {
    
    switch ( distcomp[j] ) {
      
      case GSC_KINETIC:
        sinit = vtree.getValue<String>("initv");
        nvar = icomptype->contains(distcomp[j], ivar, mvar);
        comp.resize(nvar);
        for ( GINT i=0; i<nvar; i++ ) comp[i] = ivar[i];
        bret = doinitv(sinit, grid, time, utmp, ub, comp);
        break;
      case GSC_MAGNETIC:
        sinit = vtree.getValue<String>("initb");
        nvar = icomptype->contains(distcomp[j], ivar, mvar);
        comp.resize(nvar);
        for ( GINT i=0; i<nvar; i++ ) comp[i] = ivar[i];
        bret = doinitb(sinit, grid, time, utmp, ub, comp);
        break;
      case GSC_ACTIVE_SCALAR:
        sinit = vtree.getValue<String>("inits");
        nvar = icomptype->contains(distcomp[j], ivar, mvar);
        comp.resize(nvar);
        for ( GINT i=0; i<nvar; i++ ) comp[i] = ivar[i];
        bret = doinits(sinit, grid, time, utmp, ub, comp);
        break;
      case GSC_PASSIVE_SCALAR:
        sinit = vtree.getValue<String>("initp");
        nvar = icomptype->contains(distcomp[j], ivar, mvar);
        comp.resize(nvar);
        for ( GINT i=0; i<nvar; i++ ) comp[i] = ivar[i];
        bret = doinits(sinit, grid, time, utmp, ub, comp);
        break;
      case GSC_PRESCRIBED:
        sinit = vtree.getValue<String>("initc");
        nvar = icomptype->contains(distcomp[j], ivar, mvar);
        comp.resize(nvar);
        for ( GINT i=0; i<nvar; i++ ) comp[i] = ivar[i];
        bret = doinitc(sinit, grid, time, utmp, ub, comp);
        break;
      case GSC_NONE:
        break;
      default:
        assert(FALSE && "Invalid component type");
    } // end, switch

  } // end, j loop

  if ( distcomp  != NULLPTR ) delete [] distcomp;
  if ( ivar   != NULLPTR ) delete [] ivar ;

  return bret;

} // end, set_by_grp method


//**********************************************************************************
//**********************************************************************************
// METHOD : doinitv
// DESC   : Do init of kinetic components. Full list of available
//          kinetic initializations are contained here. Only
//          kinetic components are passed in.
// ARGS   : ptree  : main property tree
//          grid   : grid object
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          u      : state to be initialized. 
// RETURNS: none.
//**********************************************************************************
template<typename EquationType>
GBOOL GInitStateFactory<EquationType>::doinitv(GString &sinit, GGrid &grid, Time &time, State &utmp, State &ub, State &u)
{
  GBOOL           bret    = TRUE;

  if      ( "null"   == sinit
       ||   ""             == sinit ) {
    for ( GINT i=0; i<u.size(); i++ ) *u[i] = 0.0;
    bret = TRUE;
  }
  else if ( "zer0" == sinit ) {
    bret = ginitv::random(grid, time, utmp, ub, u);
    for ( GINT i=0; i<u.size(); i++ ) *u[i] = 0.0;
  }
//else if ( "random" == sinit ) {
//  bret = ginitv::random(grid, time, utmp, ub, u);
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
// ARGS   : ptree  : main property tree
//          grid   : grid object
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          u      : state to be initialized. 
// RETURNS: none.
//**********************************************************************************
template<typename EquationType>
GBOOL GInitStateFactory<EquationType>::doinitb(GString &sinit, GGrid &grid, Time &time, State &utmp, State &ub, State &u)
{
  GBOOL           bret    = FALSE;

  if      ( "null"   == sinit
       ||   ""             == sinit ) {
    bret = TRUE;
  }
  else if ( "zero" == sinit ) {
    for ( GINT i=0; i<u.size(); i++ ) *u[i] = 0.0;
    bret = TRUE;
  }
//else if ( "random" == sinit ) {
//  bret = ginitb::random(grid, time, utmp, ub, u);
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
// ARGS   : ptree  : main property tree
//          grid   : grid object
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          u      : state to be initialized. 
// RETURNS: none.
//**********************************************************************************
template<typename EquationType>
GBOOL GInitStateFactory<EquationType>::doinits(GString &sinit, GGrid &grid, Time &time, State &utmp, State &ub, State &u)
{
  GBOOL           bret    = TRUE;

  if      ( "null"   == sinit
       ||   ""             == sinit ) {
    bret = TRUE;
  }
  else if ( "zero" == sinit ) {
    for ( GINT i=0; i<u.size(); i++ ) *u[i] = 0.0;
    bret = TRUE;
  }
//else if ( "random" == sinit ) {
//  bret = ginits::random(grid, time, utmp, ub, u);
//} 
  else {
    assert(FALSE && "Unknown b-field initialization method");
  }

  return bret;
} // end, doinits method



//**********************************************************************************
//**********************************************************************************
// METHOD : doinitc
// DESC   : Do init of prescribed components. Full list of available
//          prescribed initializations are contained here. Only
//          prescribed components are passed in.
// ARGS   : ptree  : main property tree
//          grid   : grid object
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          u      : state to be initialized. 
// RETURNS: none.
//**********************************************************************************
template<typename EquationType>
GBOOL GInitStateFactory<EquationType>::doinitc(GString &sinit, GGrid &grid, Time &time, State &utmp, State &ub, State &u)
{
  GBOOL           bret    = TRUE;

  if      ( "null"   == sinit
       ||   ""             == sinit ) {
    bret = TRUE;
  }
  else if ( "zero" == sinit ) {
    for ( GINT i=0; i<u.size(); i++ ) *u[i] = 0.0;
    bret = TRUE;
  }
//else if ( "random" == sinit ) {
//  bret = ginitc::random(grid, time, utmp, ub, u);
//} 
  else {
    assert(FALSE && "Unknown b-field initialization method");
  }

  return bret;
} // end, doinits method




