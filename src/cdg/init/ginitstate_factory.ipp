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
GBOOL GInitStateFactory<EquationType>::init(const PropertyTree& ptree, GGrid &grid, EqnBasePtr &peqn, Time &time, State &utmp, State &ub, State &u)
{
  GBOOL         bret    = FALSE;
  GString       stype;  

  // Get type of initialization: by-name or by-block:
  stype = ptree.getValue<GString>("initstate_type","");
  if ( "name"   == stype 
    || ""       == stype ) {
    bret = set_by_name(ptree, grid, peqn, time, utmp, ub, u);
  }
  else if ( "block" == stype ) {
    bret = set_by_blk (ptree, grid, peqn, time, utmp, ub, u);
  }
  else {
    assert(FALSE && "Invalid state initialization type");
  }

  return bret;

} // end, init method


//**********************************************************************************
//**********************************************************************************
// METHOD : set_by_name
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
GBOOL GInitStateFactory<EquationType>::set_by_name(const PropertyTree& ptree, GGrid &grid, EqnBasePtr &peqn, Time &time, State &utmp, State &ub, State &u)
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
    bret = ginitstate::impl_icosgauss  (vtree, eqn_ptr, grid, time, utmp, ub, u);
  }
  else if ( "initstate_boxdirgauss"       == sinit ) {
    bret = ginitstate::impl_boxdirgauss(vtree, eqn_ptr, grid, time, utmp, ub, u);
  }
  else if ( "initstate_boxpergauss"       == sinit ) {
    bret = ginitstate::impl_boxpergauss(vtree, eqn_ptr, grid, time, utmp, ub, u);
  }
  else if ( "initstate_nwave"             == sinit ) {
    bret = ginitstate::impl_nwave      (vtree, eqn_ptr, grid, time, utmp, ub, u);
  }
  else                                        {
    assert(FALSE & "Specified state initialization method unknown");
  }

  return bret;

} // end, set_by_name method


//**********************************************************************************
//**********************************************************************************
// METHOD : set_by_blk
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
// RETURNS: none.
//**********************************************************************************
template<typename EquationType>
GBOOL GInitStateFactory<EquationType>::set_by_blk(const PropertyTree& ptree, GGrid &grid, EqnBasePtr &peqn, Time &time, State &utmp, State &ub, State &u)
{
  GBOOL           bret    = TRUE;
  GSISET          ndistcomp, mvar, nvar;
  GString         sblk    = ptree.getValue<GString>("initstate_block");
  GString         sinit;
  GStateCompType *distcomp, *ivar;
  PropertyTree    vtree     = ptree.getPropertyTree(sblk);
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
        for ( GINT i=0; i<nvar; i++ ) comp[i] = u[ivar[i]];
        bret = doinitv(vtree, grid, peqn, time, utmp, ub, comp);
        break;
      case GSC_MAGNETIC:
        sinit = vtree.getValue<String>("initb");
        nvar = icomptype->contains(distcomp[j], ivar, mvar);
        comp.resize(nvar);
        for ( GINT i=0; i<nvar; i++ ) comp[i] = u[ivar[i]];
        bret = doinitb(vtree, grid, peqn, time, utmp, ub, comp);
        break;
      case GSC_ACTIVE_SCALAR:
        sinit = vtree.getValue<String>("inits");
        nvar = icomptype->contains(distcomp[j], ivar, mvar);
        comp.resize(nvar);
        for ( GINT i=0; i<nvar; i++ ) comp[i] = u[ivar[i]];
        bret = doinits(vtree, grid, peqn, time, utmp, ub, comp);
        break;
      case GSC_PASSIVE_SCALAR:
        sinit = vtree.getValue<String>("initp");
        nvar = icomptype->contains(distcomp[j], ivar, mvar);
        comp.resize(nvar);
        for ( GINT i=0; i<nvar; i++ ) comp[i] = u[ivar[i]];
        bret = doinits(vtree, grid, peqn, time, utmp, ub, comp);
        break;
      case GSC_PRESCRIBED:
        sinit = vtree.getValue<String>("initc");
        nvar = icomptype->contains(distcomp[j], ivar, mvar);
        comp.resize(nvar);
        for ( GINT i=0; i<nvar; i++ ) comp[i] = u[ivar[i]];
        bret = doinitc(vtree, grid, peqn, time, utmp, ub, comp);
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

} // end, set_by_blk method


//**********************************************************************************
//**********************************************************************************
// METHOD : doinitv
// DESC   : Do init of kinetic components. Full list of available
//          kinetic initializations are contained here. Only
//          kinetic components are passed in.
// ARGS   : vtree  : initial condition property tree
//          grid   : grid object
//          peqn   : pointer to EqnBase
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          u      : state to be initialized. 
// RETURNS: none.
//**********************************************************************************
template<typename EquationType>
GBOOL GInitStateFactory<EquationType>::doinitv(const PropertyTree &vtree, GGrid &grid, EqnBasePtr &peqn,  Time &time, State &utmp, State &ub, State &u)
{
  GBOOL           bret    = TRUE;
  GString         sinit = vtree.getValue<GString>("name");

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
//  bret = ginitv::random(vtree, grid, peqn, time, utmp, ub, u);
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
// RETURNS: none.
//**********************************************************************************
template<typename EquationType>
GBOOL GInitStateFactory<EquationType>::doinitb(const PropertyTree &vtree, GGrid &grid, EqnBasePtr &peqn, Time &time, State &utmp, State &ub, State &u)
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
//  bret = ginitb::random(vtree, grid, peqn, time, utmp, ub, u);
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
//          peqn   : EqnBase pointer
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          u      : state to be initialized. 
// RETURNS: none.
//**********************************************************************************
template<typename EquationType>
GBOOL GInitStateFactory<EquationType>::doinits(const PropteryTree &vtree, GGrid &grid,  EqnBasePtr &peqn, Time &time, State &utmp, State &ub, State &u)
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
//  bret = ginits::random(vtree, grid, peqn, time, utmp, ub, u);
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
//          peqn   : EqnBase pointer
//          time   : initialization time
//          utmp   : tmp arrays
//          ub     : boundary state (also initialized here)
//          u      : state to be initialized. 
// RETURNS: none.
//**********************************************************************************
template<typename EquationType>
GBOOL GInitStateFactory<EquationType>::doinitps(const PropteryTree &vtree, GGrid &grid,  EqnBasePtr &peqn, Time &time, State &utmp, State &ub, State &u)
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
//  bret = ginits::random(vtree, grid, peqn, time, utmp, ub, u);
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
// RETURNS: none.
//**********************************************************************************
template<typename EquationType>
GBOOL GInitStateFactory<EquationType>::doinitc(const PropteryTree &vtree, GGrid &grid, EqnBaseePtr &peqn, Time &time, State &utmp, State &ub, State &u)
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
//  bret = ginitc::random(vtree, grid, peqn, time, utmp, ub, u);
//} 
  else {
    assert(FALSE && "Unknown b-field initialization method");
  }

  return bret;
} // end, doinitc method


