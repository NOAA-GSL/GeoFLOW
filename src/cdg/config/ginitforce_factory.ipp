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
//          eqn    : equation implementation
//          grid   : grid object
//          time   : initialization time
//          utmp   : tmp arrays
//          u      : state to be initialized. 
//          uf     : resultsing force
// RETURNS: TRUE on success; else FALSE
//**********************************************************************************
template<typename Types>
GBOOL GInitForceFactory<Types>::init(const PropertyTree& ptree, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u, State &uf)
{
  GBOOL         bret    = FALSE;
  GBOOL         bforced = ptree.getValue<GBOOL>("use_forcing");
  GString       stype;  
  PropertyTree  vtree;

  if ( !bforced ) return TRUE;

  // Get type of initialization: direct or by-var:
  stype = ptree.getValue<GString>("initforce_type","");
  if ( "direct"   == stype 
    || ""         == stype ) {
    bret = set_by_direct(ptree, eqn, grid, time, utmp, u, uf);
  }
  else if ( "component" == stype ) {
    bret = set_by_comp  (ptree, eqn, grid, time, utmp, u, uf);
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
//          eqn    : equation implementation
//          grid   : grid object
//          time   : initialization time
//          utmp   : tmp arrays
//          u      : state to be initialized. 
//          uf     : resulting forces
// RETURNS: TRUE on success; else FALSE
//**********************************************************************************
template<typename Types>
GBOOL GInitForceFactory<Types>::set_by_direct(const PropertyTree& ptree, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u, State &uf)
{
  GBOOL         bret    = FALSE;
  GString       sinit   = ptree.getValue<GString>("initforce_block");
  PropertyTree  vtree   = ptree.getPropertyTree(sinit);

  if      ( "zero"                        == sinit ) {
    for ( GINT i=0; i<uf.size(); i++ ) {
      if ( uf[i] != NULLPTR ) *uf[i] = 0.0;
    }
  }
  else if ( "initforce_rand"            == sinit ) {
    bret = ginitforce<Types>::impl_rand        (ptree, sinit, eqn, grid, time, utmp, u, uf);
  }
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
//          eqn    : equation implementation
//          grid   : grid object
//          time   : initialization time
//          utmp   : tmp arrays
//          u      : state array
//          uf     : force to be initialized. May be NULLPTR. 
// RETURNS: TRUE on success; else FALSE
//**********************************************************************************
template<typename Types>
GBOOL GInitForceFactory<Types>::set_by_comp(const PropertyTree& ptree, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u, State &uf)
{

#if 0
  GBOOL           bret    = TRUE;
  GSIZET          mvar, ndist, nvar;
  GSIZET         *idist, *ivar, *pisz;
  GString         sblk    = ptree.getValue<GString>("initforce_block");
  GString         sinit;
  GStateCompType  itype;
  PropertyTree    vtree   = ptree.getPropertyTree(sblk);
  State           comp;

  GStateCompType *pct;
  GTVector<GSIZET>   itmp;

  itmp .resize(icomptype->size());
  pisz = itmp.data();


  ndist = 0; // # distinct comp types
  idist = NULLPTR; // list of ids for  distinct comp types

  mvar  = 0; // # max of specific comp types
  ivar  = NULLPTR;  // list of ids for specific comp types


  // Get distinct component types:
  icomptype->distinct(idist, ndist, pct, pisz);

  // Cycle over all types required, get components of that
  // type, and initialize all of them. There should be a member
  // function for each GStateCompType:
  for ( GSIZET j=0; j<ndist && bret; j++ ) {

    itype = (*icomptype)[idist[j]];
    switch ( itype ) {
      
      case GSC_KINETIC:
        sinit = vtree.getValue<GString>("initfv");
        nvar  = icomptype->contains(itype, ivar, mvar);
        comp.resize(nvar);
        for ( GINT i=0; i<nvar; i++ ) comp[i] = uf[ivar[i]];
        bret = doinitfv(ptree, sinit, eqn, grid, time, utmp, comp);
        break;
      case GSC_MAGNETIC:
        sinit = vtree.getValue<GString>("initfb");
        nvar  = icomptype->contains(itype, ivar, mvar);
        comp.resize(nvar);
        for ( GINT i=0; i<nvar; i++ ) comp[i] = uf[ivar[i]];
        bret = doinitfb(ptree, sinit, eqn, grid, time, utmp, comp);
        break;
      case GSC_TEMPERATURE:
        sinit = vtree.getValue<GString>("initftemp");
        nvar  = icomptype->contains(itype, ivar, mvar);
        comp.resize(nvar);
        for ( GINT i=0; i<nvar; i++ ) comp[i] = uf[ivar[i]];
        bret = doinitftemp(ptree, sinit, eqn, grid, time, utmp, comp);
        break;
      case GSC_PRESCRIBED:
      case GSC_NONE:
        break;
      default:
        assert(FALSE && "Invalid component type");
    } // end, switch

  } // end, j loop

  if ( idist  != NULLPTR ) delete [] idist;
  if ( ivar   != NULLPTR ) delete [] ivar ;

  return bret;
#endif

  return FALSE;

} // end, set_by_comp method


//**********************************************************************************
//**********************************************************************************
// METHOD : doinitfv
// DESC   : Do init of kinetic force components. Full list of available
//          kinetic initializations are contained here. Only
//          kinetic components are passed in.
// ARGS   : ptree  : main property tree
//          sconfig: ptree block name containing variable config
//          eqn    : equation implementation
//          grid   : grid object
//          time   : initialization time
//          utmp   : tmp arrays
//          uf     : force to be initialized. 
// RETURNS: TRUE on success; else FALSE
//**********************************************************************************
template<typename Types>
GBOOL GInitForceFactory<Types>::doinitfv(const PropertyTree &ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u, State &uf)
{
  GBOOL             bret    = TRUE;
  GString           sinit;
  GridIcos         *icos;
  GridBox          *box;
  PropertyTree      vtree = ptree.getPropertyTree(sconfig);

  sinit = vtree.getValue<GString>("name");

  icos = dynamic_cast<GridIcos*>(&grid);
  box  = dynamic_cast<GridBox*> (&grid);
  if      ( "null"   == sinit
       ||   ""             == sinit ) {     // set to 0
    for ( GINT i=0; i<uf.size(); i++ ) *uf[i] = 0.0;
    bret = TRUE;
  }
  else if ( "zero" == sinit ) {            // set to 0
    for ( GINT i=0; i<uf.size(); i++ ) *uf[i] = 0.0;
  }
  else if ( "abc" == sinit ) {             // ABC forcing

    if       ( icos != NULLPTR ) {
      bret = ginitfv<Types>::impl_abc_icos(ptree, sconfig, eqn, grid, time, utmp, u, uf);
    }
    else {
      bret = ginitfv<Types>::impl_abc_box (ptree, sconfig, eqn, grid, time, utmp, u, uf);
    }

  } 
  else if ( "random" == sinit ) {
    bret = ginitfv<Types>::impl_rand   (ptree, sinit, eqn, grid, time, utmp, u, uf);
  } 
  else {
    assert(FALSE && "Unknown velocity force initialization method");
  }

  return bret;

} // end, doinitfv method


//**********************************************************************************
//**********************************************************************************
// METHOD : doinitftemp
// DESC   : Do init of temperature force components. Full list of available
//          initializations are contained here.
//          Only temperature components are passed in.
// ARGS   : ptree  : initial condition property tree
//          sconfig: ptree block name containing variable config
//          eqn    : equation implementation
//          grid   : grid object
//          time   : initialization time
//          utmp   : tmp arrays
//          u      : state arrays
//          uf     : force to be initialized. 
// RETURNS: TRUE on success; else FALSE
//**********************************************************************************
template<typename Types>
GBOOL GInitForceFactory<Types>::doinitftemp(const PropertyTree &ptree, GString &sconfig, EqnBasePtr &eqn, Grid &grid, Time &time, State &utmp, State &u, State &uf)
{
  GBOOL           bret    = TRUE;
  GString         sinit;
  PropertyTree    vtree = ptree.getPropertyTree(sconfig);

  sinit = vtree.getValue<GString>("name");

  if      ( "null"   == sinit
       ||   ""             == sinit ) {
    bret = TRUE;
  }
  else if ( "zero" == sinit ) {
    for ( GINT i=0; i<uf.size(); i++ ) *uf[i] = 0.0;
    bret = TRUE;
  }
  else if ( "random" == sinit ) {
    bret = ginitfs<Types>::impl_rand(ptree, sinit, eqn, grid, time, utmp, u, uf);
  } 
  else {
    assert(FALSE && "Unknown scalar force  initialization method");
  }

  return bret;
} // end, doinitftemp method


