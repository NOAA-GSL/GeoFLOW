//==================================================================================
// Module       : ginitbdy_factory
// Date         : 7/11/19 (DLR)
// Description  : GeoFLOW state initialization factory. Note: This should be 
//                called before GInitStateFactory is called, as the state, u, 
//                is used here for temp space.
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================


//**********************************************************************************
//**********************************************************************************
// METHOD : init
// DESC   : Make call to to do bdy initialization
// ARGS   : ptree  : main property tree
//          time   : initialization time
//          utmp   : tmp arrays
//          u      : state to be initialized. 
//          ub     : boundary state 
// RETURNS: none.
//**********************************************************************************
template<typename EquationType>
GBOOL GInitBdyFactory<EquationType>::init(const geoflow::tbox::PropertyTree& ptree, GGrid &grid, Time &time, State &utmp, State &u, State &ub)
{
  GBOOL         bret=FALSE;
  GBOOL         use_inits; // use state init method to set bdy?
  GFTYPE        tt;
  GString       sgrid, supdate;
  PropertyTree  gtree;  

  sgrid     = ptree.getValue<GString>("grid_type");
  gtree     = ptree.getPropertyTree(sgrid);
  supdate   = gtree.getValue<GStrig>("bdy_init_method","none");
  use_inits = gtree.getValue<GStrig>("use_state_init",FALSE);
  tt        = time;

  if ( "initb_none" == sinitb
    || "none"       == sinitb 
    || ""           == sinitb ) {
    bret = TRUE;
  }
  else if ( use_inits ) { // initialize from state initialization
    bret = GInitStateFactory<EquationType>::GInitStateFactory::init(ptree, grid, tt, utmp, ub, u);
    if ( bret ) {
      bret = setbdy_from_state(ptree, grid, tt, utmp, u, ub);
    }
  }
  else if ( "mybdyinit" == sinit ) { // initialize from user-specified fcn
    bret = ginitbdy::impl_mybdyinit  (ptree, grid, tt, utmp, u, ub);
  }
  else                                        {
    assert(FALSE && "Specified bdy initialization method unknown");
  }

  return bret;
} // end, init method init


//**********************************************************************************
//**********************************************************************************
// METHOD : set_bdy_from_state
// DESC   : use state var, u, to set bdy, ub
// ARGS   : ptree  : main property tree
//          time   : initialization time
//          utmp   : tmp arrays
//          u      : state to be initialized. 
//          ub     : boundary state 
// RETURNS: none.
//**********************************************************************************
template<typename EquationType>
void GInitBdyFactory<EquationType>::set_bdy_from_state(const geoflow::tbox::PropertyTree& ptree, GGrid &grid, Time &time, State &utmp, State &u, State &ub)
{
  GBOOL         bret=FALSE;
  GBOOL         use_inits; // use state init method to set bdy?
  GFTYPE        tt;
  GString       sgrid, supdate;

  GTVector<GTVector<GSIZET>> *igbdy = &grid.igbdy_binned();

  // Set from State vector, u and others that we _can_ set:
  for ( auto k=0; k<u.size(); k++ ) { 
    for ( auto j=0; j<(*igbdy)[GBDY_DIRICHLET].size()
       && ub[k] != NULLPTR; j++ ) {
      (*ub[k])[j] = (*uu[k])[(*igbdy)[GBDY_DIRICHLET][j]];
    }
    for ( auto j=0; j<(*igbdy)[GBDY_INFLOWT].size()
       && ub[k] != NULLPTR; j++ ) {
      (*ub[k])[j] = (*uu[k])[(*igbdy)[GBDY_INFLOWT][j]];
    }
    for ( auto j=0; j<(*igbdy)[GBDY_NOSLIP].size()
       && ub[k] != NULLPTR; j++ ) {
      (*ub[k])[j] = 0.0;
    }
  }

} // end, set_bdy_from_state
