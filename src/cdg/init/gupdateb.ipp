
namespace gupdateb {




//**********************************************************************************
//**********************************************************************************
// METHOD : impl_bystateinit
// DESC   : Update bdy using state initialization method 
//          specified by ptree
// ARGS   : ptree: main prop tree
//          grid : grid
//          t    : time
//          utmp : tmp arrays; require at least ub.size arrays
//          u    : current state
//          ub   : bdy vectors (one for each state element, unless NULL)
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_bystateinit(const PropteryTree &ptree, GGrid &grid, Time &time, State &utmp, const State &u, State &ub)
{
  Time             tt = t;
  State            uu(u.size());
  GString          serr = "impl_bystateinit: ";

  Time  tt = t;

  // Use tmp from back end, so that 'init' isn't 
  // disturbed. NOTE: this could still be a problem!
  assert(utmp.size() >= 2*u.size() && "Insufficient temp space");
  for ( auto j=0; j<u.size(); j++ ) {
    uu[j] = utmp[utmp.size()-1-j];
  }
  GInitSFactory::init(ptree, grid, tt, utmp, ub, uu);

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

  return TRUE;
} // end, impl_bystateinit


//**********************************************************************************
//**********************************************************************************
// METHOD : impl_simple_outflow
// DESC   : Update bdy using simple outflow. Applies
//          only to GBDY_OUTFLOW bdy types. Called each time step.
// ARGS   : ptree: main prop tree
//          grid : grid
//          t    : time
//          utmp : tmp arrays; require at least ub.size arrays
//          u    : current state, unchanged here
//          ub   : bdy vectors (one for each state element, unless NULL)
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_simple_outflow(const PropteryTree &ptree, GGrid &grid, Time &time, State &utmp, State &u, State &ub)
{
  Time             tt = t;
  State            uu(u.size());
  GString          serr = "impl_simple_outflow: ";

  Time  tt = t;

  GTVector<GTVector<GSIZET>> *igbdy = &grid.igbdy();

  // Set from State vector, u:
  for ( auto k=0; k<u.size(); k++ ) { 
    for ( auto j=0; j<(*igbdy)[GBDY_OUTFLOW].size() 
       && ub[k] != NULLPTR; j++ ) {
      (*ub[k])[j] = (*u[k])[(*igbdy)[GBDY_OUTFLOW][j]];
    }
  }

  return TRUE;

} // end, impl_simple_outflow


//**********************************************************************************
//**********************************************************************************
// METHOD : impl_noslip
// DESC   : Update bdy to be no-slip. Applies
//          only to GBDY_NOSLIP bdy types. Called each time step.
// ARGS   : ptree: main prop tree
//          grid : grid
//          t    : time
//          utmp : tmp arrays; require at least ub.size arrays
//          u    : current state, unchanged here
//          ub   : bdy vectors (one for each state element, unless NULL)
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_noslip(const PropteryTree &ptree, GGrid &grid, Time &time, State &utmp, State &u, State &ub)
{
  Time             tt = t;
  State            uu(u.size());
  GString          serr = "impl_noslip: ";

  Time  tt = t;

  GTVector<GTVector<GSIZET>> *igbdy = &grid.igbdy();

  // Set from State vector, u:
  for ( auto k=0; k<u.size(); k++ ) { 
    for ( auto j=0; j<(*igbdy)[GBDY_NOSLIP].size() 
       && ub[k] != NULLPTR; j++ ) {
      (*ub[k])[j] = 0.0
    }
  }

  return TRUE;

} // end, impl_noslip


} // end, gupdateb namespace
