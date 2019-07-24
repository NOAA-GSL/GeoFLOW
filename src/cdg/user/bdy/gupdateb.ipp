
namespace gupdateb {




//**********************************************************************************
//**********************************************************************************
// METHOD : impl_bystateinit
// DESC   : Update bdy using state initialization method 
//          specified by ptree
// ARGS   : ptree: state prop tree
//          grid : grid
//          t    : time
//          utmp : tmp arrays
//          u    : current state, overwritten here
//          ub   : bdy vectors (one for each state element)
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_bystateinit(const PropteryTree &ptree, GGrid &grid, Time &time, State &utmp, State &u, State &ub)
{
  Time             tt = t;
  GString          serr = "impl_bystateinit: ";

  Time  tt = t;

  GInitSFactory::init(ptree, grid, tt, utmp, ub, u);


  // Set from State vector, u and others that we _can_ set:
  for ( auto k=0; k<u.size(); k++ ) { 
    for ( auto j=0; j<(*igbdy)[GBDY_DIRICHLET].size() 
       && ub[k] != NULLPTR; j++ ) {
      (*ub[k])[j] = (*u[k])[(*igbdy)[GBDY_DIRICHLET][j]];
    }
    for ( auto j=0; j<(*igbdy)[GBDY_INFLOW].size() 
       && ub[k] != NULLPTR; j++ ) {
      (*ub[k])[j] = (*u[k])[(*igbdy)[GBDY_INFLOW][j]];
    }
    for ( auto j=0; j<(*igbdy)[GBDY_0FLUX].size() 
       && ub[k] != NULLPTR; j++ ) {
      (*ub[k])[j] = 0.0;
    }
    for ( auto j=0; j<(*igbdy)[GBDY_NOSLIP].size() 
       && ub[k] != NULLPTR; j++ ) {
      (*ub[k])[j] = 0.0;
    }
  }


  

} // end, impl_bystateinit




} // end, gupdateb namespace
