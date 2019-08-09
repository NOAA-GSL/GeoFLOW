
namespace gupdateb {


//**********************************************************************************
//**********************************************************************************
// METHOD : impl_mybdyupdate
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
GBOOL impl_mybdyupdate(const PropteryTree &ptree, GGrid &grid, Time &time, State &utmp, const State &u, State &ub)
{

#if 0
  Time             tt = t;
  State            uu(u.size());
  GString          serr = "impl_mybdyupdate: ";

  Time  tt = t;

  // Use tmp from back end, so that 'init' isn't 
  // disturbed:
  for ( auto j=0; j<u.size(); j++ ) {
    uu[j] = utmp[utmp.size()-1-j];
  }
  GInitSFactory::init(ptree, grid, tt, utmp, ub, uu);

  GTVector<GTVector<GSIZET>> *igbdy = &grid.igbdy();
  ...

  return TRUE;
#else
  return FALSE;
#endif

} // end, method impl_mybdyupdate


} // end, gupdateb namespace
