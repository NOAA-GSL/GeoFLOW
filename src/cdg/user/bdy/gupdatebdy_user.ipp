
namespace gupdatebdy {



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
GBOOL impl_simple_outflow(const PropertyTree &ptree, GGrid &grid, Time &time, State &utmp, State &u, State &ub)
{
  Time             tt = time;
  State            uu(u.size());
  GString          serr = "impl_simple_outflow: ";

  GTVector<GTVector<GSIZET>> *igbdy = &grid.igbdy_binned();

  // Set from State vector, u:
  for ( auto k=0; k<u.size(); k++ ) {
    for ( auto j=0; j<(*igbdy)[GBDY_OUTFLOW].size()
       && ub[k] != NULLPTR; j++ ) {
      (*ub[k])[j] = (*u[k])[(*igbdy)[GBDY_OUTFLOW][j]];
    }
  }

  return TRUE;

} // end, impl_simple_outflow


} // end, gupdatebdy namespace
