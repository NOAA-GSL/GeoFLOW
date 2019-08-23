
namespace ginitstate {




//**********************************************************************************
//**********************************************************************************
// METHOD : impl_mystateinit
// DESC   : Initialize bdy using state initialization method 
//          specified by ptree
// ARGS   : ptree: main prop tree
//          grid : grid
//          t    : time
//          utmp : tmp arrays
//          ub   : bdy vectors (one for each state element)
//          u    : current state, overwritten here
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_mystateinit(const PropertyTree &ptree, GGrid &grid, Time &time, State &utmp, State &ub, State &u)
{

  Time             tt = time;
  GString          serr = "impl_bystateinit: ";


  return FALSE;

} // end, impl_mystateinit




} // end, ginitstate namespace
