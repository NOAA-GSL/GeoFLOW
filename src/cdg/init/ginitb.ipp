
namespace ginitb {




//**********************************************************************************
//**********************************************************************************
// METHOD : impl_bystateinit
// DESC   : Initialize bdy using state initialization method 
//          specified by ptree
// ARGS   : ptree: main prop tree
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

  bret = gupdateb::impl_bystateinit(ptree, grid, tt, utmp, u, ub);

  return TRUE;

} // end, impl_bystateinit




} // end, ginitb namespace
